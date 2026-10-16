version 1.0

import "../utils/Structs.wdl"

workflow CreateCNVDepthPlots {
    meta {
        description: [
            "This utility draws read-depth plots for the CNVs of a `LongReadCNVs` `genotyped_depth_vcf` that fall within a set of regions, so calls can be reviewed by eye. It writes one plot per non-reference sample per selected CNV, and packs them into one tarball with a folder per variant ID holding one PNG per non-reference sample, named by sample ID. A table lists every selected CNV with its plot window and non-reference samples.",
            "A record is selected when it overlaps any region in `regions_bed`, or lies wholly within one when `require_full_containment` is set, and passes the SV type and length filters. A sample is non-reference when its `FORMAT/RD_CN` is below its `FORMAT/ECN` for a DEL, above it for a DUP, or differs from it for a multiallelic CNV; genotypes are not consulted, since CNV records carry none. Each plot covers the call plus a window of `padding_fraction` of its length on either side. It shows every cohort sample's copy number estimate, 2 x depth / sample median coverage, from `merged_bincov` averaged to about 120 points whatever the window's length:the plotted sample in bold orange, the variant's other non-reference samples in light orange, and all remaining samples in grey with their median in black. Dashed lines mark the call. A sample in the VCF but missing from `merged_bincov` or `median_coverage` fails the workflow.",
            "Every non-reference sample is plotted by default, which can mean thousands of plots for a large region or a common CNV. Set `max_carriers_per_variant` to skip variants with more non-reference samples than that; skipped variants are still listed in the table, with `plotted` set to 0."
        ]
    }

    parameter_meta {
        genotyped_depth_vcf: "Genotyped CNV VCF from `LongReadCNVs`, carrying `FORMAT/RD_CN` and `FORMAT/ECN`."
        genotyped_depth_vcf_idx: "Index for `genotyped_depth_vcf`."
        regions_bed: "BED of regions whose CNVs are plotted; overlapping regions select a CNV once."
        merged_bincov: "Cohort binned read-depth matrix, the `merged_bincov` input to `LongReadCNVs`."
        merged_bincov_idx: "Index for `merged_bincov`."
        median_coverage: "Per-sample median coverage table, the `median_coverage` input to `LongReadCNVs`."
        svtypes: "SV types to plot, from `DEL`, `DUP` and `CNV`."
        min_svlen: "Minimum CNV length, in bp."
        max_svlen: "Maximum CNV length, in bp; unset keeps every length."
        require_full_containment: "Whether a CNV must lie wholly within a region, rather than overlap it, to be selected."
        padding_fraction: "Window added on each side of a CNV, as a share of its length."
        max_carriers_per_variant: "Maximum non-reference samples for a variant to be plotted; unset plots every variant."
        records_per_shard: "Number of variants plotted within a single shard."
        cnv_variants_tsv: "One row per selected CNV, with its `variant_id`, `contig`, `pos`, `end`, `svtype`, `svlen`, plot window `window_start` and `window_end`, `n_carriers` non-reference samples, `plotted` (1, or 0 when skipped by `max_carriers_per_variant`) and `carriers` as comma-separated `sample:RD_CN:ECN`."
        cnv_depth_plots_tar: "Tarball of `cnv_depth_plots/<variant_id>/<sample>.png`, one plot per non-reference sample of each plotted CNV."
    }

    input {
        File genotyped_depth_vcf
        File genotyped_depth_vcf_idx
        File regions_bed
        File merged_bincov
        File merged_bincov_idx
        File median_coverage
        String prefix

        Array[String] svtypes = ["DEL", "DUP", "CNV"]
        Int min_svlen = 0
        Int? max_svlen
        Boolean require_full_containment = false
        Float padding_fraction = 0.25
        Int? max_carriers_per_variant

        String utils_docker
        Int records_per_shard = 100

        RuntimeAttr? runtime_attr_select_cnvs_in_regions
        RuntimeAttr? runtime_attr_slice_depth_matrix
        RuntimeAttr? runtime_attr_create_depth_plots
        RuntimeAttr? runtime_attr_merge_plot_tarballs
    }

    call SelectCNVsInRegions {
        input:
            vcf = genotyped_depth_vcf,
            vcf_idx = genotyped_depth_vcf_idx,
            regions_bed = regions_bed,
            svtypes = svtypes,
            min_svlen = min_svlen,
            max_svlen = max_svlen,
            require_full_containment = require_full_containment,
            padding_fraction = padding_fraction,
            max_carriers_per_variant = max_carriers_per_variant,
            records_per_shard = records_per_shard,
            prefix = prefix,
            docker = utils_docker,
            runtime_attr_override = runtime_attr_select_cnvs_in_regions
    }

    call SliceDepthMatrix {
        input:
            merged_bincov = merged_bincov,
            merged_bincov_idx = merged_bincov_idx,
            windows_bed = SelectCNVsInRegions.windows_bed,
            prefix = prefix,
            docker = utils_docker,
            runtime_attr_override = runtime_attr_slice_depth_matrix
    }

    scatter (shard_i in range(length(SelectCNVsInRegions.variant_shards))) {
        call CreateDepthPlots {
            input:
                variant_shard = SelectCNVsInRegions.variant_shards[shard_i],
                rd_slice = SliceDepthMatrix.rd_slice,
                rd_slice_idx = SliceDepthMatrix.rd_slice_idx,
                median_coverage = median_coverage,
                prefix = "~{prefix}.shard_~{shard_i}",
                docker = utils_docker,
                runtime_attr_override = runtime_attr_create_depth_plots
        }
    }

    call MergePlotTarballs {
        input:
            plot_tars = CreateDepthPlots.plots_tar,
            prefix = prefix,
            docker = utils_docker,
            runtime_attr_override = runtime_attr_merge_plot_tarballs
    }

    output {
        File cnv_variants_tsv = SelectCNVsInRegions.variants_tsv
        File cnv_depth_plots_tar = MergePlotTarballs.merged_tar
    }
}

task SelectCNVsInRegions {
    input {
        File vcf
        File vcf_idx
        File regions_bed
        Array[String] svtypes
        Int min_svlen
        Int? max_svlen
        Boolean require_full_containment
        Float padding_fraction
        Int? max_carriers_per_variant
        Int records_per_shard
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        if [[ "~{vcf_idx}" != "~{vcf}.tbi" ]]; then
            ln -sf "~{vcf_idx}" "~{vcf}.tbi"
        fi

        mkdir shards

        python3 <<'CODE'
import csv
import pysam

SVTYPES = set("~{sep=',' svtypes}".split(','))
MIN_SVLEN = ~{min_svlen}
MAX_SVLEN = ~{default="None" max_svlen}
FULL_CONTAINMENT = ~{if require_full_containment then "True" else "False"}
PADDING_FRACTION = ~{padding_fraction}
MAX_CARRIERS = ~{default="None" max_carriers_per_variant}
PER_SHARD = ~{records_per_shard}
COLUMNS = ['variant_id', 'contig', 'pos', 'end', 'svtype', 'svlen', 'window_start', 'window_end', 'n_carriers', 'plotted',
           'carriers']
NONREF = {'DEL': lambda cn, ecn: cn < ecn, 'DUP': lambda cn, ecn: cn > ecn, 'CNV': lambda cn, ecn: cn != ecn}

vcf = pysam.VariantFile("~{vcf}")
indexed = set(vcf.index)
# Key by variant ID so a CNV hit by several regions is kept once
selected = {}
with open("~{regions_bed}") as handle:
    for line in handle:
        if not line.strip() or line.startswith(('#', 'track', 'browser')):
            continue
        contig, start, end = line.split('\t')[:3]
        start, end = int(start), int(end)
        if contig not in indexed:
            continue
        for rec in vcf.fetch(contig, start, end):
            svtype = rec.info['SVTYPE']
            svlen = rec.stop - rec.pos + 1
            if rec.id in selected or svtype not in SVTYPES or svlen < MIN_SVLEN:
                continue
            if MAX_SVLEN is not None and svlen > MAX_SVLEN:
                continue
            if FULL_CONTAINMENT and (rec.pos - 1 < start or rec.stop > end):
                continue
            carriers = []
            for sample, fmt in rec.samples.items():
                cn, ecn = fmt.get('RD_CN'), fmt.get('ECN')
                if cn is not None and ecn is not None and NONREF[svtype](cn, ecn):
                    carriers.append(f'{sample}:{cn}:{ecn}')
            if not carriers:
                continue
            padding = int(PADDING_FRACTION * svlen)
            plotted = int(MAX_CARRIERS is None or len(carriers) <= MAX_CARRIERS)
            selected[rec.id] = [rec.id, contig, rec.pos, rec.stop, svtype, svlen, max(0, rec.pos - 1 - padding),
                                rec.stop + padding, len(carriers), plotted, ','.join(carriers)]

rows = sorted(selected.values(), key=lambda row: (row[1], row[2], row[0]))
with open("~{prefix}.cnv_variants.tsv", 'w', newline='') as handle:
    writer = csv.writer(handle, delimiter='\t')
    writer.writerow(COLUMNS)
    writer.writerows(rows)

plotted_rows = [row for row in rows if row[9]]
with open('windows.bed', 'w') as handle:
    for row in plotted_rows:
        handle.write(f'{row[1]}\t{row[6]}\t{row[7]}\n')
for shard_i in range(0, len(plotted_rows), PER_SHARD):
    with open(f'shards/shard_{shard_i // PER_SHARD:06d}.tsv', 'w', newline='') as handle:
        writer = csv.writer(handle, delimiter='\t')
        writer.writerow(COLUMNS)
        writer.writerows(plotted_rows[shard_i:shard_i + PER_SHARD])

print(f'Selected {len(rows)} CNVs; plotting {len(plotted_rows)} with {sum(row[8] for row in plotted_rows)} plots')
CODE

        # Merge windows within 1 kb so no depth bin is fetched twice from neighbouring windows
        sort -k1,1 -k2,2n windows.bed | bedtools merge -d 1000 -i - > ~{prefix}.windows.bed
    >>>

    output {
        File variants_tsv = "~{prefix}.cnv_variants.tsv"
        File windows_bed = "~{prefix}.windows.bed"
        Array[File] variant_shards = glob("shards/shard_*.tsv")
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 4,
        disk_gb: ceil(size(vcf, "GB")) + 10,
        boot_disk_gb: 10,
        preemptible_tries: 1,
        max_retries: 0
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
        memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
        disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
        bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
        docker: docker
        preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
    }
}

task SliceDepthMatrix {
    input {
        File merged_bincov
        File merged_bincov_idx
        File windows_bed
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    parameter_meta {
        merged_bincov: { localization_optional: true }
        merged_bincov_idx: { localization_optional: true }
    }

    command <<<
        set -euo pipefail

        export GCS_OAUTH_TOKEN=$(gcloud auth application-default print-access-token)

        tabix -h -R ~{windows_bed} "~{merged_bincov}##idx##~{merged_bincov_idx}" | bgzip -c > ~{prefix}.windows.RD.txt.gz
        tabix -p bed ~{prefix}.windows.RD.txt.gz
    >>>

    output {
        File rd_slice = "~{prefix}.windows.RD.txt.gz"
        File rd_slice_idx = "~{prefix}.windows.RD.txt.gz.tbi"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 2,
        disk_gb: ceil(size(merged_bincov, "GB")) + 20,
        boot_disk_gb: 10,
        preemptible_tries: 1,
        max_retries: 0
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
        memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
        disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
        bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
        docker: docker
        preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
    }
}

task CreateDepthPlots {
    input {
        File variant_shard
        File rd_slice
        File rd_slice_idx
        File median_coverage
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        if [[ "~{rd_slice_idx}" != "~{rd_slice}.tbi" ]]; then
            ln -sf "~{rd_slice_idx}" "~{rd_slice}.tbi"
        fi

        python3 <<'CODE'
import csv
import os
import sys

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pysam

# Average every window to about this many points, so plots look alike whatever the CNV length
POINTS = 120
OUT_DIR = 'cnv_depth_plots'
ORANGE = '#eb6834'
GREY = '#52514e'
BLACK = '#0b0b0b'

tabix = pysam.TabixFile("~{rd_slice}")
samples = tabix.header[-1].lstrip('#').rstrip('\n').split('\t')[3:]
column = {sample: i for i, sample in enumerate(samples)}
with open("~{median_coverage}") as handle:
    median_samples, median_values = (line.split() for line in handle)
median_by_sample = dict(zip(median_samples, map(float, median_values)))
missing = [sample for sample in samples if sample not in median_by_sample]
if missing:
    sys.exit(f'{len(missing)} samples in merged_bincov are missing from median_coverage, e.g. {missing[0]}')
medians = np.array([median_by_sample[sample] for sample in samples], dtype=np.float32)

n_plots = 0
with open("~{variant_shard}", newline='') as handle:
    for row in csv.DictReader(handle, delimiter='\t'):
        variant_id, contig, svtype = row['variant_id'], row['contig'], row['svtype']
        pos, end, window_start, window_end = (int(row[key]) for key in ('pos', 'end', 'window_start', 'window_end'))
        carriers = [token.split(':') for token in row['carriers'].split(',')]
        missing = [sample for sample, _, _ in carriers if sample not in column]
        if missing:
            sys.exit(f'{variant_id}: non-reference sample {missing[0]} is missing from merged_bincov')

        # Average the bins of the window down to about POINTS points of copy number, 2 x depth / median
        bins = [line.split('\t') for line in tabix.fetch(contig, window_start, window_end)]
        starts = np.array([int(fields[1]) for fields in bins])
        depth = np.array([fields[3:] for fields in bins], dtype=np.float32) / medians
        step = max(1, len(bins) // POINTS)
        n_bins = len(bins) // step * step
        bin_width = int(bins[0][2]) - int(bins[0][1])
        x = (starts[:n_bins:step] + bin_width * step / 2) / 1e6
        cn = 2 * depth[:n_bins].reshape(-1, step, len(samples)).mean(axis=1)

        carrier_mask = np.zeros(len(samples), dtype=bool)
        carrier_mask[[column[sample] for sample, _, _ in carriers]] = True
        reference_median = np.median(cn[:, ~carrier_mask], axis=1)
        os.makedirs(f'{OUT_DIR}/{variant_id}', exist_ok=True)
        for sample, rd_cn, ecn in carriers:
            others = carrier_mask.copy()
            others[column[sample]] = False
            fig, ax = plt.subplots(figsize=(6.4, 3.2))
            ax.plot(x, cn[:, ~carrier_mask], color=GREY, lw=0.3, alpha=0.12)
            ax.plot(x, cn[:, others], color=ORANGE, lw=0.5, alpha=0.3)
            ax.plot(x, reference_median, color=BLACK, lw=1.6)
            ax.plot(x, cn[:, column[sample]], color=ORANGE, lw=1.6)
            for edge in (pos - 1, end):
                ax.axvline(edge / 1e6, color=BLACK, ls='--', lw=0.8)
            top = max(4, int(rd_cn) + 1)
            ax.set_ylim(0, top)
            ax.set_yticks(range(top + 1))
            ax.set_xlim(window_start / 1e6, window_end / 1e6)
            ax.ticklabel_format(axis='x', useOffset=False)
            ax.set_xlabel(f'{contig} (Mb)')
            ax.set_ylabel('Copy number')
            ax.set_title(f'{variant_id} {svtype} {contig}:{pos:,}-{end:,} ({(end - pos + 1) / 1000:,.1f} kb)\n'
                         f'{sample}: RD_CN {rd_cn}, ECN {ecn}; {len(carriers)} non-reference samples', fontsize=8)
            fig.savefig(f'{OUT_DIR}/{variant_id}/{sample}.png', dpi=150, bbox_inches='tight')
            plt.close(fig)
            n_plots += 1

print(f'Wrote {n_plots} plots')
CODE

        tar -czf ~{prefix}.cnv_depth_plots.tar.gz cnv_depth_plots
    >>>

    output {
        File plots_tar = "~{prefix}.cnv_depth_plots.tar.gz"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 4,
        disk_gb: ceil(size(rd_slice, "GB")) * 2 + 10,
        boot_disk_gb: 10,
        preemptible_tries: 1,
        max_retries: 0
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
        memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
        disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
        bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
        docker: docker
        preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
    }
}

task MergePlotTarballs {
    input {
        Array[File] plot_tars
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        mkdir cnv_depth_plots
        for TAR in ~{sep=' ' plot_tars}; do
            tar -xzf "$TAR"
        done

        tar -czf ~{prefix}.cnv_depth_plots.tar.gz cnv_depth_plots
    >>>

    output {
        File merged_tar = "~{prefix}.cnv_depth_plots.tar.gz"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 2,
        disk_gb: ceil(size(plot_tars, "GB")) * 3 + 10,
        boot_disk_gb: 10,
        preemptible_tries: 1,
        max_retries: 0
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
        memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
        disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
        bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
        docker: docker
        preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
    }
}
