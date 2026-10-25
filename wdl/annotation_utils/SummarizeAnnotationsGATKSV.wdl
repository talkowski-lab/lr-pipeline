# Derived from wdl/QcAnnotations.wdl of broadinstitute/gatk-sv (branch kj_project_gnomad_lr)
# https://github.com/broadinstitute/gatk-sv/blob/kj_project_gnomad_lr/wdl/QcAnnotations.wdl
# Derived from wdl/CollectQcVcfWide.wdl of broadinstitute/gatk-sv (branch kj_project_gnomad_lr)
# https://github.com/broadinstitute/gatk-sv/blob/main/wdl/CollectQcVcfWide.wdl
# Derived from wdl/CollectQcPerSample.wdl of broadinstitute/gatk-sv (branch kj_project_gnomad_lr)
# https://github.com/broadinstitute/gatk-sv/blob/main/wdl/CollectQcPerSample.wdl

version 1.0

import "../utils/Helpers.wdl"
import "../utils/Structs.wdl"

workflow SummarizeAnnotationsGATKSV {
    meta {
        description: [
            "This utility draws the GATK-SV VCF QC plots for one or more VCFs: site-level count, size, frequency and Hardy-Weinberg distributions by variant class, and, when `run_per_sample` is enabled, per-genome variant counts and per-family inheritance and de novo rates from the PED. It covers the same ground as `SummarizeAnnotations`, but uses the GATK-SV QC scripts in place of the summary tables.",
            "Each VCF is sharded by record count. Every shard has its multiallelic records split and its alleles converted to symbolic classes (SNV, INS_SHORT/INS_SV, DEL_SHORT/DEL_SV, DUP_SHORT/DUP_SV at a 50 bp cutoff, and TR or VNTR for TRV alleles by shortest motif) before its site statistics are collected. Benchmarking against external callsets is not ported."
        ]
    }

    parameter_meta {
        vcfs: "VCFs to summarize."
        vcf_idxs: "Indexes for `vcfs`."
        ped: "Cohort pedigree, used for the per-family plots."
        ref_fa: "From references."
        ref_fai: "From references."
        records_per_shard: "Number of variants to keep within a single shard."
        samples_per_shard: "Number of samples per shard when collecting per-sample variant lists."
        create_variant_attributes: "Whether to derive `allele_type` and `allele_length` from REF and ALT, for VCFs that lack them."
        run_per_sample: "Whether to additionally draw the per-sample and per-family plots."
        subset_sample_count: "Number of randomly chosen samples to restrict the per-sample and per-family plots to."
        subset_vcf_string: "`bcftools view` arguments used to pre-subset each shard."
        random_seed: "Seed for the random sample subset."
        max_gq: "Maximum GQ value on the per-sample and per-family GQ axes."
        downsample_qc_per_sample: "Number of samples to downsample the per-sample plots to."
        vcf_qc_tarball: "Tarball of the QC plots, the site statistics BED with its index, the samples analyzed and, when `run_per_sample` is enabled, the cleaned PED of trios."
        vcf2bed_bed: "Every record of the preprocessed shards with its called samples and INFO fields, as a gzipped BED."
    }

    input {
        Array[File] vcfs
        Array[File] vcf_idxs
        String prefix

        File ped
        File ref_fa
        File ref_fai

        Int records_per_shard
        Int samples_per_shard
        Boolean create_variant_attributes = false
        Boolean run_per_sample = true
        Int? subset_sample_count
        String subset_vcf_string = ""
        Int random_seed = 0
        Int max_gq = 99
        Int? downsample_qc_per_sample

        String utils_docker

        RuntimeAttr? runtime_attr_shard
        RuntimeAttr? runtime_attr_preprocess
        RuntimeAttr? runtime_attr_merge_stats
        RuntimeAttr? runtime_attr_merge_vcf2bed
        RuntimeAttr? runtime_attr_plot_vcf_wide
        RuntimeAttr? runtime_attr_subset_samples
        RuntimeAttr? runtime_attr_split_samples
        RuntimeAttr? runtime_attr_collect_vids
        RuntimeAttr? runtime_attr_merge_vids
        RuntimeAttr? runtime_attr_plot_per_sample
        RuntimeAttr? runtime_attr_plot_per_family
        RuntimeAttr? runtime_attr_sanitize
    }

    scatter (i in range(length(vcfs))) {
        call Helpers.ShardVcfByRecords {
            input:
                vcf = vcfs[i],
                vcf_idx = vcf_idxs[i],
                records_per_shard = records_per_shard,
                prefix = "~{prefix}.~{i}",
                docker = utils_docker,
                runtime_attr_override = runtime_attr_shard
        }

    }

    # Preprocess every shard in one flat scatter, since Cromwell can't read unset optional inputs inside a nested scatter
    Array[File] shards = flatten(ShardVcfByRecords.shards)

    scatter (j in range(length(shards))) {
        call PreprocessCollectAndConvert {
            input:
                vcf = shards[j],
                subset_vcf_string = subset_vcf_string,
                ref_fa = ref_fa,
                ref_fai = ref_fai,
                create_variant_attributes = create_variant_attributes,
                prefix = "~{prefix}.shard_~{j}",
                docker = utils_docker,
                runtime_attr_override = runtime_attr_preprocess
        }
    }

    File samples_list = PreprocessCollectAndConvert.samples_list[0]

    call ConcatBeds as MergeVcfWideStats {
        input:
            beds = PreprocessCollectAndConvert.vcf_stats,
            index_output = true,
            prefix = "~{prefix}.VCF_sites.stats",
            docker = utils_docker,
            runtime_attr_override = runtime_attr_merge_stats
    }

    call ConcatBeds as MergeVcf2Bed {
        input:
            beds = PreprocessCollectAndConvert.vcf2bed_out,
            index_output = false,
            prefix = "~{prefix}.vcf2bed",
            docker = utils_docker,
            runtime_attr_override = runtime_attr_merge_vcf2bed
    }

    call PlotQcVcfWide {
        input:
            vcf_stats = MergeVcfWideStats.merged_bed,
            samples_list = samples_list,
            prefix = prefix,
            docker = utils_docker,
            runtime_attr_override = runtime_attr_plot_vcf_wide
    }

    if (run_per_sample) {
        if (defined(subset_sample_count)) {
            call SubsetSamplesList {
                input:
                    samples_list = samples_list,
                    subset_count = select_first([subset_sample_count]),
                    random_seed = random_seed,
                    prefix = prefix,
                    docker = utils_docker,
                    runtime_attr_override = runtime_attr_subset_samples
            }
        }

        File per_sample_samples_list = select_first([SubsetSamplesList.subset_samples_list, samples_list])

        call SplitSamplesList {
            input:
                samples_list = per_sample_samples_list,
                samples_per_shard = samples_per_shard,
                prefix = "~{prefix}.list_shard.",
                docker = utils_docker,
                runtime_attr_override = runtime_attr_split_samples
        }

        scatter (vcf_and_samples in cross(vcfs, SplitSamplesList.shards)) {
            call CollectVidsPerSample {
                input:
                    vcf = vcf_and_samples.left,
                    samples_list = vcf_and_samples.right,
                    prefix = prefix,
                    docker = utils_docker,
                    runtime_attr_override = runtime_attr_collect_vids
            }
        }

        call MergeShardedPerSampleVidLists {
            input:
                tarballs = CollectVidsPerSample.vid_lists_tarball,
                samples_list = per_sample_samples_list,
                prefix = prefix,
                docker = utils_docker,
                runtime_attr_override = runtime_attr_merge_vids
        }

        call PlotQcPerSample {
            input:
                vcf_stats = MergeVcfWideStats.merged_bed,
                samples_list = per_sample_samples_list,
                per_sample_tarball = MergeShardedPerSampleVidLists.merged_tarball,
                max_gq = max_gq,
                downsample_qc_per_sample = downsample_qc_per_sample,
                prefix = prefix,
                docker = utils_docker,
                runtime_attr_override = runtime_attr_plot_per_sample
        }

        call PlotQcPerFamily {
            input:
                vcf_stats = MergeVcfWideStats.merged_bed,
                samples_list = per_sample_samples_list,
                ped = ped,
                per_sample_tarball = MergeShardedPerSampleVidLists.merged_tarball,
                max_gq = max_gq,
                prefix = prefix,
                docker = utils_docker,
                runtime_attr_override = runtime_attr_plot_per_family
        }
    }

    call SanitizeOutputs {
        input:
            samples_list = samples_list,
            vcf_stats = MergeVcfWideStats.merged_bed,
            vcf_stats_idx = MergeVcfWideStats.merged_bed_idx,
            plot_qc_vcfwide_tarball = PlotQcVcfWide.plots_tarball,
            plot_qc_per_sample_tarball = PlotQcPerSample.per_sample_plots_tarball,
            plot_qc_per_family_tarball = PlotQcPerFamily.per_family_plots_tarball,
            cleaned_fam_file = PlotQcPerFamily.cleaned_fam_file,
            prefix = prefix,
            docker = utils_docker,
            runtime_attr_override = runtime_attr_sanitize
    }

    output {
        File vcf_qc_tarball = SanitizeOutputs.vcf_qc_tarball
        File vcf2bed_bed = MergeVcf2Bed.merged_bed
    }
}

task PreprocessCollectAndConvert {
    input {
        File vcf
        String subset_vcf_string
        File ref_fa
        File ref_fai
        Boolean create_variant_attributes
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    Int sort_max_mem_gb = floor(0.9 * select_first([runtime_attr.mem_gb, default_attr.mem_gb]))

    command <<<
        set -euo pipefail

        # Optionally pre-filter records on the shard
        if [ -n "~{subset_vcf_string}" ]; then
            bcftools view --no-version ~{subset_vcf_string} -Oz -o filtered_input.vcf.gz ~{vcf}
        else
            ln -s ~{vcf} filtered_input.vcf.gz
        fi

        bcftools norm \
            -m-any \
            --check-ref w \
            --fasta-ref ~{ref_fa} \
            -Oz -o ~{prefix}.split.vcf.gz \
            filtered_input.vcf.gz

        rm -f filtered_input.vcf.gz

        # Convert to symbolic alleles, optionally annotating allele_length / allele_type
        python3 <<'PYCODE'
import pysam

VCF_IN = "~{prefix}.split.vcf.gz"
VCF_OUT = "~{prefix}.unsorted.vcf.gz"
CREATE_VARIANT_ATTRIBUTES = "~{create_variant_attributes}".lower() == "true"

group_carrier_counts = []
vcf_in = pysam.VariantFile(VCF_IN, "r")
current_locus = None
current_carrier_set = set()
for rec in vcf_in:
    if CREATE_VARIANT_ATTRIBUTES:
        ref_len = len(rec.ref)
        first_alt = rec.alts[0] if rec.alts else ""
        alt_len = len(first_alt)
        if ref_len == 1 and alt_len == 1:
            at = "SNV"
        elif alt_len > ref_len:
            at = "INS"
        else:
            at = "DEL"
    else:
        at = rec.info["allele_type"].upper() if "allele_type" in rec.info else ""

    if at == "TRV":
        key = (rec.chrom, rec.id)
        if key != current_locus:
            if current_locus is not None:
                group_carrier_counts.append(len(current_carrier_set))
            current_locus = key
            current_carrier_set = set()
        for s_name in rec.samples:
            gt = rec.samples[s_name].get("GT", (None,))
            if gt and any(a is not None and a > 0 for a in gt):
                current_carrier_set.add(s_name)
    elif current_locus is not None:
        group_carrier_counts.append(len(current_carrier_set))
        current_locus = None
        current_carrier_set = set()

if current_locus is not None:
    group_carrier_counts.append(len(current_carrier_set))
vcf_in.close()

vcf_in = pysam.VariantFile(VCF_IN, "r")

if CREATE_VARIANT_ATTRIBUTES:
    if "allele_length" not in vcf_in.header.info:
        vcf_in.header.info.add("allele_length", 1, "Integer", "Allele length (alt_len - ref_len)")
    if "allele_type" not in vcf_in.header.info:
        vcf_in.header.info.add("allele_type", 1, "String", "Allele type (snv/ins/del)")

if "SVTYPE" not in vcf_in.header.info:
    vcf_in.header.info.add("SVTYPE", 1, "String", "Type of variant.")
if "SVLEN" not in vcf_in.header.info:
    vcf_in.header.info.add("SVLEN", 1, "Integer", "Length of variant.")
if "END" not in vcf_in.header.info:
    vcf_in.header.info.add("END", 1, "Integer", "End position of variant.")
if "ORIG_ALT" not in vcf_in.header.info:
    vcf_in.header.info.add("ORIG_ALT", 1, "String", "Original ALT nucleotide for SNVs before symbolic conversion.")
if "TRV_EXPANSION_RATIO" not in vcf_in.header.info:
    vcf_in.header.info.add("TRV_EXPANSION_RATIO", 1, "Float", "Proportion of non-ref alleles called as expansions (alt allele length > ref allele length) among non-neutral alleles for TRV variants.")
if "TR_ALLELE_CLASS" not in vcf_in.header.info:
    vcf_in.header.info.add("TR_ALLELE_CLASS", 1, "String", "Sub-classification of TRV biallelic allele: TR_INS, TR_DEL, or TR_SNV.")
if "TRV_LOCUS_CARRIERS" not in vcf_in.header.info:
    vcf_in.header.info.add("TRV_LOCUS_CARRIERS", 1, "Integer", "Number of unique carriers across all biallelic alleles at this TRV locus.")
vcf_out = pysam.VariantFile(VCF_OUT, "w", header=vcf_in.header)

prev_trv_key = None
trv_counter = 0

prev_nontrv_key = None
nontrv_counter = 1

group_idx = -1
current_group_carriers = 0

for rec in vcf_in:
    if CREATE_VARIANT_ATTRIBUTES:
        ref_len = len(rec.ref)
        first_alt = rec.alts[0] if rec.alts else ""
        alt_len = len(first_alt)
        if ref_len == 1 and alt_len == 1:
            inferred_type = "snv"
        elif alt_len > ref_len:
            inferred_type = "ins"
        else:
            inferred_type = "del"
        rec.info["allele_length"] = alt_len - ref_len
        rec.info["allele_type"] = inferred_type

    allele_type = rec.info["allele_type"].upper()

    if allele_type == "TRV":
        current_key = (rec.chrom, rec.id)
        if current_key != prev_trv_key:
            group_idx += 1
            current_group_carriers = group_carrier_counts[group_idx]
            trv_counter = 1
            prev_trv_key = current_key
        else:
            trv_counter += 1
        rec.id = f"{rec.id}_{trv_counter}"

        motifs = rec.info.get("MOTIFS", None)
        if motifs:
            motif_list = list(motifs) if isinstance(motifs, (list, tuple)) else [str(motifs)]
            min_motif_len = min(len(str(m)) for m in motif_list)
        else:
            min_motif_len = 0
        svtype = "TR" if min_motif_len <= 6 else "VNTR"

        ref_seq = rec.ref
        alt_seq = rec.alts[0] if rec.alts else ""
        if len(alt_seq) > len(ref_seq):
            allele_class = "TR_INS"
            allele_length = len(alt_seq) - len(ref_seq)
        elif len(ref_seq) > len(alt_seq):
            allele_class = "TR_DEL"
            allele_length = len(ref_seq) - len(alt_seq)
        else:
            allele_class = "TR_SNV"
            allele_length = max(len(ref_seq), 1)

        rec.info["SVTYPE"] = svtype
        rec.info["SVLEN"] = allele_length
        rec.info["TR_ALLELE_CLASS"] = allele_class
        rec.info["TRV_LOCUS_CARRIERS"] = current_group_carriers
        rec.stop = rec.pos + len(rec.ref)

        rec.alts = (f"<{allele_type}>",)

        for s_name in rec.samples:
            try:
                rec.samples[s_name].phased = False
            except (AttributeError, TypeError):
                pass

        vcf_out.write(rec)
        continue

    prev_trv_key = None

    nontrv_key = (rec.chrom, rec.pos, rec.id)
    if nontrv_key == prev_nontrv_key:
        nontrv_counter += 1
        rec.id = f"{rec.id}_{nontrv_counter}"
    else:
        prev_nontrv_key = nontrv_key
        nontrv_counter = 1

    allele_length = abs(rec.info["allele_length"]) if "allele_length" in rec.info else len(rec.ref)

    if "DEL" in allele_type:
        svtype = "DEL_SHORT" if allele_length < 50 else "DEL_SV"
    elif "INS" in allele_type:
        svtype = "INS_SHORT" if allele_length < 50 else "INS_SV"
    elif "DUP" in allele_type:
        svtype = "DUP_SHORT" if allele_length < 50 else "DUP_SV"
    else:
        svtype = allele_type

    if allele_type == "SNV" and len(rec.alts[0]) == 1:
        rec.info["ORIG_ALT"] = rec.alts[0]
    rec.info["SVTYPE"] = svtype
    rec.info["SVLEN"] = allele_length
    rec.stop = rec.pos + len(rec.ref)

    rec.alts = (f"<{allele_type}>",)

    for s_name in rec.samples:
        try:
            rec.samples[s_name].phased = False
        except (AttributeError, TypeError):
            pass

    vcf_out.write(rec)

vcf_out.close()
vcf_in.close()
PYCODE

        rm -f ~{prefix}.split.vcf.gz

        bcftools sort \
            --max-mem ~{sort_max_mem_gb}G \
            -T . \
            -Oz -o ~{prefix}.sorted.vcf.gz \
            ~{prefix}.unsorted.vcf.gz

        rm -f ~{prefix}.unsorted.vcf.gz

        tabix ~{prefix}.sorted.vcf.gz

        # Keep collectQC's scratch on the sized local disk rather than the boot disk's /tmp
        mkdir -p collectQC_tmp
        TMPDIR="$(pwd)/collectQC_tmp" bash /opt/scripts/qc/vcf_qc/collectQC.vcf_wide.sh \
            ~{prefix}.sorted.vcf.gz \
            /opt/scripts/qc/vcf_qc/SV_colors.txt \
            collectQC_vcfwide_output/

        cp collectQC_vcfwide_output/data/VCF_sites.stats.bed.gz ~{prefix}.VCF_sites.stats.bed.gz
        cp collectQC_vcfwide_output/analysis_samples.list ~{prefix}.analysis_samples.list
        bgzip -c collectQC_vcfwide_output/vcf2bed_raw.bed > ~{prefix}.vcf2bed.bed.gz

        rm -f ~{prefix}.sorted.vcf.gz ~{prefix}.sorted.vcf.gz.tbi
        rm -rf collectQC_vcfwide_output/ collectQC_tmp/
    >>>

    output {
        File vcf_stats = "~{prefix}.VCF_sites.stats.bed.gz"
        File samples_list = "~{prefix}.analysis_samples.list"
        File vcf2bed_out = "~{prefix}.vcf2bed.bed.gz"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 16,
        # Split TR alleles each carry every INFO field and collectQC keeps several uncompressed copies, so a 37 MB shard can exceed 25 GB
        disk_gb: 50 * ceil(size(vcf, "GB")) + 100,
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

task ConcatBeds {
    input {
        Array[File] beds
        Boolean index_output
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        # Keep the first BED's header, then merge and sort every BED's records
        (zcat ~{beds[0]} || true) | head -n1 > header.txt

        while read -r BED; do
            zcat "$BED"
        done < ~{write_lines(beds)} \
            | (grep -Ev "^#" || printf "") \
            | sort -Vk1,1 -k2,2n -k3,3n \
            | cat header.txt - \
            | bgzip -c \
            > ~{prefix}.bed.gz

        if ~{index_output}; then
            tabix -p bed ~{prefix}.bed.gz
        else
            touch ~{prefix}.bed.gz.tbi
        fi
    >>>

    output {
        File merged_bed = "~{prefix}.bed.gz"
        File merged_bed_idx = "~{prefix}.bed.gz.tbi"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 2,
        disk_gb: 7 * ceil(size(beds, "GB")) + 10,
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

task PlotQcVcfWide {
    input {
        File vcf_stats
        File samples_list
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        Rscript /opt/scripts/qc/vcf_qc/plot_sv_vcf_distribs.R \
            --skip_supporting \
            -N $(sort -u ~{samples_list} | wc -l) \
            -S /opt/scripts/qc/vcf_qc/SV_colors.txt \
            ~{vcf_stats} \
            plotQC_vcfwide_output/

        rm -rf plotQC_vcfwide_output/supporting_plots/

        tar -czf ~{prefix}.plotQC_vcfwide_output.tar.gz plotQC_vcfwide_output
    >>>

    output {
        File plots_tarball = "~{prefix}.plotQC_vcfwide_output.tar.gz"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 8,
        disk_gb: 5 * ceil(size(vcf_stats, "GB")) + 10,
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

task SubsetSamplesList {
    input {
        File samples_list
        Int subset_count
        Int random_seed
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        sort -u ~{samples_list} > all_samples.uniq.list
        n_total=$(wc -l < all_samples.uniq.list)

        if [ "$n_total" -le ~{subset_count} ]; then
            cp all_samples.uniq.list ~{prefix}.subset_samples.list
        else
            sort -R --random-source <(yes ~{random_seed}) all_samples.uniq.list \
                | head -n ~{subset_count} \
                | sort \
                > ~{prefix}.subset_samples.list
        fi
    >>>

    output {
        File subset_samples_list = "~{prefix}.subset_samples.list"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 1,
        disk_gb: 10,
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

task SplitSamplesList {
    input {
        File samples_list
        Int samples_per_shard
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        N_LINES=$(wc -l < ~{samples_list})
        N_CHUNKS=$((N_LINES / ~{samples_per_shard}))
        if [ "$N_CHUNKS" -eq 0 ]; then N_CHUNKS=1; fi
        N_DIGITS=${#N_CHUNKS}

        split -d -a "$N_DIGITS" -n l/"$N_CHUNKS" \
            --numeric-suffixes=$(printf "%0${N_DIGITS}d" 1) \
            ~{samples_list} \
            ~{prefix}
    >>>

    output {
        Array[File] shards = glob("~{prefix}*")
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 2,
        disk_gb: 2 * ceil(size(samples_list, "GB")) + 10,
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

task CollectVidsPerSample {
    input {
        File vcf
        File samples_list
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    String outdirprefix = "~{prefix}_perSample_VIDs"

    command <<<
        set -euo pipefail

        mkdir -p ~{outdirprefix}

        # Write each carrier's variant IDs, collapsed genotypes and GQs to a per-sample list
        gq_fmt=$(bcftools view -h ~{vcf} | grep -q '##FORMAT=<ID=GQ' && echo '%GQ' || echo 'NA')
        bcftools view -S ~{samples_list} ~{vcf} \
            | bcftools +fill-tags -- -t AC \
            | bcftools view -i 'AC>0' \
            | bcftools query -f "[%SAMPLE\t%ID\t%GT\t${gq_fmt}\n]" \
            | awk '{OFS="\t"; gt=$3; gsub(/\|/,"/",gt); gq=$4; if(gq=="."){gq=99}; split(gt,a,"/"); x_miss=(a[1]=="."); y_miss=(a[2]=="."); if(x_miss&&y_miss){gt="./."}else if(x_miss){if(a[2]+0>0) gt="0/1"; else gt="0/0"}else if(y_miss){if(a[1]+0>0) gt="0/1"; else gt="0/0"}else{x=a[1]+0;y=a[2]+0; if(x==0&&y==0) gt="0/0"; else if(x==0||y==0) gt="0/1"; else gt="1/1"}; print $1,$2,gt,gq}' \
            | awk -v outprefix="~{outdirprefix}" '$3 != "0/0" && $3 != "./." {OFS="\t"; print $2, $3, $4 >> outprefix"/"$1".VIDs_genotypes.txt" }'

        for FILE in ~{outdirprefix}/*.VIDs_genotypes.txt; do
            gzip -f "$FILE"
        done

        cd ~{outdirprefix}
        tar -czf ../~{outdirprefix}.tar.gz *.VIDs_genotypes.txt.gz
    >>>

    output {
        File vid_lists_tarball = "~{outdirprefix}.tar.gz"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 4,
        disk_gb: ceil(size([vcf, samples_list], "GB")) + 10,
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

task MergeShardedPerSampleVidLists {
    input {
        Array[File] tarballs
        File samples_list
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        mkdir "~{prefix}_perSample_VID_lists"
        mkdir shards

        while read -r i tarball_path; do
            mkdir "shards/shard_$i"
            tar -xzf "$tarball_path" --directory "shards/shard_$i"/
        done < <(awk -v OFS="\t" '{ print NR, $1 }' ~{write_lines(tarballs)})

        N_VID_FILES=$(find shards/ -name "*.VIDs_genotypes.txt.gz" | wc -l)
        if [ "$N_VID_FILES" -eq 0 ]; then
            echo "ERROR: No VID files found in extracted tarballs" >&2
            exit 1
        fi

        OUTDIR="~{prefix}_perSample_VID_lists"
        export OUTDIR

        merge_sample() {
            local sample="$1"
            find shards/ -name "${sample}.VIDs_genotypes.txt.gz" \
                | xargs -I {} zcat {} \
                | sort -Vk1,1 -k2,2n -k3,3n \
                | gzip -c \
                > "${OUTDIR}/${sample}.VIDs_genotypes.txt.gz"
        }
        export -f merge_sample

        xargs -P "$(nproc)" -n 1 bash -c 'set -euo pipefail; merge_sample "$1"' _ < ~{samples_list}

        N_NONEMPTY=$(find "${OUTDIR}" -name "*.VIDs_genotypes.txt.gz" -size +20c | wc -l)
        if [ "$N_NONEMPTY" -eq 0 ]; then
            echo "ERROR: All merged VID files are empty" >&2
            exit 1
        fi

        tar -czf "~{prefix}_perSample_VID_lists.tar.gz" "~{prefix}_perSample_VID_lists"
    >>>

    output {
        File merged_tarball = "~{prefix}_perSample_VID_lists.tar.gz"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 2,
        mem_gb: 4,
        disk_gb: 5 * ceil(size(tarballs, "GB")) + 10,
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

task PlotQcPerSample {
    input {
        File vcf_stats
        File samples_list
        File per_sample_tarball
        Int max_gq
        Int? downsample_qc_per_sample
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        mkdir ~{prefix}_perSample/ tmp_untar/
        tar -xzf ~{per_sample_tarball} --directory tmp_untar/
        find tmp_untar/ -name "*.VIDs_genotypes.txt.gz" -exec mv {} ~{prefix}_perSample/ \;

        Rscript /opt/scripts/qc/vcf_qc/plot_sv_perSample_distribs.R \
            --skip_supporting \
            -S /opt/scripts/qc/vcf_qc/SV_colors.txt \
            ~{vcf_stats} \
            ~{samples_list} \
            ~{prefix}_perSample/ \
            ~{prefix}_perSample_plots/ \
            --maxgq ~{max_gq} \
            ~{"--downsample " + downsample_qc_per_sample}

        rm -rf ~{prefix}_perSample_plots/supporting_plots/

        tar -czf ~{prefix}.plotQC_perSample.tar.gz ~{prefix}_perSample_plots
    >>>

    output {
        File per_sample_plots_tarball = "~{prefix}.plotQC_perSample.tar.gz"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 4,
        disk_gb: 5 * ceil(size([vcf_stats, per_sample_tarball], "GB")) + 10,
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

task PlotQcPerFamily {
    input {
        File vcf_stats
        File samples_list
        File ped
        File per_sample_tarball
        Int max_gq
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        bash /opt/scripts/qc/vcf_qc/cleanFamFile.sh \
            ~{samples_list} \
            ~{ped} \
            ~{prefix}.cleaned.fam

        n_fams=$(grep -Evc "^#" ~{prefix}.cleaned.fam || true)
        echo "DETECTED $n_fams FAMILIES"
        if [ "$n_fams" -gt 0 ]; then
            mkdir ~{prefix}_perSample/ tmp_untar/
            tar -xzf ~{per_sample_tarball} --directory tmp_untar/
            find tmp_untar/ -name "*.VIDs_genotypes.txt.gz" -exec mv {} ~{prefix}_perSample/ \;

            Rscript /opt/scripts/qc/vcf_qc/analyze_fams.R \
                --skip_supporting \
                -S /opt/scripts/qc/vcf_qc/SV_colors.txt \
                ~{vcf_stats} \
                ~{prefix}.cleaned.fam \
                ~{prefix}_perSample/ \
                ~{prefix}_perFamily_plots/ \
                --maxgq ~{max_gq}
        else
            mkdir ~{prefix}_perFamily_plots/
        fi

        rm -rf ~{prefix}_perFamily_plots/supporting_plots/

        tar -czf ~{prefix}.plotQC_perFamily.tar.gz ~{prefix}_perFamily_plots
    >>>

    output {
        File per_family_plots_tarball = "~{prefix}.plotQC_perFamily.tar.gz"
        File cleaned_fam_file = "~{prefix}.cleaned.fam"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 4,
        disk_gb: 5 * ceil(size([vcf_stats, per_sample_tarball], "GB")) + 10,
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

task SanitizeOutputs {
    input {
        File samples_list
        File vcf_stats
        File vcf_stats_idx
        File plot_qc_vcfwide_tarball
        File? plot_qc_per_sample_tarball
        File? plot_qc_per_family_tarball
        File? cleaned_fam_file
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    String outdir = "~{prefix}_SV_VCF_QC_output"

    command <<<
        set -euo pipefail

        mkdir -p ~{outdir}/data ~{outdir}/plots

        cp ~{vcf_stats} ~{outdir}/data/~{prefix}.VCF_sites.stats.bed.gz
        cp ~{vcf_stats_idx} ~{outdir}/data/~{prefix}.VCF_sites.stats.bed.gz.tbi
        tar -xzf ~{plot_qc_vcfwide_tarball}
        cp plotQC_vcfwide_output/main_plots/* ~{outdir}/plots/

        if ~{defined(plot_qc_per_sample_tarball)}; then
            tar -xzf ~{plot_qc_per_sample_tarball}
            cp ~{prefix}_perSample_plots/main_plots/* ~{outdir}/plots/
        fi

        if ~{defined(plot_qc_per_family_tarball)}; then
            tar -xzf ~{plot_qc_per_family_tarball}
            cp ~{prefix}_perFamily_plots/main_plots/* ~{outdir}/plots/ || true
        fi

        if ~{defined(cleaned_fam_file)}; then
            cp ~{cleaned_fam_file} ~{outdir}/data/~{prefix}.cleaned_trios.fam
        fi

        cp ~{samples_list} ~{outdir}/data/~{prefix}.samples_analyzed.list

        tar -czf ~{outdir}.tar.gz ~{outdir}
    >>>

    output {
        File vcf_qc_tarball = "~{outdir}.tar.gz"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 2,
        disk_gb: 5 * ceil(size([plot_qc_vcfwide_tarball, plot_qc_per_sample_tarball, plot_qc_per_family_tarball, vcf_stats], "GB")) + 10,
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
