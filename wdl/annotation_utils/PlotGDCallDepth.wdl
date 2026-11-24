version 1.0

import "../utils/Structs.wdl"
import "SummarizeGDCallDepth.wdl"

workflow PlotGDCallDepth {
    input {
        File depth_traces
        File call_windows
        File depth_bin_stats
        String prefix
        Array[File]? gd_vcfs
        Array[File]? gd_vcf_idxs

        Boolean nahr_check = true
        Int outside_bins = 5
        Float max_dev = 0.5
        File? manual_labels
        File? keep_samples
        Int raster_px = 2000
        File? sample_pop_assignments
        File? ped_file

        File? sd_pairs_bed
        File? segdup_bed
        File? genes_bed
        File? gene_info
        File? par_bed
        File? allosomes_list

        String plot_docker
        String? sv_pipeline_docker

        RuntimeAttr? runtime_attr_plot_region_depth
        RuntimeAttr? runtime_attr_apply_carrier_labels
    }

    call SummarizeGDCallDepth.PlotRegionDepth {
        input:
            traces = depth_traces,
            windows = call_windows,
            bin_stats = depth_bin_stats,
            sd_pairs_bed = sd_pairs_bed,
            segdup_bed = segdup_bed,
            genes_bed = genes_bed,
            gene_info = gene_info,
            manual_labels = manual_labels,
            keep_samples = keep_samples,
            nahr_check = nahr_check,
            outside_bins = outside_bins,
            max_dev = max_dev,
            raster_px = raster_px,
            prefix = prefix,
            docker = plot_docker,
            runtime_attr_override = runtime_attr_plot_region_depth
    }

    if (defined(gd_vcfs)) {
        scatter (i in range(length(select_first([gd_vcfs])))) {
            call ApplyCarrierLabelsToVcf {
                input:
                    vcf = select_first([gd_vcfs])[i],
                    vcf_idx = select_first([gd_vcf_idxs])[i],
                    carrier_labels = PlotRegionDepth.carrier_labels,
                    manual_labels = manual_labels,
                    sample_pop_assignments = sample_pop_assignments,
                    ped_file = ped_file,
                    par_bed = par_bed,
                    allosomes_list = allosomes_list,
                    prefix = "~{prefix}.relabeled.~{i}",
                    docker = select_first([sv_pipeline_docker]),
                    runtime_attr_override = runtime_attr_apply_carrier_labels
            }
        }
    }

    output {
        File region_depth_pdf = PlotRegionDepth.region_depth_pdf
        File region_counts = PlotRegionDepth.region_counts
        File? carrier_flags = PlotRegionDepth.carrier_flags
        File carrier_labels = PlotRegionDepth.carrier_labels
        Array[File]? relabeled_vcfs = ApplyCarrierLabelsToVcf.relabeled_vcf
        Array[File]? relabeled_vcf_idxs = ApplyCarrierLabelsToVcf.relabeled_vcf_idx
        Array[File]? relabel_summaries = ApplyCarrierLabelsToVcf.relabel_summary
    }
}

task ApplyCarrierLabelsToVcf {
    input {
        File vcf
        File vcf_idx
        File carrier_labels
        File? manual_labels
        File? sample_pop_assignments
        File? ped_file
        File? par_bed
        File? allosomes_list
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        cat > relabel.py <<'PYSCRIPT'
        import argparse
        import gzip

        import pysam

        p = argparse.ArgumentParser()
        p.add_argument("--vcf", required=True)
        p.add_argument("--labels", required=True)
        p.add_argument("--manual", default=None)
        p.add_argument("--out", required=True)
        p.add_argument("--summary", required=True)
        a = p.parse_args()

        # carriers whose final status is non_NAHR, keyed on (chrom, start, end, svtype) of the SD-pair region (= POS - 1, END, SVTYPE)
        reset = {}
        with gzip.open(a.labels, "rt") as f:
            cols = f.readline().rstrip("\n").split("\t")
            for line in f:
                r = dict(zip(cols, line.rstrip("\n").split("\t")))
                if r["final"] == "non_NAHR":
                    reset.setdefault((r["chrom"], int(r["start"]), int(r["end"]), r["svtype"]), set()).add(r["sample"])
        blacklist = set()
        if a.manual:
            with (gzip.open(a.manual, "rt") if a.manual.endswith(".gz") else open(a.manual)) as f:
                cols = f.readline().rstrip("\n").split("\t")
                for line in f:
                    r = dict(zip(cols, line.rstrip("\n").split("\t")))
                    if r["label"] == "blacklist":
                        blacklist.add(r["region"])

        vin = pysam.VariantFile(a.vcf)
        hdr = vin.header
        hdr.add_line('##FORMAT=<ID=NON_NAHR,Number=1,Type=Integer,Description="1: carrier relabelled as not NAHR between the '
                     'flanking SD pair (CNV beyond the SDs, by depth review); GT set to 0/0">')
        hdr.add_line('##FILTER=<ID=BLACKLIST,Description="SD-pair region blacklisted after depth review">')
        vout = pysam.VariantFile(a.out, "w", header=hdr)
        with open(a.summary, "w") as s:
            s.write("id\tregion\tsvtype\tn_reset\tn_reset_missing\tblacklist\n")
            for rec in vin:
                svtype = rec.info.get("SVTYPE")
                region = f"{rec.chrom}:{rec.pos - 1}-{rec.stop}"
                samples = reset.get((rec.chrom, rec.pos - 1, rec.stop, svtype), set())
                n, miss = 0, 0
                for sm in samples:
                    if sm not in rec.samples:
                        miss += 1
                        continue
                    rec.samples[sm]["GT"] = (0, 0)
                    rec.samples[sm]["NON_NAHR"] = 1
                    n += 1
                bl = region in blacklist
                if bl:
                    if "PASS" in rec.filter.keys():
                        rec.filter.clear()
                    rec.filter.add("BLACKLIST")
                if samples or bl:
                    s.write(f"{rec.id}\t{region}\t{svtype}\t{n}\t{miss}\t{bl}\n")
                vout.write(rec)
        vout.close()
        PYSCRIPT

        python3 relabel.py \
            --vcf ~{vcf} \
            --labels ~{carrier_labels} \
            ~{"--manual " + manual_labels} \
            --out relabeled.vcf.gz \
            --summary ~{prefix}.relabel_summary.tsv

        python3 /opt/sv-pipeline/05_annotation/scripts/compute_AFs.py relabeled.vcf.gz stdout \
            ~{"-p " + sample_pop_assignments} \
            ~{"-f " + ped_file} \
            ~{"--par " + par_bed} \
            ~{"--allosomes-list " + allosomes_list} \
            | bgzip -c > ~{prefix}.vcf.gz
        tabix -p vcf ~{prefix}.vcf.gz
    >>>

    output {
        File relabeled_vcf = "~{prefix}.vcf.gz"
        File relabeled_vcf_idx = "~{prefix}.vcf.gz.tbi"
        File relabel_summary = "~{prefix}.relabel_summary.tsv"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 8,
        disk_gb: 3 * ceil(size(vcf, "GB")) + 10,
        boot_disk_gb: 10,
        preemptible_tries: 2,
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
