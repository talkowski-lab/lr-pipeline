version 1.0

import "../wdl/utils/Structs.wdl"

# Tandem-repeat variant (TRV) analysis of per-contig cohort VCFs.
#
# `vcfs` / `vcf_idxs` / `contigs` are parallel per-contig arrays. All analysis logic lives in analyze_TRV_vcf.py
# (passed in as `analysis_script`; standard library only). Per contig (scattered):
#   1. ExtractTRV: TRV records (INFO/allele_type == "trv") to a VCF, the sample list, and a bcftools query table of
#      site fields (CHROM, POS, ID, REF, ALT, FILTER, TRID, MOTIFS, AC, AF, AN, MC_allele) plus per-sample GT.
#   2. AnalyzeTRV: site BED (repeat span, REF repeat sequence/length/motif count, per-ALT length, size and motif-count
#      differences vs REF, AC/AF/AN, genic context coding > UTR > intronic > intergenic with gene names, and for
#      intergenic TRs the distance to and gene of the closest 5'UTR and 3'UTR); per-sample size-difference and
#      motif-count-difference matrices (site x sample, "d1,d2" in GT order); per-sample summary by genic context.
# Then the site BEDs are concatenated (bgzipped + tabix-indexed) and the per-sample summaries are summed.

workflow AnalyzeTRVariants {
    input {
        Array[File] vcfs
        Array[File] vcf_idxs
        Array[String] contigs
        String prefix

        Boolean all_filters = false

        File gtf

        File analysis_script

        String bcftools_docker
        String python_docker

        RuntimeAttr? runtime_attr_extract_trv
        RuntimeAttr? runtime_attr_analyze_trv
        RuntimeAttr? runtime_attr_concat_sites
        RuntimeAttr? runtime_attr_merge_summaries
    }

    scatter (i in range(length(contigs))) {
        call ExtractTRV {
            input:
                vcf = vcfs[i],
                vcf_idx = vcf_idxs[i],
                prefix = "~{prefix}.~{contigs[i]}",
                docker = bcftools_docker,
                runtime_attr_override = runtime_attr_extract_trv
        }

        call AnalyzeTRV {
            input:
                query_tsv = ExtractTRV.query_tsv,
                samples = ExtractTRV.samples,
                gtf = gtf,
                all_filters = all_filters,
                analysis_script = analysis_script,
                prefix = "~{prefix}.~{contigs[i]}",
                docker = python_docker,
                runtime_attr_override = runtime_attr_analyze_trv
        }
    }

    call ConcatSites {
        input:
            sites_beds = AnalyzeTRV.sites_bed,
            prefix = "~{prefix}.TRV.sites",
            docker = bcftools_docker,
            runtime_attr_override = runtime_attr_concat_sites
    }

    call MergeSummaries {
        input:
            summaries = AnalyzeTRV.per_sample_summary,
            analysis_script = analysis_script,
            prefix = "~{prefix}.TRV",
            docker = python_docker,
            runtime_attr_override = runtime_attr_merge_summaries
    }

    output {
        File sites_bed = ConcatSites.sites_bed
        File sites_bed_idx = ConcatSites.sites_bed_idx
        File per_sample_summary = MergeSummaries.per_sample_summary
        Array[File] trv_vcfs = ExtractTRV.trv_vcf
        Array[File] trv_vcf_idxs = ExtractTRV.trv_vcf_idx
        Array[File] sample_size_diff_matrices = AnalyzeTRV.sample_size_diff
        Array[File] sample_motif_count_diff_matrices = AnalyzeTRV.sample_motif_count_diff
        Array[File] per_contig_sample_summaries = AnalyzeTRV.per_sample_summary
    }
}

task ExtractTRV {
    input {
        File vcf
        File vcf_idx
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        bcftools view -i 'INFO/allele_type="trv"' -Oz -o ~{prefix}.TRV.vcf.gz ~{vcf}
        bcftools index -t ~{prefix}.TRV.vcf.gz
        bcftools query -l ~{prefix}.TRV.vcf.gz > ~{prefix}.samples.txt
        bcftools query \
            -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\t%FILTER\t%INFO/TRID\t%INFO/MOTIFS\t%INFO/AC\t%INFO/AF\t%INFO/AN\t%INFO/MC_allele[\t%GT]\n' \
            ~{prefix}.TRV.vcf.gz \
        | gzip > ~{prefix}.TRV.query.tsv.gz
    >>>

    output {
        File trv_vcf = "~{prefix}.TRV.vcf.gz"
        File trv_vcf_idx = "~{prefix}.TRV.vcf.gz.tbi"
        File samples = "~{prefix}.samples.txt"
        File query_tsv = "~{prefix}.TRV.query.tsv.gz"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 4,
        disk_gb: 2 * ceil(size(vcf, "GB")) + 10,
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

task AnalyzeTRV {
    input {
        File query_tsv
        File samples
        File gtf
        Boolean all_filters
        File analysis_script
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        python3 ~{analysis_script} analyze \
            --query-tsv ~{query_tsv} \
            --samples ~{samples} \
            --gtf ~{gtf} \
            --prefix ~{prefix} \
            ~{if all_filters then "--all-filters" else ""}
    >>>

    output {
        File sites_bed = "~{prefix}.TRV.sites.bed.gz"
        File sample_size_diff = "~{prefix}.TRV.sample_size_diff.tsv.gz"
        File sample_motif_count_diff = "~{prefix}.TRV.sample_motif_count_diff.tsv.gz"
        File per_sample_summary = "~{prefix}.TRV.per_sample_summary.tsv"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 8,
        disk_gb: 4 * ceil(size(query_tsv, "GB")) + ceil(size(gtf, "GB")) + 10,
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

task ConcatSites {
    input {
        Array[File] sites_beds
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        first=1
        for f in ~{sep=" " sites_beds}; do
            if [ "$first" -eq 1 ]; then
                gzip -dc "$f"
                first=0
            else
                gzip -dc "$f" | tail -n +2
            fi
        done | bgzip > ~{prefix}.bed.gz
        tabix -p bed ~{prefix}.bed.gz
    >>>

    output {
        File sites_bed = "~{prefix}.bed.gz"
        File sites_bed_idx = "~{prefix}.bed.gz.tbi"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 2,
        disk_gb: 3 * ceil(size(sites_beds, "GB")) + 10,
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

task MergeSummaries {
    input {
        Array[File] summaries
        File analysis_script
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        python3 ~{analysis_script} merge-summaries \
            --summaries ~{sep=" " summaries} \
            --out ~{prefix}.per_sample_summary.tsv
    >>>

    output {
        File per_sample_summary = "~{prefix}.per_sample_summary.tsv"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 2,
        disk_gb: 2 * ceil(size(summaries, "GB")) + 10,
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
