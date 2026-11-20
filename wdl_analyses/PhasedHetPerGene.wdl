version 1.0

# Phased heterozygous coding (LoF, missense, synonymous) and UTR (5', 3') variants per gene and per sample. Per contig:
# count_phased_het_per_gene.py on the genotype VCF (with the per-variant table); then the contigs are merged
# (merge_phased_het_counts.py) and summarized per gene: samples with >= 1 / 3 / 5 phased hets of each category, with
# gene locations from a GTF (summarize_phased_het_per_gene.py). Scripts are File inputs.

import "Structs.wdl"

workflow PhasedHetPerGene {
    input {
        Array[File] vcfs
        Array[File] vcf_idxs
        Array[String] contigs
        String prefix

        String transcripts = "canonical"

        File gtf
        File count_script
        File merge_script
        File summarize_script

        String docker

        RuntimeAttr? runtime_attr_count_phased_het
        RuntimeAttr? runtime_attr_merge_phased_het
    }

    parameter_meta {
        vcfs: "Genotype VCF per contig (phased GT, PS), with INFO vep, PREDICTED_LOF, allele_type, SOURCE"
        transcripts: "VEP transcripts used: canonical, mane or all"
        gtf: "GTF with gene records (gene_name) for gene locations"
        count_script: "count_phased_het_per_gene.py"
        merge_script: "merge_phased_het_counts.py"
        summarize_script: "summarize_phased_het_per_gene.py"
    }

    scatter (i in range(length(contigs))) {
        call CountPhasedHet {
            input:
                vcf = vcfs[i],
                vcf_idx = vcf_idxs[i],
                prefix = prefix + "." + contigs[i],
                transcripts = transcripts,
                count_script = count_script,
                docker = docker,
                runtime_attr_override = runtime_attr_count_phased_het
        }
    }

    call MergePhasedHet {
        input:
            per_gene_sample = CountPhasedHet.per_gene_sample,
            per_sample = CountPhasedHet.per_sample,
            variants = CountPhasedHet.variants,
            prefix = prefix,
            gtf = gtf,
            merge_script = merge_script,
            summarize_script = summarize_script,
            docker = docker,
            runtime_attr_override = runtime_attr_merge_phased_het
    }

    output {
        File per_gene_sample = MergePhasedHet.merged_per_gene_sample
        File per_sample = MergePhasedHet.merged_per_sample
        File variants = MergePhasedHet.merged_variants
        File per_gene_summary = MergePhasedHet.per_gene_summary
    }
}

task CountPhasedHet {
    input {
        File vcf
        File vcf_idx
        String prefix
        String transcripts
        File count_script
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        python3 ~{count_script} --vcf ~{vcf} --transcripts ~{transcripts} --variants --out-prefix ~{prefix}
    >>>

    output {
        File per_gene_sample = "~{prefix}.per_gene_sample.tsv.gz"
        File per_sample = "~{prefix}.per_sample.tsv"
        File variants = "~{prefix}.variants.tsv.gz"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 8,
        disk_gb: ceil(size(vcf, "GB")) + 10,
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

task MergePhasedHet {
    input {
        Array[File] per_gene_sample
        Array[File] per_sample
        Array[File] variants
        String prefix
        File gtf
        File merge_script
        File summarize_script
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        mkdir -p contigs
        for f in ~{sep=" " per_gene_sample} ~{sep=" " per_sample} ~{sep=" " variants}; do
            ln -s "$f" contigs/
        done
        inputs=$(for f in ~{sep=" " per_sample}; do b=$(basename "$f"); echo "contigs/${b%.per_sample.tsv}"; done)
        python3 ~{merge_script} --inputs ${inputs} --out-prefix ~{prefix}
        python3 ~{summarize_script} --in ~{prefix}.per_gene_sample.tsv.gz --gtf ~{gtf} \
            --out ~{prefix}.per_gene_summary.tsv
    >>>

    output {
        File merged_per_gene_sample = "~{prefix}.per_gene_sample.tsv.gz"
        File merged_per_sample = "~{prefix}.per_sample.tsv"
        File merged_variants = "~{prefix}.variants.tsv.gz"
        File per_gene_summary = "~{prefix}.per_gene_summary.tsv"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 16,
        disk_gb: 3 * ceil(size(variants, "GB") + size(per_gene_sample, "GB")) + 10,
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
