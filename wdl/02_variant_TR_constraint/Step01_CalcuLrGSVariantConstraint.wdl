version 1.0

# Analyses_01 lrGS variant constraint. Per contig: a sites-only BCF of the genotype VCF, then
# (1) per-gene variant counts by type (SNV, indel 1-49 bp, SV), TR context (noTR / STR / VNTR), frequency (all / rare)
# and genic region / consequence (Analyses_01.lrGS_Variant_Constraint.count_gene_variants.py) and (2) counts of non-TR
# variants by type and size with the TR-enveloped share (Analyses_01.lrGS_Variant_Constraint.tr_enveloped_by_size.py).
# Then the per-contig counts are merged into gene x variant tables (merge_gene_tables.py) and the TR-enveloped tables
# summed (tr_enveloped_by_size.py --tables). Scripts are File inputs.

import "../utils/Structs.wdl"

workflow Step01_CalcuLrGSVariantConstraint {
    input {
        Array[File] vcfs
        Array[File] vcf_idxs
        Array[String] contigs
        String prefix

        Float rare_af = 0.01

        File gtf
        File count_script
        File tr_script
        File merge_script

        String docker

        RuntimeAttr? runtime_attr_count_contig
        RuntimeAttr? runtime_attr_merge_gene_tables
        RuntimeAttr? runtime_attr_sum_tr_enveloped
    }

    parameter_meta {
        vcfs: "Genotype (or sites) VCF per contig, with INFO vep, TRID, TR_ENVELOPED, allele_type, allele_length, SOURCE, AC, AF, PREDICTED_*"
        gtf: "GTF with gene records (gene_name) for gene rows and locations"
        count_script: "Analyses_01.lrGS_Variant_Constraint.count_gene_variants.py"
        tr_script: "Analyses_01.lrGS_Variant_Constraint.tr_enveloped_by_size.py"
        merge_script: "Analyses_01.lrGS_Variant_Constraint.merge_gene_tables.py"
    }

    scatter (i in range(length(contigs))) {
        call CountContig {
            input:
                vcf = vcfs[i],
                vcf_idx = vcf_idxs[i],
                prefix = prefix + "." + contigs[i],
                rare_af = rare_af,
                count_script = count_script,
                tr_script = tr_script,
                docker = docker,
                runtime_attr_override = runtime_attr_count_contig
        }
    }

    call MergeGeneTables {
        input:
            counts = CountContig.counts,
            prefix = prefix,
            gtf = gtf,
            merge_script = merge_script,
            docker = docker,
            runtime_attr_override = runtime_attr_merge_gene_tables
    }

    call SumTrEnveloped {
        input:
            tables = CountContig.tr_enveloped,
            prefix = prefix,
            tr_script = tr_script,
            docker = docker,
            runtime_attr_override = runtime_attr_sum_tr_enveloped
    }

    output {
        Array[File] gene_tables = MergeGeneTables.gene_tables
        File tr_enveloped_by_size = SumTrEnveloped.table
        Array[File] contig_counts = CountContig.counts
        Array[File] contig_tr_enveloped = CountContig.tr_enveloped
    }
}

task CountContig {
    input {
        File vcf
        File vcf_idx
        String prefix
        Float rare_af
        File count_script
        File tr_script
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        bcftools view -G -Ob -o ~{prefix}.sites.bcf ~{vcf}
        bcftools index ~{prefix}.sites.bcf
        python3 ~{count_script} --vcf ~{prefix}.sites.bcf --rare-af ~{rare_af} --out ~{prefix}.counts.tsv.gz
        python3 ~{tr_script} --vcf ~{prefix}.sites.bcf --out ~{prefix}.tr_enveloped_by_size.tsv
    >>>

    output {
        File counts = "~{prefix}.counts.tsv.gz"
        File tr_enveloped = "~{prefix}.tr_enveloped_by_size.tsv"
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

task MergeGeneTables {
    input {
        Array[File] counts
        String prefix
        File gtf
        File merge_script
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        python3 ~{merge_script} --inputs ~{sep=" " counts} --gtf ~{gtf} --out-prefix ~{prefix}
    >>>

    output {
        Array[File] gene_tables = glob("~{prefix}.*.gene_table.tsv")
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 8,
        disk_gb: ceil(size(counts, "GB")) + ceil(size(gtf, "GB")) + 10,
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

task SumTrEnveloped {
    input {
        Array[File] tables
        String prefix
        File tr_script
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        python3 ~{tr_script} --tables ~{sep=" " tables} --out ~{prefix}.tr_enveloped_by_size.tsv
    >>>

    output {
        File table = "~{prefix}.tr_enveloped_by_size.tsv"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 2,
        disk_gb: 10,
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
