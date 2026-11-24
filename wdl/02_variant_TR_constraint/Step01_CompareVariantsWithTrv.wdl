version 1.0

# Compare, per tandem-repeat variant record (TRV, allele_type trv), the variants enveloped by it (INFO/TR_ENVELOPED,
# INFO/TRID = TRV record ID) with the variants obtained by splitting the TRV's own REF / ALT alleles with bcftools
# (INFO dropped, multi-allelics split into one record per ALT allele and left-aligned against the reference; no
# atomizing, so each TRV allele stays one variant). Per contig: a sites-only BCF, the TRV records split, the enveloped
# records taken as they are in the VCF, then
# Analyses_01.lrGS_Variant_Constraint.compare_variants_with_trv.py writes one row per TRV (repeat size, motifs, allele
# lengths, counts of SNV / Indel_del / Indel_ins / SV_del / SV_ins / other in each set and in both). The per-contig tables
# are concatenated. The script is a File input.

import "../utils/Structs.wdl"

workflow Step01_CompareVariantsWithTrv {
    input {
        Array[File] vcfs
        Array[File] vcf_idxs
        Array[String] contigs
        String prefix

        Boolean pass_only = false

        File ref_fa
        File ref_fai
        File compare_script

        String docker

        RuntimeAttr? runtime_attr_compare_contig
        RuntimeAttr? runtime_attr_concat_tables
    }

    parameter_meta {
        vcfs: "Callset VCF per contig with TRV records (allele_type trv, MOTIFS, TRID) and TR_ENVELOPED / TRID on other records"
        pass_only: "Use only FILTER PASS records (TRVs and enveloped variants)"
        ref_fa: "Reference FASTA used for bcftools norm (left-alignment, REF check)"
        compare_script: "Analyses_01.lrGS_Variant_Constraint.compare_variants_with_trv.py"
    }

    scatter (i in range(length(contigs))) {
        call CompareContig {
            input:
                vcf = vcfs[i],
                vcf_idx = vcf_idxs[i],
                prefix = prefix + "." + contigs[i],
                pass_only = pass_only,
                ref_fa = ref_fa,
                ref_fai = ref_fai,
                compare_script = compare_script,
                docker = docker,
                runtime_attr_override = runtime_attr_compare_contig
        }
    }

    call ConcatTables {
        input:
            tables = CompareContig.table,
            prefix = prefix,
            docker = docker,
            runtime_attr_override = runtime_attr_concat_tables
    }

    output {
        File trv_vs_enveloped = ConcatTables.table
        Array[File] contig_tables = CompareContig.table
        Array[File] contig_logs = CompareContig.log
    }
}

task CompareContig {
    input {
        File vcf
        File vcf_idx
        String prefix
        Boolean pass_only
        File ref_fa
        File ref_fai
        File compare_script
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        bcftools view -G -Ob -o sites.bcf ~{vcf}
        bcftools index sites.bcf
        bcftools view -i 'INFO/allele_type="trv"' -Ob -o trv.bcf sites.bcf
        bcftools annotate -x INFO -Ou trv.bcf \
            | bcftools norm -m-any -f ~{ref_fa} -Ob -o trv.norm.bcf 2> norm.log
        bcftools view -i 'INFO/TR_ENVELOPED=1 && INFO/allele_type!="trv"' -Ob -o env.bcf sites.bcf
        python3 ~{compare_script} --trv trv.bcf --trv-norm trv.norm.bcf --env env.bcf \
            --out ~{prefix}.trv_vs_enveloped.tsv.gz ~{if pass_only then "--pass-only" else ""} > compare.log
        cat norm.log compare.log > ~{prefix}.trv_vs_enveloped.log
    >>>

    output {
        File table = "~{prefix}.trv_vs_enveloped.tsv.gz"
        File log = "~{prefix}.trv_vs_enveloped.log"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 8,
        disk_gb: 2 * ceil(size(vcf, "GB")) + ceil(size(ref_fa, "GB")) + 10,
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

task ConcatTables {
    input {
        Array[File] tables
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        python3 - ~{prefix}.trv_vs_enveloped.tsv.gz ~{sep=" " tables} <<'EOF'
import gzip
import sys

out, tables = sys.argv[1], sys.argv[2:]
with gzip.open(out, "wt") as o:
    for i, t in enumerate(tables):
        with gzip.open(t, "rt") as f:
            head = f.readline()
            if i == 0:
                o.write(head)
            for line in f:
                o.write(line)
EOF
    >>>

    output {
        File table = "~{prefix}.trv_vs_enveloped.tsv.gz"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 2,
        disk_gb: 2 * ceil(size(tables, "GB")) + 10,
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
