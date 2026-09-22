version 1.0

import "../utils/Helpers.wdl"
import "../utils/Structs.wdl"

workflow ExtractDisparateTRLoci {
    meta {
        description: [
            "This utility subsets two VCFs to tandem-repeat variants (`INFO/allele_type=trv`) on one contig, then compares their loci. It produces one TSV for identities present in only one VCF, where identity is `CHROM`, `POS` and `len(REF)`, and another TSV for positive-base overlaps with distinct identities. Overlaps are identified with `bedtools intersect`; the overlapping TSV includes the `INFO/TRID` value from both VCFs."
        ]
    }

    parameter_meta {
        vcf_a: "First VCF to compare."
        vcf_a_idx: "Index for `vcf_a`."
        vcf_b: "Second VCF to compare."
        vcf_b_idx: "Index for `vcf_b`."
        contig: "Contig to compare within both VCFs."
        missing_variants_tsv: "Locus identities present in one VCF but missing from the other."
        overlapping_variants_tsv: "Overlapping locus pairs with distinct identities and their `TRID` values."
    }

    input {
        File vcf_a
        File vcf_a_idx
        File vcf_b
        File vcf_b_idx
        String contig
        String prefix

        String utils_docker

        RuntimeAttr? runtime_attr_subset_vcf
        RuntimeAttr? runtime_attr_extract_disparate_loci
    }

    call Helpers.SubsetVcfToContig as SubsetVcfA {
        input:
            vcf = vcf_a,
            vcf_idx = vcf_a_idx,
            contig = contig,
            extra_args = "-i 'INFO/allele_type=\"trv\"'",
            prefix = "~{prefix}.~{contig}.a",
            docker = utils_docker,
            runtime_attr_override = runtime_attr_subset_vcf
    }

    call Helpers.SubsetVcfToContig as SubsetVcfB {
        input:
            vcf = vcf_b,
            vcf_idx = vcf_b_idx,
            contig = contig,
            extra_args = "-i 'INFO/allele_type=\"trv\"'",
            prefix = "~{prefix}.~{contig}.b",
            docker = utils_docker,
            runtime_attr_override = runtime_attr_subset_vcf
    }

    call ExtractDisparateTRLociForContig {
        input:
            vcf_a = SubsetVcfA.subset_vcf,
            vcf_b = SubsetVcfB.subset_vcf,
            prefix = "~{prefix}.~{contig}",
            docker = utils_docker,
            runtime_attr_override = runtime_attr_extract_disparate_loci
    }

    output {
        File missing_variants_tsv = ExtractDisparateTRLociForContig.missing_variants_tsv
        File overlapping_variants_tsv = ExtractDisparateTRLociForContig.overlapping_variants_tsv
    }
}

task ExtractDisparateTRLociForContig {
    input {
        File vcf_a
        File vcf_b
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        bcftools query \
            -f '%CHROM\t%POS0\t%POS\t%REF\t%INFO/TRID\n' \
            ~{vcf_a} \
            | awk 'BEGIN {OFS="\t"} {print $1, $2, $2 + length($4), $3, length($4), $5}' \
            > a.bed

        bcftools query \
            -f '%CHROM\t%POS0\t%POS\t%REF\t%INFO/TRID\n' \
            ~{vcf_b} \
            | awk 'BEGIN {OFS="\t"} {print $1, $2, $2 + length($4), $3, length($4), $5}' \
            > b.bed

        {
            bedtools intersect \
                -a a.bed \
                -b b.bed \
                -f 1.0 \
                -r \
                -v \
                | awk 'BEGIN {OFS="\t"} {print "vcf_a", "vcf_b", $1, $4, $5}'
            bedtools intersect \
                -a b.bed \
                -b a.bed \
                -f 1.0 \
                -r \
                -v \
                | awk 'BEGIN {OFS="\t"} {print "vcf_b", "vcf_a", $1, $4, $5}'
        } | sort -u > missing_variants.rows.tsv

        {
            printf 'present_in\tmissing_from\tchrom\tpos\tref_length\n'
            cat missing_variants.rows.tsv
        } > ~{prefix}.missing_variants.tsv

        bedtools intersect \
            -a a.bed \
            -b b.bed \
            -wa \
            -wb \
            | awk 'BEGIN {OFS="\t"} $1 != $7 || $4 != $10 || $5 != $11 {
                print $1, $4, $5, $6, $7, $10, $11, $12
            }' \
            | sort -u > overlapping_variants.rows.tsv

        {
            printf 'chrom_a\tpos_a\tref_length_a\ttrid_a\tchrom_b\tpos_b\tref_length_b\ttrid_b\n'
            cat overlapping_variants.rows.tsv
        } > ~{prefix}.overlapping_variants.tsv
    >>>

    output {
        File missing_variants_tsv = "~{prefix}.missing_variants.tsv"
        File overlapping_variants_tsv = "~{prefix}.overlapping_variants.tsv"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 4,
        disk_gb: 4 * ceil(size([vcf_a, vcf_b], "GB")) + 10,
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
