version 1.0

import "../utils/Helpers.wdl"
import "../utils/Structs.wdl"

workflow TRGTLPS {
    input {
        File vcf
        File vcf_idx
        Array[String] contigs
        String prefix

        Boolean normalize_chry_haploid_genotypes

        String trgt_lps_docker
        String utils_docker
        String stranalysis_docker

        RuntimeAttr? runtime_attr_subset_vcf
        RuntimeAttr? runtime_attr_normalize_chry_haploid_genotypes
        RuntimeAttr? runtime_attr_add_end
        RuntimeAttr? runtime_attr_trgt_lps
        RuntimeAttr? runtime_attr_extract_trid_metadata
        RuntimeAttr? runtime_attr_concat
    }

    scatter (contig in contigs) {
        Boolean normalize_chry_haploid_genotypes_for_contig = normalize_chry_haploid_genotypes && contig == "chrY"

        call Helpers.SubsetVcfToContig {
            input:
                vcf = vcf,
                vcf_idx = vcf_idx,
                contig = contig,
                prefix = "~{prefix}.~{contig}",
                docker = utils_docker,
                runtime_attr_override = runtime_attr_subset_vcf
        }

        if (normalize_chry_haploid_genotypes_for_contig) {
            call Helpers.NormalizeTRGTHaploidGenotypes {
                input:
                    vcf = SubsetVcfToContig.subset_vcf,
                    vcf_idx = SubsetVcfToContig.subset_vcf_idx,
                    prefix = "~{prefix}.~{contig}.haploid_genotypes",
                    docker = utils_docker,
                    runtime_attr_override = runtime_attr_normalize_chry_haploid_genotypes
            }
        }

        File trgt_lps_vcf = select_first([NormalizeTRGTHaploidGenotypes.normalized_vcf, SubsetVcfToContig.subset_vcf])
        File trgt_lps_vcf_idx = select_first([NormalizeTRGTHaploidGenotypes.normalized_vcf_idx, SubsetVcfToContig.subset_vcf_idx])

        call RunTRGTLPS {
            input:
                vcf = trgt_lps_vcf,
                vcf_idx = trgt_lps_vcf_idx,
                prefix = "~{prefix}.~{contig}.lps",
                docker = trgt_lps_docker,
                runtime_attr_override = runtime_attr_trgt_lps
        }

        call Helpers.AddTREndTag {
            input:
                vcf = SubsetVcfToContig.subset_vcf,
                vcf_idx = SubsetVcfToContig.subset_vcf_idx,
                prefix = "~{prefix}.~{contig}.with_end",
                docker = utils_docker,
                runtime_attr_override = runtime_attr_add_end
        }

        call ExtractTridMetadata {
            input:
                vcf = AddTREndTag.vcf_with_end,
                vcf_idx = AddTREndTag.vcf_with_end_idx,
                contig = contig,
                prefix = "~{prefix}.~{contig}.trid_metadata",
                docker = stranalysis_docker,
                runtime_attr_override = runtime_attr_extract_trid_metadata
        }
    }

    call Helpers.ConcatTsvs {
        input:
            tsvs = RunTRGTLPS.lps_tsv,
            sort_output = false,
            preserve_header = true,
            prefix = "~{prefix}.concat",
            docker = utils_docker,
            runtime_attr_override = runtime_attr_concat
    }

    output {
        File trgt_lps_tsv = ConcatTsvs.concatenated_tsv
        Array[File] vcf_trid_metadata_tsvs = ExtractTridMetadata.trid_metadata_tsv
    }
}

task RunTRGTLPS {
    input {
        File vcf
        File vcf_idx
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        trgt-lps \
            --vcf ~{vcf} \
            --threads ~{select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])} \
            > ~{prefix}.tsv
    >>>

    output {
        File lps_tsv = "~{prefix}.tsv"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 8,
        mem_gb: 4,
        disk_gb: 2 * ceil(size(vcf, "GB")) + 20,
        boot_disk_gb: 10,
        preemptible_tries: 1,
        max_retries: 0
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
        memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GB"
        disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " SSD"
        bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
        docker: docker
        preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
    }
}

task ExtractTridMetadata {
    input {
        File vcf
        File vcf_idx
        String contig
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        ln -s ~{vcf} input.vcf.gz
        ln -s ~{vcf_idx} input.vcf.gz.tbi

        python3 -m str_analysis.extract_trid_metadata_from_TRGT_vcf \
            --input-vcf input.vcf.gz \
            --output-tsv ~{prefix}.tsv.gz \
            --contig ~{contig}
    >>>

    output {
        File trid_metadata_tsv = "~{prefix}.tsv.gz"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 2,
        mem_gb: 8,
        disk_gb: ceil(size(vcf, "GB")) + 20,
        boot_disk_gb: 10,
        preemptible_tries: 1,
        max_retries: 0
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
        memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GB"
        disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
        bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
        docker: docker
        preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
    }
}
