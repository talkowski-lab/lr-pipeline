version 1.0

import "../utils/Helpers.wdl"
import "../utils/Structs.wdl"

workflow SubsetVcfToContigs {
    meta {
        description: [
            "This utility subsets a VCF to a chosen set of contigs and concatenates the result. It outputs the subset VCF."
        ]
    }

    parameter_meta {
        vcf: "VCF to subset."
        vcf_idx: "Index for VCF."
        contigs: "Contigs to retain."
        subset_contigs_vcf: "Contig-subset VCF."
        subset_contigs_vcf_idx: "Index for the subset VCF."
    }

    input {
        File vcf
        File vcf_idx
        Array[String] contigs
        String prefix

        String utils_docker

        RuntimeAttr? runtime_attr_subset_vcf
        RuntimeAttr? runtime_attr_concat_vcf
    }

    scatter (contig in contigs) {
        call Helpers.SubsetVcfToContig {
            input:
                vcf = vcf,
                vcf_idx = vcf_idx,
                contig = contig,
                prefix = "~{prefix}.~{contig}",
                docker = utils_docker,
                runtime_attr_override = runtime_attr_subset_vcf
        }
    }

    call Helpers.ConcatVcfs {
        input:
            vcfs = SubsetVcfToContig.subset_vcf,
            vcf_idxs = SubsetVcfToContig.subset_vcf_idx,
            allow_overlaps = false,
            naive = true,
            prefix = "~{prefix}.subset_contigs",
            docker = utils_docker,
            runtime_attr_override = runtime_attr_concat_vcf
    }

    output {
        File subset_contigs_vcf = ConcatVcfs.concat_vcf
        File subset_contigs_vcf_idx = ConcatVcfs.concat_vcf_idx
    }
}
