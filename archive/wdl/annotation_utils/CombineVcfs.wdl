version 1.0

import "../utils/Helpers.wdl"
import "../utils/Structs.wdl"

workflow CombineVcfs {
    meta {
        description: [
            "This utility concatenates two VCFs holding different variants for the same samples. The sample lists are checked for a match, each contig is subset from both inputs, and the per-contig results are concatenated into one VCF."
        ]
    }

    parameter_meta {
        a_vcf: "First VCF to combine."
        a_vcf_idx: "Index for `a_vcf`."
        b_vcf: "Second VCF to combine."
        b_vcf_idx: "Index for `b_vcf`."
        contigs: "Contigs to process."
        concat_vcf: "Combined VCF."
        concat_vcf_idx: "Index for the combined VCF."
    }

    input {
        File a_vcf
        File a_vcf_idx
        File b_vcf
        File b_vcf_idx
        Array[String] contigs
        String prefix

        String utils_docker

        RuntimeAttr? runtime_attr_check_samples
        RuntimeAttr? runtime_attr_subset_a
        RuntimeAttr? runtime_attr_subset_b
        RuntimeAttr? runtime_attr_merge
        RuntimeAttr? runtime_attr_concat
        RuntimeAttr? runtime_attr_concat_merged
        RuntimeAttr? runtime_attr_concat_concatenated
    }

    call Helpers.CheckSampleMatch {
        input:
            vcf_a = a_vcf,
            vcf_a_idx = a_vcf_idx,
            vcf_b = b_vcf,
            vcf_b_idx = b_vcf_idx,
            docker = utils_docker,
            runtime_attr_override = runtime_attr_check_samples
    }

    scatter (contig in contigs) {
        call Helpers.SubsetVcfToContig as SubsetA {
            input:
                vcf = a_vcf,
                vcf_idx = a_vcf_idx,
                contig = contig,
                prefix = "~{prefix}.~{contig}.a",
                docker = utils_docker,
                runtime_attr_override = runtime_attr_subset_a
        }

        call Helpers.SubsetVcfToContig as SubsetB {
            input:
                vcf = b_vcf,
                vcf_idx = b_vcf_idx,
                contig = contig,
                prefix = "~{prefix}.~{contig}.b",
                docker = utils_docker,
                runtime_attr_override = runtime_attr_subset_b
        }

        call Helpers.ConcatVcfs {
            input:
                vcfs = [SubsetA.subset_vcf, SubsetB.subset_vcf],
                vcf_idxs = [SubsetA.subset_vcf_idx, SubsetB.subset_vcf_idx],
                prefix = "~{prefix}.~{contig}.concatenated",
                docker = utils_docker,
                runtime_attr_override = runtime_attr_concat
        }
    }

    call Helpers.ConcatVcfs as ConcatContigs {
        input:
            vcfs = ConcatVcfs.concat_vcf,
            vcf_idxs = ConcatVcfs.concat_vcf_idx,
            prefix = "~{prefix}.concatenated",
            docker = utils_docker,
            runtime_attr_override = runtime_attr_concat_concatenated
    }

    output {
        File concat_vcf = ConcatContigs.concat_vcf
        File concat_vcf_idx = ConcatContigs.concat_vcf_idx
    }
}
