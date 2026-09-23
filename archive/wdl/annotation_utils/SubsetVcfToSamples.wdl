version 1.0

import "../utils/Structs.wdl"
import "../utils/Helpers.wdl"

workflow SubsetVcfToSamples {
    meta {
        description: [
            "This utility subsets a cohort VCF to a list of samples, one contig at a time, and concatenates the results."
        ]
    }

    parameter_meta {
        vcf: "Cohort VCF to subset."
        vcf_idx: "Index for `vcf`."
        samples: "Sample IDs to retain."
        contigs: "Contigs to process."
        subset_samples_vcf: "VCF containing only the requested samples."
        subset_samples_vcf_idx: "Index for `subset_samples_vcf`."
    }

    input {
        File vcf
        File vcf_idx
        Array[String] samples
        Array[String] contigs
        String prefix

        String utils_docker

        RuntimeAttr? runtime_attr_subset_vcf
        RuntimeAttr? runtime_attr_concat_vcfs
    }

    scatter (contig in contigs) {
        call Helpers.SubsetVcfToSampleList {
            input:
                vcf = vcf,
                vcf_idx = vcf_idx,
                samples = samples,
                extra_args = "--regions ~{contig}",
                prefix = "~{prefix}.~{contig}",
                docker = utils_docker,
                runtime_attr_override = runtime_attr_subset_vcf
        }
    }

    call Helpers.ConcatVcfs {
        input:
            vcfs = SubsetVcfToSampleList.subset_vcf,
            vcfs_idx = SubsetVcfToSampleList.subset_vcf_idx,
            merge_sort = true,
            prefix = "~{prefix}.subset_samples",
            docker = utils_docker,
            runtime_attr_override = runtime_attr_concat_vcfs
    }

    output {
        File subset_samples_vcf = ConcatVcfs.concat_vcf
        File subset_samples_vcf_idx = ConcatVcfs.concat_vcf_idx
    }
}
