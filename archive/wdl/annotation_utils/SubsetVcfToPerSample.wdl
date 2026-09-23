version 1.0

import "../utils/Helpers.wdl"
import "../utils/Structs.wdl"

workflow SubsetVcfToPerSample {
    meta {
        description: [
            "This utility extracts a separate single-sample VCF for each requested sample from a set of cohort VCFs, optionally dropping specified fields first. It outputs the per-sample VCFs."
        ]
    }

    parameter_meta {
        cohort_vcfs: "Cohort VCFs to extract from."
        cohort_vcf_idxs: "Indexes for `cohort_vcfs`."
        contigs: "Contigs to process."
        sample_ids: "Samples to extract."
        drop_fields: "Fields to drop from each VCF before extraction."
        subset_vcfs: "Per-sample VCFs."
        subset_vcf_idxs: "Indexes for the per-sample VCFs."
    }

    input {
        Array[File] cohort_vcfs
        Array[File] cohort_vcf_idxs
        Array[String] contigs
        String prefix

        Array[String] sample_ids
        String? drop_fields

        String utils_docker

        RuntimeAttr? runtime_attr_drop_fields
        RuntimeAttr? runtime_attr_extract_sample
        RuntimeAttr? runtime_attr_concat_vcfs
    }

    scatter (i in range(length(contigs))) {
        String contig = contigs[i]

        if (defined(drop_fields)) {
            call Helpers.DropVcfFields {
                input:
                    vcf = cohort_vcfs[i],
                    vcf_idx = cohort_vcf_idxs[i],
                    drop_fields = select_first([drop_fields]),
                    prefix = "~{prefix}.~{contig}.dropped",
                    docker = utils_docker,
                    runtime_attr_override = runtime_attr_drop_fields
            }
        }

        File base_vcf = select_first([DropVcfFields.dropped_vcf, cohort_vcfs[i]])
        File base_vcf_idx = select_first([DropVcfFields.dropped_vcf_idx, cohort_vcf_idxs[i]])

        scatter (sample_id in sample_ids) {
            call Helpers.ExtractSample {
                input:
                    vcf = base_vcf,
                    vcf_idx = base_vcf_idx,
                    sample = sample_id,
                    normalize_output = true,
                    prefix = "~{prefix}.~{sample_id}.~{contig}",
                    docker = utils_docker,
                    runtime_attr_override = runtime_attr_extract_sample
            }
        }
    }

    # Transpose from [contig][sample] to [sample][contig] for per-sample concat
    Array[Array[File]] vcfs_by_sample = transpose(ExtractSample.subset_vcf)
    Array[Array[File]] vcf_idxs_by_sample = transpose(ExtractSample.subset_vcf_idx)

    scatter (i in range(length(sample_ids))) {
        call Helpers.ConcatVcfs {
            input:
                vcfs = vcfs_by_sample[i],
                vcf_idxs = vcf_idxs_by_sample[i],
                allow_overlaps = true,
                naive = false,
                prefix = "~{prefix}.~{sample_ids[i]}",
                docker = utils_docker,
                runtime_attr_override = runtime_attr_concat_vcfs
        }
    }

    output {
        Array[File] subset_vcfs = ConcatVcfs.concat_vcf
        Array[File] subset_vcf_idxs = ConcatVcfs.concat_vcf_idx
    }
}
