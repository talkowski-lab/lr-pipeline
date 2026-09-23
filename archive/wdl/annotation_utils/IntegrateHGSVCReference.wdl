version 1.0

import "../utils/Helpers.wdl"
import "../utils/Structs.wdl"

workflow IntegrateHGSVCReference {
    meta {
        description: [
            "This utility merges the HGSVC SNV, indel and SV reference callsets into one VCF. Each input is optionally sample-swapped, checked for a consistent sample list, subset per contig and tagged with its own source label before the three are merged."
        ]
    }

    parameter_meta {
        snv_vcf: "HGSVC SNV callset."
        snv_vcf_idx: "Index for `snv_vcf`."
        indel_vcf: "HGSVC indel callset."
        indel_vcf_idx: "Index for `indel_vcf`."
        sv_vcf: "HGSVC SV callset."
        sv_vcf_idx: "Index for `sv_vcf`."
        contigs: "Contigs to process."
        sample_ids: "Sample IDs expected in every input callset."
        sample_swap_list: "Two-column file mapping original sample IDs to their replacements."
        snv_source_tag: "Source label applied to variants from `snv_vcf`."
        snv_source_tag_description: "Header description for `snv_source_tag`."
        indel_source_tag: "Source label applied to variants from `indel_vcf`."
        indel_source_tag_description: "Header description for `indel_source_tag`."
        sv_source_tag: "Source label applied to variants from `sv_vcf`."
        sv_source_tag_description: "Header description for `sv_source_tag`."
        integrated_reference_vcf: "Merged reference callset."
        integrated_reference_vcf_idx: "Index for `integrated_reference_vcf`."
    }

    input {
        File snv_vcf
        File snv_vcf_idx
        File indel_vcf
        File indel_vcf_idx
        File sv_vcf
        File sv_vcf_idx
        Array[String] contigs
        String prefix

        Array[String] sample_ids
        File? sample_swap_list
        String snv_source_tag
        String snv_source_tag_description
        String indel_source_tag
        String indel_source_tag_description
        String sv_source_tag
        String sv_source_tag_description
        
        String utils_docker

        RuntimeAttr? runtime_attr_swap
        RuntimeAttr? runtime_attr_check_samples
        RuntimeAttr? runtime_attr_subset
        RuntimeAttr? runtime_attr_annotate
        RuntimeAttr? runtime_attr_add_info
        RuntimeAttr? runtime_attr_merge
        RuntimeAttr? runtime_attr_concat
    }

    if (defined(sample_swap_list)) {
        call Helpers.SwapSampleIds as SwapSnv {
            input:
                vcf = snv_vcf,
                vcf_idx = snv_vcf_idx,
                sample_swap_list = select_first([sample_swap_list]),
                prefix = "~{prefix}.snv.swapped",
                docker = utils_docker,
                runtime_attr_override = runtime_attr_swap
        }

        call Helpers.SwapSampleIds as SwapIndel {
            input:
                vcf = indel_vcf,
                vcf_idx = indel_vcf_idx,
                sample_swap_list = select_first([sample_swap_list]),
                prefix = "~{prefix}.indel.swapped",
                docker = utils_docker,
                runtime_attr_override = runtime_attr_swap
        }
        
        call Helpers.SwapSampleIds as SwapSv {
            input:
                vcf = sv_vcf,
                vcf_idx = sv_vcf_idx,
                sample_swap_list = select_first([sample_swap_list]),
                prefix = "~{prefix}.sv.swapped",
                docker = utils_docker,
                runtime_attr_override = runtime_attr_swap
        }
    }

    File snv_vcf_final = select_first([SwapSnv.swapped_vcf, snv_vcf])
    File snv_vcf_idx_final = select_first([SwapSnv.swapped_vcf_idx, snv_vcf_idx])
    File indel_vcf_final = select_first([SwapIndel.swapped_vcf, indel_vcf])
    File indel_vcf_idx_final = select_first([SwapIndel.swapped_vcf_idx, indel_vcf_idx])
    File sv_vcf_final = select_first([SwapSv.swapped_vcf, sv_vcf])
    File sv_vcf_idx_final = select_first([SwapSv.swapped_vcf_idx, sv_vcf_idx])

    call Helpers.CheckSampleConsistency {
        input:
            vcfs = [snv_vcf_final, indel_vcf_final, sv_vcf_final],
            vcfs_idx = [snv_vcf_idx_final, indel_vcf_idx_final, sv_vcf_idx_final],
            sample_ids = sample_ids,
            docker = utils_docker,
            runtime_attr_override = runtime_attr_check_samples
    }

    scatter (contig in contigs) {
        # SNV VCF Processing
        call Helpers.SubsetVcfToSampleList as SubsetSnv {
            input:
                vcf = snv_vcf_final,
                vcf_idx = snv_vcf_idx_final,
                samples = sample_ids,
                contig = contig,
                prefix = "~{prefix}.~{contig}.snv.subset",
                docker = utils_docker,
                runtime_attr_override = runtime_attr_subset
        }

        call Helpers.AnnotateVariantAttributes as AnnotateSnv {
            input:
                vcf = SubsetSnv.subset_vcf,
                vcf_idx = SubsetSnv.subset_vcf_idx,
                prefix = "~{prefix}.~{contig}.snv.annotated",
                docker = utils_docker,
                runtime_attr_override = runtime_attr_annotate
        }

        call Helpers.AddInfo as AddInfoSnv {
            input:
                vcf = AnnotateSnv.annotated_vcf,
                vcf_idx = AnnotateSnv.annotated_vcf_idx,
                tag_id = "SOURCE",
                tag_value = snv_source_tag,
                tag_description = snv_source_tag_description,
                prefix = "~{prefix}.~{contig}.snv.source",
                docker = utils_docker,
                runtime_attr_override = runtime_attr_add_info
        }

        # Indel VCF Processing
        call Helpers.SubsetVcfToSampleList as SubsetIndel {
            input:
                vcf = indel_vcf_final,
                vcf_idx = indel_vcf_idx_final,
                samples = sample_ids,
                contig = contig,
                prefix = "~{prefix}.~{contig}.indel.subset",
                docker = utils_docker,
                runtime_attr_override = runtime_attr_subset
        }

        call Helpers.AnnotateVariantAttributes as AnnotateIndel {
            input:
                vcf = SubsetIndel.subset_vcf,
                vcf_idx = SubsetIndel.subset_vcf_idx,
                prefix = "~{prefix}.~{contig}.indel.annotated",
                docker = utils_docker,
                runtime_attr_override = runtime_attr_annotate
        }

        call Helpers.AddInfo as AddInfoIndel {
            input:
                vcf = AnnotateIndel.annotated_vcf,
                vcf_idx = AnnotateIndel.annotated_vcf_idx,
                tag_id = "SOURCE",
                tag_value = indel_source_tag,
                tag_description = indel_source_tag_description,
                prefix = "~{prefix}.~{contig}.indel.source",
                docker = utils_docker,
                runtime_attr_override = runtime_attr_add_info
        }

        # SV VCF Processing
        call Helpers.SubsetVcfToSampleList as SubsetSv {
            input:
                vcf = sv_vcf_final,
                vcf_idx = sv_vcf_idx_final,
                samples = sample_ids,
                contig = contig,
                prefix = "~{prefix}.~{contig}.sv.subset",
                docker = utils_docker,
                runtime_attr_override = runtime_attr_subset
        }

        call Helpers.AnnotateVariantAttributes as AnnotateSv {
            input:
                vcf = SubsetSv.subset_vcf,
                vcf_idx = SubsetSv.subset_vcf_idx,
                prefix = "~{prefix}.~{contig}.sv.annotated",
                docker = utils_docker,
                runtime_attr_override = runtime_attr_annotate
        }

        call Helpers.AddInfo as AddInfoSv {
            input:
                vcf = AnnotateSv.annotated_vcf,
                vcf_idx = AnnotateSv.annotated_vcf_idx,
                tag_id = "SOURCE",
                tag_value = sv_source_tag,
                tag_description = sv_source_tag_description,
                prefix = "~{prefix}.~{contig}.sv.source",
                docker = utils_docker,
                runtime_attr_override = runtime_attr_add_info
        }

        # Merging
        call Helpers.ConcatVcfs as MergeContigVcfs {
            input:
                vcfs = [AddInfoSnv.annotated_vcf, AddInfoIndel.annotated_vcf, AddInfoSv.annotated_vcf],
                vcfs_idx = [AddInfoSnv.annotated_vcf_idx, AddInfoIndel.annotated_vcf_idx, AddInfoSv.annotated_vcf_idx],
                prefix = "~{prefix}.~{contig}.integrated",
                docker = utils_docker,
                runtime_attr_override = runtime_attr_merge
        }
    }
    
    call Helpers.ConcatVcfs {
        input:
            vcfs = MergeContigVcfs.concat_vcf,
            vcfs_idx = MergeContigVcfs.concat_vcf_idx,
            prefix = prefix,
            docker = utils_docker,
            runtime_attr_override = runtime_attr_concat
    }

    output {
        File integrated_reference_vcf = ConcatVcfs.concat_vcf
        File integrated_reference_vcf_idx = ConcatVcfs.concat_vcf_idx
    }
}