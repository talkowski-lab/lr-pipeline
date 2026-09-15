version 1.0

import "../utils/Structs.wdl"
import "../utils/Helpers.wdl"

workflow FilterDuplicateZeroDepthReferenceBlocks {
    meta {
        description: "Collapse duplicate zero-depth non-alt reference blocks from a single gVCF while preserving coverage."
    }

    parameter_meta {
        gvcf: "Input gVCF to clean."
        gvcf_idx: "Index corresponding to input gVCF."
        prefix: "Prefix for cleaned gVCF and index."
        utils_docker: "Docker image containing bcftools, bgzip, and tabix."
        runtime_attr_filter: "Override runtime attributes for duplicate reference block filtering."
    }

    input {
        File gvcf
        File gvcf_idx
        String prefix
        String utils_docker
        RuntimeAttr? runtime_attr_filter
    }

    call Helpers.FilterDuplicateZeroDepthReferenceBlocks as FilterDuplicateZeroDepthReferenceBlocksTask {
        input:
            gvcf = gvcf,
            gvcf_idx = gvcf_idx,
            prefix = prefix,
            create_indexes = true,
            docker = utils_docker,
            runtime_attr_override = runtime_attr_filter
    }

    output {
        File cleaned_vcf = FilterDuplicateZeroDepthReferenceBlocksTask.cleaned_gvcfs[0]
        File cleaned_vcf_idx = FilterDuplicateZeroDepthReferenceBlocksTask.cleaned_gvcf_idxs[0]
    }
}
