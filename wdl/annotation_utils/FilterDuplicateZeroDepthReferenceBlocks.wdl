version 1.0

import "../utils/Structs.wdl"
import "../utils/Helpers.wdl"

workflow FilterDuplicateZeroDepthReferenceBlocks {
    meta {
        description: [
            "This utility cleans a single-sample gVCF by removing exact duplicate zero-depth, non-alt records, except that it retains one representative when removal would leave its start uncovered. It preserves gVCF coverage: a retained duplicate block is shortened by updating its `END` to one base before the next non-duplicate record when that record begins inside the block. This prevents cleanup from overlapping a distinct record or creating a coverage gap that GLNexus would genotype as `./.`. Singleton records, distinct records at the same coordinate, alternate genotypes, and records with non-zero or missing `MIN_DP` are retained unchanged."
        ]
    }

    parameter_meta {
        gvcf: "Single-sample gVCF to clean."
        gvcf_idx: "Index for `gvcf`."
        cleaned_vcf: "Cleaned gVCF."
        cleaned_vcf_idx: "Index for `cleaned_vcf`."
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
