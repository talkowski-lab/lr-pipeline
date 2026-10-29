version 1.0

import "../../../wdl/utils/Structs.wdl"
import "../../../wdl/utils/Helpers.wdl"

workflow SetMissingFiltersToPass {
    meta {
        description: [
            "This one-off utility sets the FILTER of every record without a FILTER value to `PASS`, for callsets built before `PreprocessVcfs` did so for the GLnexus records."
        ]
    }

    parameter_meta {
        vcf: "VCF to update."
        vcf_idx: "Index for `vcf`."
        pass_filled_vcf: "`vcf` with every missing FILTER value set to `PASS`."
        pass_filled_vcf_idx: "Index for `pass_filled_vcf`."
    }

    input {
        File vcf
        File vcf_idx
        String prefix

        String utils_docker

        RuntimeAttr? runtime_attr_set_missing_filters
    }

    call Helpers.SetMissingFiltersToPass as SetMissingFilters {
        input:
            vcf = vcf,
            vcf_idx = vcf_idx,
            prefix = "~{prefix}.annotated",
            docker = utils_docker,
            runtime_attr_override = runtime_attr_set_missing_filters
    }

    output {
        File pass_filled_vcf = SetMissingFilters.filtered_vcf
        File pass_filled_vcf_idx = SetMissingFilters.filtered_vcf_idx
    }
}
