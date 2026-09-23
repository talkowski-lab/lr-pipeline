version 1.0

import "../utils/Helpers.wdl"
import "../utils/Structs.wdl"

workflow FilterTRGTCalls {
    meta {
        description: [
            "This utility filters a TRGT tandem-repeat VCF, optionally dropping calls below a minimum repeat-unit length or length difference, or above a maximum catalog length. It outputs the filtered VCF."
        ]
    }

    parameter_meta {
        trgt_vcf: "TRGT VCF to filter."
        trgt_vcf_idx: "Index for the TRGT VCF."
        min_repeat_unit: "Minimum repeat-unit length to retain a call."
        min_length_diff: "Minimum length difference from the reference to retain a call."
        max_catalog_length: "Maximum catalog locus length to retain a call."
        trgt_filtered_vcf: "Filtered TRGT VCF."
        trgt_filtered_vcf_idx: "Index for the filtered VCF."
    }

    input {
        File trgt_vcf
        File trgt_vcf_idx
        String prefix

        Int? min_repeat_unit
        Int? min_length_diff
        Int? max_catalog_length

        String utils_docker

        RuntimeAttr? runtime_attr_filter
    }

    call Helpers.FilterTRGTVcf {
        input:
            vcf = trgt_vcf,
            vcf_idx = trgt_vcf_idx,
            min_repeat_unit = min_repeat_unit,
            min_length_diff = min_length_diff,
            max_catalog_length = max_catalog_length,
            prefix = prefix,
            docker = utils_docker,
            runtime_attr_override = runtime_attr_filter
    }

    output {
        File trgt_filtered_vcf = FilterTRGTVcf.processed_vcf
        File trgt_filtered_vcf_idx = FilterTRGTVcf.processed_vcf_idx
    }
}
