version 1.0

import "../utils/Helpers.wdl"
import "../utils/Structs.wdl"

workflow AnnotateTREndTags {
    meta {
        description: [
            "This utility adds an `END` INFO tag to the tandem-repeat records of a VCF, computed per contig, so that downstream tools correctly interpret the span of each TR call. It outputs the updated VCF."
        ]
    }

    parameter_meta {
        vcf: "VCF to update."
        vcf_idx: "Index for VCF to update."
        contigs: "Contigs to process within the input VCF."
        vcf_with_end: "VCF with `END` tags added to tandem-repeat records."
        vcf_with_end_idx: "Index for the updated VCF."
    }

    input {
        File vcf
        File vcf_idx
        Array[String] contigs
        String prefix

        String utils_docker

        RuntimeAttr? runtime_attr_subset_contig
        RuntimeAttr? runtime_attr_add_end
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
                runtime_attr_override = runtime_attr_subset_contig
        }

        call Helpers.AddTREndTag as AddEndTagContig {
            input:
                vcf = SubsetVcfToContig.subset_vcf,
                vcf_idx = SubsetVcfToContig.subset_vcf_idx,
                prefix = "~{prefix}.~{contig}.with_end",
                docker = utils_docker,
                runtime_attr_override = runtime_attr_add_end
        }
    }

    call Helpers.ConcatVcfs {
        input:
            vcfs = AddEndTagContig.vcf_with_end,
            vcf_idxs = AddEndTagContig.vcf_with_end_idx,
            allow_overlaps = false,
            naive = true,
            prefix = "~{prefix}.with_end",
            docker = utils_docker,
            runtime_attr_override = runtime_attr_concat_vcf
    }

    output {
        File vcf_with_end = ConcatVcfs.concat_vcf
        File vcf_with_end_idx = ConcatVcfs.concat_vcf_idx
    }
}
