version 1.0

import "../utils/Helpers.wdl"
import "../utils/Structs.wdl"

workflow ExtractSampleVcfs {
    meta {
        description: [
            "This utility extracts one sample's VCFs from a cohort VCF, splitting its variants into a SNV/indel VCF and an SV VCF based on a minimum SV length and the `INFO/SOURCE` tags set by PreprocessVcfs. It is run once per sample, and outputs the sample's SNV/indel and SV VCFs."
        ]
    }

    parameter_meta {
        cohort_vcf: "Cohort VCF to extract from."
        cohort_vcf_idx: "Index for the cohort VCF."
        contigs: "Contigs to process within the cohort VCF."
        sample_id: "Sample to extract."
        min_sv_length: "Minimum length at which a variant is routed to the SV VCF rather than the SNV/indel VCF."
        snv_indel_source_tag: "`INFO/SOURCE` of the calls routed to the SNV/indel VCF."
        sv_source_tag: "`INFO/SOURCE` of the SV calls routed to the SV VCF."
        long_sv_source_tag: "`INFO/SOURCE` of the long-SV calls routed to the SV VCF."
        snv_indel_vcf: "The sample's SNV/indel VCF."
        snv_indel_vcf_idx: "Index for the SNV/indel VCF."
        sv_vcf: "The sample's SV VCF."
        sv_vcf_idx: "Index for the SV VCF."
    }

    input {
        File cohort_vcf
        File cohort_vcf_idx
        Array[String] contigs
        String prefix

        String sample_id
        Int min_sv_length
        String snv_indel_source_tag = "DeepVariant"
        String sv_source_tag = "SV_Integration"
        String long_sv_source_tag = "LongSV_Integration"

        String utils_docker

        RuntimeAttr? runtime_attr_subset_contig
        RuntimeAttr? runtime_attr_extract_sample
        RuntimeAttr? runtime_attr_subset_snv
        RuntimeAttr? runtime_attr_subset_sv
        RuntimeAttr? runtime_attr_concat_snv
        RuntimeAttr? runtime_attr_concat_sv
    }

    # Cromwell rejects an undefined optional entering a nested scatter, so these cross as one-element arrays

    Boolean single_contig = length(contigs) == 1

    scatter (contig in contigs) {
        if (!single_contig) {
            call Helpers.SubsetVcfToContig as SubsetContig {
                input:
                    vcf = cohort_vcf,
                    vcf_idx = cohort_vcf_idx,
                    contig = contig,
                    prefix = "~{prefix}.~{contig}",
                    docker = utils_docker,
                    runtime_attr_override = runtime_attr_subset_contig
            }
        }

        File contig_vcf = select_first([SubsetContig.subset_vcf, cohort_vcf])
        File contig_vcf_idx = select_first([SubsetContig.subset_vcf_idx, cohort_vcf_idx])

        call Helpers.ExtractSample {
            input:
                vcf = contig_vcf,
                vcf_idx = contig_vcf_idx,
                sample = sample_id,
                prefix = "~{prefix}.~{contig}",
                docker = utils_docker,
                runtime_attr_override = runtime_attr_extract_sample
        }

        call Helpers.SubsetVcfByArgs as SubsetSnvIndel {
            input:
                vcf = ExtractSample.subset_vcf,
                vcf_idx = ExtractSample.subset_vcf_idx,
                include_args = 'abs(INFO/allele_length) < ~{min_sv_length} && INFO/SOURCE = \"~{snv_indel_source_tag}\"',
                prefix = "~{prefix}.~{contig}.snv",
                docker = utils_docker,
                runtime_attr_override = runtime_attr_subset_snv
        }

        call Helpers.SubsetVcfByArgs as SubsetSv {
            input:
                vcf = ExtractSample.subset_vcf,
                vcf_idx = ExtractSample.subset_vcf_idx,
                include_args = 'abs(INFO/allele_length) >= ~{min_sv_length} && (INFO/SOURCE = \"~{sv_source_tag}\" || INFO/SOURCE = \"~{long_sv_source_tag}\")',
                prefix = "~{prefix}.~{contig}.sv",
                docker = utils_docker,
                runtime_attr_override = runtime_attr_subset_sv
        }
    }

    call Helpers.ConcatVcfs as ConcatSnv {
        input:
            vcfs = SubsetSnvIndel.subset_vcf,
            vcf_idxs = SubsetSnvIndel.subset_vcf_idx,
            allow_overlaps = false,
            naive = true,
            prefix = "~{prefix}.snv",
            docker = utils_docker,
            runtime_attr_override = runtime_attr_concat_snv
    }

    call Helpers.ConcatVcfs as ConcatSv {
        input:
            vcfs = SubsetSv.subset_vcf,
            vcf_idxs = SubsetSv.subset_vcf_idx,
            allow_overlaps = false,
            naive = true,
            prefix = "~{prefix}.sv",
            docker = utils_docker,
            runtime_attr_override = runtime_attr_concat_sv
    }

    output {
        File snv_indel_vcf = ConcatSnv.concat_vcf
        File snv_indel_vcf_idx = ConcatSnv.concat_vcf_idx
        File sv_vcf = ConcatSv.concat_vcf
        File sv_vcf_idx = ConcatSv.concat_vcf_idx
    }
}
