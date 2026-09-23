version 1.0

import "../utils/Structs.wdl"
import "../utils/Helpers.wdl" as Helpers

workflow MergeVEPAF {
    meta {
        description: [
            "This utility combines a VEP-annotated VCF with an allele-frequency-annotated VCF, transferring the VEP consequence field onto the allele-frequency callset one contig at a time and merging the results."
        ]
    }

    parameter_meta {
        af_annotation_vcf: "VCF carrying the allele-frequency annotations."
        af_annotation_vcf_idx: "Index for `af_annotation_vcf`."
        vep_annotation_vcf: "VCF carrying the VEP annotations."
        vep_annotation_vcf_idx: "Index for `vep_annotation_vcf`."
        contigs: "Contigs to process."
        vep_info_field_name: "INFO field holding the VEP consequence string."
        merged_vcf: "VCF carrying both annotation sets."
        merged_vcf_idx: "Index for `merged_vcf`."
    }

    input {
        File af_annotation_vcf
        File af_annotation_vcf_idx
        File vep_annotation_vcf
        File vep_annotation_vcf_idx
        Array[String] contigs
        String prefix
        
        String vep_info_field_name
        
        String utils_docker
        
        RuntimeAttr? runtime_attr_subset_af
        RuntimeAttr? runtime_attr_subset_vep
        RuntimeAttr? runtime_attr_merge_annotations
        RuntimeAttr? runtime_attr_merge_vcfs
    }

    scatter (contig in contigs) {
        call Helpers.SubsetVcfToContig as SubsetAF {
            input:
                vcf = af_annotation_vcf,
                vcf_idx = af_annotation_vcf_idx,
                contig = contig,
                prefix = "~{prefix}.~{contig}.af",
                docker = utils_docker,
                runtime_attr_override = runtime_attr_subset_af
        }

        call Helpers.SubsetVcfToContig as SubsetVEP {
            input:
                vcf = vep_annotation_vcf,
                vcf_idx = vep_annotation_vcf_idx,
                contig = contig,
                prefix = "~{prefix}.~{contig}.vep",
                docker = utils_docker,
                runtime_attr_override = runtime_attr_subset_vep
        }

        call MergeVEPAnnotation {
            input:
                af_vcf = SubsetAF.subset_vcf,
                af_vcf_idx = SubsetAF.subset_vcf_idx,
                vep_vcf = SubsetVEP.subset_vcf,
                vep_vcf_idx = SubsetVEP.subset_vcf_idx,
                vep_info_field_name = vep_info_field_name,
                prefix = "~{prefix}.~{contig}",
                docker = utils_docker,
                runtime_attr_override = runtime_attr_merge_annotations
        }
    }

    call Helpers.ConcatVcfs as MergeVcfs {
        input:
            vcfs = MergeVEPAnnotation.merged_vcf,
            vcfs_idx = MergeVEPAnnotation.merged_vcf_idx,
            prefix = prefix,
            docker = utils_docker,
            runtime_attr_override = runtime_attr_merge_vcfs
    }

    output {
        File merged_vcf = MergeVcfs.concat_vcf
        File merged_vcf_idx = MergeVcfs.concat_vcf_idx
    }
}

task MergeVEPAnnotation {
    input {
        File af_vcf
        File af_vcf_idx
        File vep_vcf
        File vep_vcf_idx
        String vep_info_field_name
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        bcftools view -h ~{vep_vcf} | grep "^##INFO=<ID=~{vep_info_field_name}," > vep_header.txt

        bcftools query -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\t%INFO/~{vep_info_field_name}\n' ~{vep_vcf} | bgzip -c > vep_annots.tab.gz
        tabix -s 1 -b 2 -e 2 vep_annots.tab.gz

        bcftools annotate \
            -a vep_annots.tab.gz \
            -h vep_header.txt \
            -c CHROM,POS,ID,REF,ALT,~{vep_info_field_name} \
            -Oz -o ~{prefix}.merged.vcf.gz \
            ~{af_vcf}
        
        tabix -p vcf -f ~{prefix}.merged.vcf.gz
    >>>

    output {
        File merged_vcf = "~{prefix}.merged.vcf.gz"
        File merged_vcf_idx = "~{prefix}.merged.vcf.gz.tbi"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 2,
        mem_gb: 4,
        disk_gb: 2 * ceil(size(af_vcf, "GB") + size(vep_vcf, "GB")) + 10,
        boot_disk_gb: 10,
        preemptible_tries: 1,
        max_retries: 0
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
        memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
        disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
        bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
        docker: docker
        preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
    }
}
