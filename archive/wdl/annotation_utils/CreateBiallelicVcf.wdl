version 1.0

import "../utils/Helpers.wdl"
import "../utils/Structs.wdl"

workflow CreateBiallelicVcf {
    input {
        File vcf
        File vcf_idx
        String prefix

        File ref_fa
        File ref_fai

        String utils_docker

        RuntimeAttr? runtime_attr_normalize
        RuntimeAttr? runtime_attr_annotate_attributes
        RuntimeAttr? runtime_attr_rename
    }

    call Helpers.NormalizeVcf {
        input:
            vcf = vcf,
            vcf_idx = vcf_idx,
            ref_fa = ref_fa,
            ref_fai = ref_fai,
            prefix = "~{prefix}.normalized",
            docker = utils_docker,
            runtime_attr_override = runtime_attr_normalize
    }

    call Helpers.AnnotateVariantAttributes {
        input:
            vcf = NormalizeVcf.normalized_vcf,
            vcf_idx = NormalizeVcf.normalized_vcf_idx,
            prefix = "~{prefix}.annotated",
            docker = utils_docker,
            runtime_attr_override = runtime_attr_annotate_attributes
    }

    call RenameVariantIds {
        input:
            vcf = AnnotateVariantAttributes.annotated_vcf,
            vcf_idx = AnnotateVariantAttributes.annotated_vcf_idx,
            prefix = "~{prefix}.biallelic",
            docker = utils_docker,
            runtime_attr_override = runtime_attr_rename
    }

    output {
        File biallelic_vcf = RenameVariantIds.renamed_vcf
        File biallelic_vcf_idx = RenameVariantIds.renamed_vcf_idx
    }
}

task RenameVariantIds {
    input {
        File vcf
        File vcf_idx
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        python3 <<CODE
from pysam import VariantFile
from collections import defaultdict

def new_variant_id(record):
    a_type = record.info.get('allele_type').upper()
    a_len = abs(int(record.info.get('allele_length')))
    if a_type == "SNV":
        return f"{record.chrom}-{record.pos}-{record.ref}-{record.alts[0]}"
    return f"{record.chrom}-{record.pos}-{a_type}-{a_len}"

# First pass: count new IDs to detect collisions
vcf_in = VariantFile("~{vcf}")
id_counts = defaultdict(int)
for record in vcf_in:
    id_counts[new_variant_id(record)] += 1
vcf_in.close()

# Second pass: assign streamlined IDs, suffixing duplicates
vcf_in = VariantFile("~{vcf}")
vcf_out = VariantFile("~{prefix}.vcf.gz", "w", header=vcf_in.header)
id_seen = defaultdict(int)
for record in vcf_in:
    new_id = new_variant_id(record)
    if id_counts[new_id] > 1:
        id_seen[new_id] += 1
        record.id = f"{new_id}_{id_seen[new_id]}"
    else:
        record.id = new_id
    vcf_out.write(record)
vcf_in.close()
vcf_out.close()
CODE

        tabix -p vcf ~{prefix}.vcf.gz
    >>>

    output {
        File renamed_vcf = "~{prefix}.vcf.gz"
        File renamed_vcf_idx = "~{prefix}.vcf.gz.tbi"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 4,
        disk_gb: 2 * ceil(size(vcf, "GB")) + 5,
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
