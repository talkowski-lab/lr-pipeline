version 1.0

import "Structs.wdl"

# Kanpig writes a genotype-level FORMAT/FT field with integer values (0, 1, ...).
# Per the VCF spec, FT is a *reserved* per-sample genotype-filter key, so strict
# parsers (htsjdk / GATK SVAnnotate) interpret its value as a FILTER name and fail
# with e.g. "0 is an invalid filter name in vcf4". bcftools and pysam tolerate it,
# so the malformed field only surfaces once the VCF reaches GATK.
#
# This workflow rewrites the offending field. By default it *renames* FORMAT/FT to
# a non-reserved key (FTK) so Kanpig's per-genotype filter values are preserved but
# are no longer treated as genotype filters. Set drop_ft=true to strip it entirely.
workflow FixKanpigFT {
    input {
        File vcf
        File vcf_idx
        String prefix
        Boolean drop_ft = false
        String new_ft_name = "FTK"
        String docker
        RuntimeAttr? runtime_attr_override
    }

    call RenameOrDropFT {
        input:
            vcf = vcf,
            vcf_idx = vcf_idx,
            prefix = prefix,
            drop_ft = drop_ft,
            new_ft_name = new_ft_name,
            docker = docker,
            runtime_attr_override = runtime_attr_override
    }

    output {
        File fixed_vcf = RenameOrDropFT.fixed_vcf
        File fixed_vcf_idx = RenameOrDropFT.fixed_vcf_idx
    }
}

task RenameOrDropFT {
    input {
        File vcf
        File vcf_idx
        String prefix
        Boolean drop_ft
        String new_ft_name
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        # No-op safely if the VCF has no FORMAT/FT to begin with.
        if ! bcftools view -h ~{vcf} | grep -q '^##FORMAT=<ID=FT,'; then
            echo "No FORMAT/FT field present; copying input through unchanged." >&2
            cp ~{vcf} ~{prefix}.vcf.gz
        elif ~{if drop_ft then "true" else "false"}; then
            bcftools annotate -x FORMAT/FT -Oz -o ~{prefix}.vcf.gz ~{vcf}
        else
            printf 'FORMAT/FT\t~{new_ft_name}\n' > rename.txt
            bcftools annotate --rename-annots rename.txt -Oz -o ~{prefix}.vcf.gz ~{vcf}
        fi

        tabix -p vcf ~{prefix}.vcf.gz
    >>>

    output {
        File fixed_vcf = "~{prefix}.vcf.gz"
        File fixed_vcf_idx = "~{prefix}.vcf.gz.tbi"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 4,
        disk_gb: 2 * ceil(size(vcf, "GB")) + 10,
        boot_disk_gb: 10,
        preemptible_tries: 2,
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
