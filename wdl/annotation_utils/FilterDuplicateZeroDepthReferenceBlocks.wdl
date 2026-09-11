version 1.0

import "../utils/Structs.wdl"

workflow FilterDuplicateZeroDepthReferenceBlocks {
    meta {
        description: "Remove duplicate zero-depth non-alt reference blocks from a single gVCF."
    }

    parameter_meta {
        gvcf: "Input gVCF to clean."
        gvcf_idx: "Index corresponding to the input gVCF."
        prefix: "Prefix for the cleaned gVCF and its index."
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

    call FilterDuplicateZeroDepthReferenceBlocksTask {
        input:
            gvcf = gvcf,
            gvcf_idx = gvcf_idx,
            prefix = prefix,
            docker = utils_docker,
            runtime_attr_override = runtime_attr_filter
    }

    output {
        File cleaned_vcf = FilterDuplicateZeroDepthReferenceBlocksTask.cleaned_vcf
        File cleaned_vcf_idx = FilterDuplicateZeroDepthReferenceBlocksTask.cleaned_vcf_idx
    }
}

task FilterDuplicateZeroDepthReferenceBlocksTask {
    input {
        File gvcf
        File gvcf_idx
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = 1 + 2 * ceil(size(gvcf, "GB"))

    command <<<
        set -euo pipefail

        test -s ~{gvcf_idx}

        bcftools view ~{gvcf} | awk -F$'\t' '
            BEGIN { removed_count = 0 }
            function clear_group(    i) {
                for (i = 1; i <= row_count; i++) delete rows[i]
                for (i in duplicate_count) delete duplicate_count[i]
                for (i in removable) delete removable[i]
                row_count = 0
            }
            function flush_group(    i, line) {
                for (i = 1; i <= row_count; i++) {
                    line = rows[i]
                    if (removable[line] && duplicate_count[line] > 1) {
                        removed_count++
                    } else {
                        print line
                    }
                }
            }
            function is_removable_record(    format_fields, sample_fields, field_count, sample_count, i, gt_index, min_dp_index, gt) {
                field_count = split($9, format_fields, ":")
                gt_index = 0
                min_dp_index = 0
                for (i = 1; i <= field_count; i++) {
                    if (format_fields[i] == "GT") gt_index = i
                    if (format_fields[i] == "MIN_DP") min_dp_index = i
                }
                if (gt_index == 0 || min_dp_index == 0 || NF < 10) return 0
                sample_count = split($10, sample_fields, ":")
                if (sample_count < gt_index || sample_count < min_dp_index) return 0
                gt = sample_fields[gt_index]
                return sample_fields[min_dp_index] == "0" && gt ~ /^(0|\.)([\/|](0|\.))*$/
            }
            /^#/ { print; next }
            {
                coordinate = $1 SUBSEP $2
                if (row_count > 0 && coordinate != current_coordinate) {
                    flush_group()
                    clear_group()
                }
                current_coordinate = coordinate
                rows[++row_count] = $0
                duplicate_count[$0]++
                if (is_removable_record()) removable[$0] = 1
            }
            END {
                if (row_count > 0) flush_group()
                print "Removed " removed_count " duplicate zero-depth non-alt gVCF records" > "/dev/stderr"
            }
        ' | bgzip > ~{prefix}.cleaned.g.vcf.gz
        tabix -p vcf ~{prefix}.cleaned.g.vcf.gz
    >>>

    output {
        File cleaned_vcf = "~{prefix}.cleaned.g.vcf.gz"
        File cleaned_vcf_idx = "~{prefix}.cleaned.g.vcf.gz.tbi"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 1,
        disk_gb: disk_size,
        boot_disk_gb: 25,
        preemptible_tries: 1,
        max_retries: 0
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
        memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
        disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " SSD"
        bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
        preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
        docker: docker
    }
}
