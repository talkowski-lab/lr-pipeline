version 1.0

import "../utils/Helpers.wdl"
import "../utils/Structs.wdl"

workflow AnnotateTREndTags {
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

        call AddEndTagContig {
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

task AddEndTagContig {
    input {
        File vcf
        File vcf_idx
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        bcftools view ~{vcf} \
            | awk 'BEGIN{FS=OFS="\t"; has_end_header=0}
                /^##INFO=<ID=END,/ {
                    has_end_header=1
                    print
                    next
                }
                /^#CHROM/ {
                    if (!has_end_header) {
                        print "##INFO=<ID=END,Number=1,Type=Integer,Description=\"End position of the variant described in this record\">"
                    }
                    print
                    next
                }
                /^#/ {
                    print
                    next
                }
                {
                    end_val=$2 + length($4) - 1

                    if ($8=="." || $8=="") {
                        $8="END=" end_val
                    } else if ($8 ~ /(^|;)END=/) {
                        n=split($8, info_parts, ";")
                        for (i=1; i<=n; i++) {
                            if (info_parts[i] ~ /^END=/) {
                                info_parts[i]="END=" end_val
                            }
                        }
                        $8=info_parts[1]
                        for (i=2; i<=n; i++) {
                            $8=$8 ";" info_parts[i]
                        }
                    } else {
                        $8=$8 ";END=" end_val
                    }

                    print
                }' \
            | bgzip -c > ~{prefix}.vcf.gz

        tabix -f -p vcf ~{prefix}.vcf.gz
    >>>

    output {
        File vcf_with_end = "~{prefix}.vcf.gz"
        File vcf_with_end_idx = "~{prefix}.vcf.gz.tbi"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 4,
        disk_gb: 3 * ceil(size([vcf, vcf_idx], "GB")) + 5,
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
