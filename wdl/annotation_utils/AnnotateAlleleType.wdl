version 1.0

import "../utils/Helpers.wdl"
import "../utils/Structs.wdl"

workflow AnnotateAlleleType {
    input {
        File vcf
        File vcf_idx
        File med_tsv
        File mei_tsv
        File dup_tsv
        Array[String] contigs
        String prefix

        String? med_prefix
        String? med_suffix
        Boolean? med_lowercase

        String? mei_prefix
        String? mei_suffix
        Boolean? mei_lowercase

        String? dup_prefix
        String? dup_suffix
        Boolean? dup_lowercase

        String utils_docker

        RuntimeAttr? runtime_attr_subset_vcf
        RuntimeAttr? runtime_attr_subset_tsv
        RuntimeAttr? runtime_attr_annotate_med
        RuntimeAttr? runtime_attr_annotate_mei
        RuntimeAttr? runtime_attr_annotate_dup
        RuntimeAttr? runtime_attr_concat
    }

    Boolean single_contig = length(contigs) == 1

    scatter (contig in contigs) {
        if (!single_contig) {
            call Helpers.SubsetVcfToContig {
                input:
                    vcf = vcf,
                    vcf_idx = vcf_idx,
                    contig = contig,
                    prefix = "~{prefix}.~{contig}",
                    docker = utils_docker,
                    runtime_attr_override = runtime_attr_subset_vcf
            }

            call Helpers.SubsetTsvToContig as SubsetMedTsv {
                input:
                    tsv = med_tsv,
                    contig = contig,
                    prefix = "~{prefix}.~{contig}.med",
                    docker = utils_docker,
                    runtime_attr_override = runtime_attr_subset_tsv
            }

            call Helpers.SubsetTsvToContig as SubsetMeiTsv {
                input:
                    tsv = mei_tsv,
                    contig = contig,
                    prefix = "~{prefix}.~{contig}.mei",
                    docker = utils_docker,
                    runtime_attr_override = runtime_attr_subset_tsv
            }

            call Helpers.SubsetTsvToContig as SubsetDupTsv {
                input:
                    tsv = dup_tsv,
                    contig = contig,
                    prefix = "~{prefix}.~{contig}.dup",
                    docker = utils_docker,
                    runtime_attr_override = runtime_attr_subset_tsv
            }
        }

        File contig_vcf = select_first([SubsetVcfToContig.subset_vcf, vcf])
        File contig_vcf_idx = select_first([SubsetVcfToContig.subset_vcf_idx, vcf_idx])
        File contig_med_tsv = select_first([SubsetMedTsv.subset_tsv, med_tsv])
        File contig_mei_tsv = select_first([SubsetMeiTsv.subset_tsv, mei_tsv])
        File contig_dup_tsv = select_first([SubsetDupTsv.subset_tsv, dup_tsv])

        call AddHeaders {
            input:
                vcf = contig_vcf,
                vcf_idx = contig_vcf_idx,
                prefix = "~{prefix}.~{contig}.headers",
                docker = utils_docker,
                runtime_attr_override = runtime_attr_annotate_med
        }

        call AnnotateDup {
            input:
                vcf = AddHeaders.vcf_with_headers,
                vcf_idx = AddHeaders.vcf_with_headers_idx,
                dup_tsv = contig_dup_tsv,
                dup_prefix = dup_prefix,
                dup_suffix = dup_suffix,
                dup_lowercase = dup_lowercase,
                prefix = "~{prefix}.~{contig}.dup_annotated",
                docker = utils_docker,
                runtime_attr_override = runtime_attr_annotate_dup
        }

        call AnnotateMei {
            input:
                vcf = AnnotateDup.annotated_vcf,
                vcf_idx = AnnotateDup.annotated_vcf_idx,
                mei_tsv = contig_mei_tsv,
                mei_prefix = mei_prefix,
                mei_suffix = mei_suffix,
                mei_lowercase = mei_lowercase,
                prefix = "~{prefix}.~{contig}.mei_annotated",
                docker = utils_docker,
                runtime_attr_override = runtime_attr_annotate_mei
        }

        call AnnotateMed {
            input:
                vcf = AnnotateMei.annotated_vcf,
                vcf_idx = AnnotateMei.annotated_vcf_idx,
                med_tsv = contig_med_tsv,
                med_prefix = med_prefix,
                med_suffix = med_suffix,
                med_lowercase = med_lowercase,
                prefix = "~{prefix}.~{contig}.med_annotated",
                docker = utils_docker,
                runtime_attr_override = runtime_attr_annotate_med
        }
    }

    if (!single_contig) {
        call Helpers.ConcatVcfs {
            input:
                vcfs = AnnotateMed.annotated_vcf,
                vcf_idxs = AnnotateMed.annotated_vcf_idx,
                allow_overlaps = false,
                naive = true,
                prefix = "~{prefix}.allele_type_annotated",
                docker = utils_docker,
                runtime_attr_override = runtime_attr_concat
        }
    }

    output {
        File allele_type_annotated_vcf = select_first([ConcatVcfs.concat_vcf, AnnotateMed.annotated_vcf[0]])
        File allele_type_annotated_vcf_idx = select_first([ConcatVcfs.concat_vcf_idx, AnnotateMed.annotated_vcf_idx[0]])
    }
}

task AddHeaders {
    input {
        File vcf
        File vcf_idx
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        bcftools view -h "~{vcf}" | grep "^##" > header.txt

        echo '##INFO=<ID=SUB_FAMILY,Number=1,Type=String,Description="Sub-family of mobile element.">' >> header.txt
        echo '##INFO=<ID=ORIGIN,Number=.,Type=String,Description="Genomic coordinates for the source of the inserted sequence.">' >> header.txt

        bcftools view -h "~{vcf}" | grep "^#CHROM" >> header.txt

        bcftools reheader \
            -h header.txt \
            "~{vcf}" \
            > ~{prefix}.vcf.gz

        tabix -p vcf ~{prefix}.vcf.gz
    >>>

    output {
        File vcf_with_headers = "~{prefix}.vcf.gz"
        File vcf_with_headers_idx = "~{prefix}.vcf.gz.tbi"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 2,
        disk_gb: 2 * ceil(size(vcf, "GB")) + 10,
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

task AnnotateDup {
    input {
        File vcf
        File vcf_idx
        File dup_tsv
        String? dup_prefix
        String? dup_suffix
        Boolean? dup_lowercase
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        dup_prefix=~{select_first([dup_prefix, ''])}
        dup_suffix=~{select_first([dup_suffix, ''])}
        dup_lowercase=~{select_first([dup_lowercase, 'false'])}

        awk -v pre="$dup_prefix" -v suf="$dup_suffix" -v lower="$dup_lowercase" '
            BEGIN {
                FS = OFS = "\t"
            }
            {
                gsub(/[\n\r]/, "", $6)
                gsub(/[\n\r]/, "", $7)
                gsub(/[\n\r]/, "", $8)

                $6 = pre $6 suf
                if (lower == "true") {
                    $6 = tolower($6)
                }

                if ($7 != "" && $7 != "." && length($7) > 0) {
                    origin = $7
                } else {
                    origin = $8
                }

                print $1, $2, $3, $4, $5, $6, origin
            }
        ' "~{dup_tsv}" > "dup_modified.tsv"

        bgzip -c "dup_modified.tsv" > "dup_annotations.tsv.gz"
        tabix -s1 -b2 -e2 "dup_annotations.tsv.gz"

        bcftools annotate \
            -a "dup_annotations.tsv.gz" \
            -c CHROM,POS,REF,ALT,~ID,INFO/allele_type,INFO/ORIGIN \
            -Oz -o ~{prefix}.vcf.gz \
            "~{vcf}"

        tabix -p vcf ~{prefix}.vcf.gz
    >>>

    output {
        File annotated_vcf = "~{prefix}.vcf.gz"
        File annotated_vcf_idx = "~{prefix}.vcf.gz.tbi"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 4,
        disk_gb: 2 * ceil(size(vcf, "GB") + size(dup_tsv, "GB")) + 10,
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

task AnnotateMei {
    input {
        File vcf
        File vcf_idx
        File mei_tsv
        String? mei_prefix
        String? mei_suffix
        Boolean? mei_lowercase
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        mei_prefix=~{select_first([mei_prefix, ''])}
        mei_suffix=~{select_first([mei_suffix, ''])}
        mei_lowercase=~{select_first([mei_lowercase, 'false'])}

        awk -v pre="$mei_prefix" -v suf="$mei_suffix" -v lower="$mei_lowercase" '
            BEGIN {
                FS = OFS = "\t"
            }
            {
                gsub(/[\n\r]/, "", $6)
                gsub(/[\n\r]/, "", $7)
                $6 = pre $6 suf
                if (lower == "true") {
                    $6 = tolower($6)
                }
                print
            }
        ' "~{mei_tsv}" > "mei_modified.tsv"

        bgzip -c "mei_modified.tsv" > "mei_annotations.tsv.gz"
        tabix -s1 -b2 -e2 "mei_annotations.tsv.gz"

        bcftools annotate \
            -a "mei_annotations.tsv.gz" \
            -c CHROM,POS,REF,ALT,~ID,INFO/allele_type,INFO/SUB_FAMILY \
            -Oz -o ~{prefix}.vcf.gz \
            "~{vcf}"

        tabix -p vcf ~{prefix}.vcf.gz
    >>>

    output {
        File annotated_vcf = "~{prefix}.vcf.gz"
        File annotated_vcf_idx = "~{prefix}.vcf.gz.tbi"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 4,
        disk_gb: 2 * ceil(size(vcf, "GB") + size(mei_tsv, "GB")) + 10,
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

task AnnotateMed {
    input {
        File vcf
        File vcf_idx
        File med_tsv
        String? med_prefix
        String? med_suffix
        Boolean? med_lowercase
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        med_prefix=~{select_first([med_prefix, ''])}
        med_suffix=~{select_first([med_suffix, ''])}
        med_lowercase=~{select_first([med_lowercase, 'false'])}

        awk -v pre="$med_prefix" -v suf="$med_suffix" -v lower="$med_lowercase" '
            BEGIN {
                FS = OFS = "\t"
            }
            {
                gsub(/[\n\r]/, "", $6)
                gsub(/[\n\r]/, "", $7)
                $6 = pre $6 suf
                if (lower == "true") {
                    $6 = tolower($6)
                }
                print
            }
        ' "~{med_tsv}" > "med_modified.tsv"

        bgzip -c "med_modified.tsv" > "med_annotations.tsv.gz"
        tabix -s1 -b2 -e2 "med_annotations.tsv.gz"

        bcftools annotate \
            -a "med_annotations.tsv.gz" \
            -c CHROM,POS,REF,ALT,~ID,INFO/allele_type,INFO/SUB_FAMILY \
            -Oz -o ~{prefix}.vcf.gz \
            "~{vcf}"

        tabix -p vcf ~{prefix}.vcf.gz
    >>>

    output {
        File annotated_vcf = "~{prefix}.vcf.gz"
        File annotated_vcf_idx = "~{prefix}.vcf.gz.tbi"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 4,
        disk_gb: 2 * ceil(size(vcf, "GB") + size(med_tsv, "GB")) + 10,
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
