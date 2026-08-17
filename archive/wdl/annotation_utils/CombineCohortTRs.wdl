version 1.0

import "../utils/Helpers.wdl"
import "../utils/Structs.wdl"

workflow CombineCohortTRs {
    input {
        Array[File] tr_vcfs
        Array[File] tr_vcf_idxs
        Array[String] tr_callers
        Array[String] contigs
        String prefix

        Array[String] sample_ids

        String utils_docker

        RuntimeAttr? runtime_attr_check_sample_consistency
        RuntimeAttr? runtime_attr_subset_contig
        RuntimeAttr? runtime_attr_set_missing_filters
        RuntimeAttr? runtime_attr_tag_tr_vcf
        RuntimeAttr? runtime_attr_set_tr_ids
        RuntimeAttr? runtime_attr_deduplicate_trs
        RuntimeAttr? runtime_attr_priority_merge
        RuntimeAttr? runtime_attr_concat_vcf
    }

    call Helpers.CheckSampleConsistency {
        input:
            vcfs = tr_vcfs,
            vcf_idxs = tr_vcf_idxs,
            sample_ids = sample_ids,
            docker = utils_docker,
            runtime_attr_override = runtime_attr_check_sample_consistency
    }

    scatter (contig in contigs) {
        scatter (i in range(length(tr_vcfs))) {
            call Helpers.SubsetVcfToContig as SubsetContigTr {
                input:
                    vcf = tr_vcfs[i],
                    vcf_idx = tr_vcf_idxs[i],
                    contig = contig,
                    prefix = "~{prefix}.~{contig}.tr~{i}",
                    docker = utils_docker,
                    runtime_attr_override = runtime_attr_subset_contig
            }

            call Helpers.SetMissingFiltersToPass {
                input:
                    vcf = SubsetContigTr.subset_vcf,
                    vcf_idx = SubsetContigTr.subset_vcf_idx,
                    prefix = "~{prefix}.~{contig}.tr~{i}.refiltered",
                    docker = utils_docker,
                    runtime_attr_override = runtime_attr_set_missing_filters
            }

            call TagTRVcf {
                input:
                    vcf = SetMissingFiltersToPass.filtered_vcf,
                    vcf_idx = SetMissingFiltersToPass.filtered_vcf_idx,
                    tr_caller = tr_callers[i],
                    prefix = "~{prefix}.~{contig}.tr~{i}.tagged",
                    docker = utils_docker,
                    runtime_attr_override = runtime_attr_tag_tr_vcf
            }

            call SetTRVariantIds {
                input:
                    vcf = TagTRVcf.tagged_vcf,
                    vcf_idx = TagTRVcf.tagged_vcf_idx,
                    prefix = "~{prefix}.~{contig}.tr~{i}.ids",
                    docker = utils_docker,
                    runtime_attr_override = runtime_attr_set_tr_ids
            }

            call DeduplicateOverlappingVariants {
                input:
                    vcf = SetTRVariantIds.renamed_vcf,
                    vcf_idx = SetTRVariantIds.renamed_vcf_idx,
                    prefix = "~{prefix}.~{contig}.tr~{i}.dedup",
                    docker = utils_docker,
                    runtime_attr_override = runtime_attr_deduplicate_trs
            }
        }

        call PriorityMergeTRVcfs {
            input:
                vcfs = DeduplicateOverlappingVariants.dedup_vcf,
                vcf_idxs = DeduplicateOverlappingVariants.dedup_vcf_idx,
                prefix = "~{prefix}.~{contig}.tr_merged",
                docker = utils_docker,
                runtime_attr_override = runtime_attr_priority_merge
        }
    }

    call Helpers.ConcatVcfs {
        input:
            vcfs = PriorityMergeTRVcfs.merged_vcf,
            vcf_idxs = PriorityMergeTRVcfs.merged_vcf_idx,
            allow_overlaps = false,
            naive = true,
            prefix = "~{prefix}.combined_trs",
            docker = utils_docker,
            runtime_attr_override = runtime_attr_concat_vcf
    }

    output {
        File trgt_merged_vcf = ConcatVcfs.concat_vcf
        File trgt_merged_vcf_idx = ConcatVcfs.concat_vcf_idx
    }
}

task SetTRVariantIds {
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

vcf_in = VariantFile("~{vcf}")
id_counts = defaultdict(int)
for record in vcf_in:
    new_id = f"{record.chrom}-{record.pos}-TRV-{len(record.ref)}"
    id_counts[new_id] += 1
vcf_in.close()

vcf_in = VariantFile("~{vcf}")
vcf_out = VariantFile("~{prefix}.vcf.gz", "w", header=vcf_in.header)
id_seen = defaultdict(int)
for record in vcf_in:
    new_id = f"{record.chrom}-{record.pos}-TRV-{len(record.ref)}"
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

task TagTRVcf {
    input {
        File vcf
        File vcf_idx
        String tr_caller
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        touch new_headers.txt
        if ! bcftools view -h ~{vcf} | grep -q '##INFO=<ID=allele_type'; then
            echo '##INFO=<ID=allele_type,Number=1,Type=String,Description="Allele type">' >> new_headers.txt
        fi
        if ! bcftools view -h ~{vcf} | grep -q '##INFO=<ID=SOURCE'; then
            echo '##INFO=<ID=SOURCE,Number=1,Type=String,Description="Source of variant call">' >> new_headers.txt
        fi

        bcftools annotate \
            -h new_headers.txt \
            -Oz -o tr_header_added.vcf.gz \
            ~{vcf}

        rm -f new_headers.txt

        bcftools view tr_header_added.vcf.gz \
            | awk -v source="~{tr_caller}" 'BEGIN{OFS="\t"} /^#/ {print; next} {
                $8 = $8 ";allele_type=trv;SOURCE=" source
                print
            }' \
            | bgzip -c > ~{prefix}.vcf.gz

        tabix -p vcf ~{prefix}.vcf.gz

        rm -f tr_header_added.vcf.gz
    >>>

    output {
        File tagged_vcf = "~{prefix}.vcf.gz"
        File tagged_vcf_idx = "~{prefix}.vcf.gz.tbi"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 4,
        disk_gb: 3 * ceil(size(vcf, "GB")) + 5,
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

task DeduplicateOverlappingVariants {
    input {
        File vcf
        File vcf_idx
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        bcftools query \
            -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\n' \
            ~{vcf} \
            | awk 'BEGIN{OFS="\t"} {
                end = $2 + length($4)
                ref_len = length($4)
                print $1, $2, end, $3, ref_len
            }' > variants.bed

        bedtools intersect \
            -wa \
            -wb \
            -a variants.bed \
            -b variants.bed \
            | awk 'BEGIN{OFS="\t"} $4 != $9 {
                if ($5 < $10) {
                    print $4
                } else if ($5 == $10 && $4 > $9) {
                    print $4
                }
            }' \
            | sort -u > ids_to_remove.txt

        if [ -s ids_to_remove.txt ]; then
            bcftools view \
                -e "ID=@ids_to_remove.txt" \
                -Oz -o ~{prefix}.vcf.gz \
                ~{vcf}
        else
            cp ~{vcf} ~{prefix}.vcf.gz
        fi

        tabix -p vcf ~{prefix}.vcf.gz
    >>>

    output {
        File dedup_vcf = "~{prefix}.vcf.gz"
        File dedup_vcf_idx = "~{prefix}.vcf.gz.tbi"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 4,
        disk_gb: 3 * ceil(size(vcf, "GB")) + 10,
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

task PriorityMergeTRVcfs {
    input {
        Array[File] vcfs
        Array[File] vcf_idxs
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        VCF_FILES=(~{sep=' ' vcfs})

        cp "${VCF_FILES[0]}" merged.vcf.gz

        tabix -p vcf merged.vcf.gz

        for vcf_file in "${VCF_FILES[@]:1}"; do
            # BED of variants currently in the merged set
            bcftools query \
                -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\n' \
                merged.vcf.gz \
                | awk 'BEGIN{OFS="\t"} {print $1, $2, $2+length($4)}' > merged.bed

            # BED of candidate variants
            bcftools query \
                -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\n' \
                "$vcf_file" \
                | awk 'BEGIN{OFS="\t"} {print $1, $2, $2+length($4), $3}' > candidates.bed

            # Keep only candidates with no overlap in the merged set
            bedtools intersect \
                -v \
                -a candidates.bed \
                -b merged.bed \
                | cut -f4 > ids_to_add.txt

            if [ -s ids_to_add.txt ]; then
                bcftools view \
                    -i "ID=@ids_to_add.txt" \
                    -Oz -o to_add.vcf.gz \
                    "$vcf_file"

                tabix -p vcf to_add.vcf.gz

                bcftools concat \
                    --allow-overlaps \
                    merged.vcf.gz \
                    to_add.vcf.gz \
                    | bcftools sort -Oz -o new_merged.vcf.gz

                tabix -p vcf new_merged.vcf.gz

                mv new_merged.vcf.gz merged.vcf.gz
                mv new_merged.vcf.gz.tbi merged.vcf.gz.tbi

                rm -f to_add.vcf.gz to_add.vcf.gz.tbi
            fi

            rm -f merged.bed candidates.bed ids_to_add.txt
        done

        mv merged.vcf.gz ~{prefix}.vcf.gz
        mv merged.vcf.gz.tbi ~{prefix}.vcf.gz.tbi
    >>>

    output {
        File merged_vcf = "~{prefix}.vcf.gz"
        File merged_vcf_idx = "~{prefix}.vcf.gz.tbi"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 4,
        disk_gb: 5 * ceil(size(vcfs, "GB")) + 10,
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
