version 1.0

import "../utils/Helpers.wdl"
import "../utils/Structs.wdl"

workflow MergeVcfs {
    input {
        Array[File] contig_vcfs
        Array[File] contig_vcf_idxs
        String contig
        String prefix

        Int min_truvari_match = 20
        Int truvari_breakpoint_window = 500
        Float truvari_reciprocal_overlap = 0.0
        Float truvari_sample_similarity = 0.0
        Float truvari_sequence_similarity = 0.7
        Float truvari_size_similarity = 0.7
        Int truvari_size_max = 50000
        Int truvari_size_min = 20

        File ref_fa
        File ref_fai

        Int? shard_bin_size

        String utils_docker

        RuntimeAttr? runtime_attr_create_shards
        RuntimeAttr? runtime_attr_subset_shard
        RuntimeAttr? runtime_attr_subset_by_type
        RuntimeAttr? runtime_attr_merge_trv
        RuntimeAttr? runtime_attr_merge_non_trv
        RuntimeAttr? runtime_attr_consolidate_non_trv
        RuntimeAttr? runtime_attr_finalize_non_trv
        RuntimeAttr? runtime_attr_concat_merged
        RuntimeAttr? runtime_attr_concat_shards
        RuntimeAttr? runtime_attr_concat_shard_summaries
    }

    if (defined(shard_bin_size)) {
        call Helpers.CreateContigShards {
            input:
                vcfs = contig_vcfs,
                vcf_idxs = contig_vcf_idxs,
                contig = contig,
                shard_bin_size = select_first([shard_bin_size]),
                prefix = "~{prefix}.shards",
                docker = utils_docker,
                runtime_attr_override = runtime_attr_create_shards
        }

        scatter (j in range(length(CreateContigShards.shard_regions))) {
            scatter (i in range(length(contig_vcfs))) {
                call Helpers.SubsetVcfToRegion as SubsetCallsetToShard {
                    input:
                        vcf = contig_vcfs[i],
                        vcf_idx = contig_vcf_idxs[i],
                        region = CreateContigShards.shard_regions[j],
                        prefix = "~{prefix}.shard_~{j}.callset_~{i}",
                        docker = utils_docker,
                        runtime_attr_override = runtime_attr_subset_shard
                }

                call Helpers.SubsetVcfByArgs as SubsetTrvShard {
                    input:
                        vcf = SubsetCallsetToShard.subset_vcf,
                        vcf_idx = SubsetCallsetToShard.subset_vcf_idx,
                        include_args = 'INFO/allele_type=="trv"',
                        prefix = "~{prefix}.shard_~{j}.callset_~{i}.trv",
                        docker = utils_docker,
                        runtime_attr_override = runtime_attr_subset_by_type
                }

                call Helpers.SubsetVcfByArgs as SubsetNonTrvShard {
                    input:
                        vcf = SubsetCallsetToShard.subset_vcf,
                        vcf_idx = SubsetCallsetToShard.subset_vcf_idx,
                        include_args = 'INFO/allele_type!="trv"',
                        prefix = "~{prefix}.shard_~{j}.callset_~{i}.non_trv",
                        docker = utils_docker,
                        runtime_attr_override = runtime_attr_subset_by_type
                }
            }

            call MergeTrvVcfs as MergeTrvShard {
                input:
                    vcfs = SubsetTrvShard.subset_vcf,
                    vcf_idxs = SubsetTrvShard.subset_vcf_idx,
                    ref_fa = ref_fa,
                    ref_fai = ref_fai,
                    prefix = "~{prefix}.shard_~{j}.trv_merged",
                    docker = utils_docker,
                    runtime_attr_override = runtime_attr_merge_trv
            }

            call MergeNonTrvVcfs as MergeNonTrvShard {
                input:
                    vcfs = SubsetNonTrvShard.subset_vcf,
                    vcf_idxs = SubsetNonTrvShard.subset_vcf_idx,
                    ref_fa = ref_fa,
                    ref_fai = ref_fai,
                    min_truvari_match = min_truvari_match,
                    prefix = "~{prefix}.shard_~{j}.non_trv",
                    docker = utils_docker,
                    runtime_attr_override = runtime_attr_merge_non_trv
            }

            call Helpers.ConsolidateCollapsedSites as ConsolidateNonTrvShard {
                input:
                    vcf = MergeNonTrvShard.unmatched_large_vcf,
                    vcf_idx = MergeNonTrvShard.unmatched_large_vcf_idx,
                    breakpoint_window =truvari_breakpoint_window,
                    reciprocal_overlap = truvari_reciprocal_overlap,
                    sample_similarity = truvari_sample_similarity,
                    sequence_similarity = truvari_sequence_similarity,
                    size_similarity = truvari_size_similarity,
                    size_min = truvari_size_min,
                    size_max = truvari_size_max,

                    keep_strategy = "first",
                    set_merge_annotations = true,
                    strip_format_to_gt = true,
                    prefix = "~{prefix}.shard_~{j}.non_trv.truvari",
                    docker = utils_docker,
                    runtime_attr_override = runtime_attr_consolidate_non_trv
            }

            call FinalizeNonTrvMerge as FinalizeNonTrvShard {
                input:
                    matched_vcf = MergeNonTrvShard.matched_vcf,
                    matched_vcf_idx = MergeNonTrvShard.matched_vcf_idx,
                    unmatched_small_vcf = MergeNonTrvShard.unmatched_small_vcf,
                    unmatched_small_vcf_idx = MergeNonTrvShard.unmatched_small_vcf_idx,
                    consolidated_large_vcf = ConsolidateNonTrvShard.consolidated_vcf,
                    consolidated_large_vcf_idx = ConsolidateNonTrvShard.consolidated_vcf_idx,
                    n_non_trv_input = MergeNonTrvShard.n_non_trv_input,
                    n_truvari_input = MergeNonTrvShard.n_truvari_input,
                    contig = contig,
                    prefix = "~{prefix}.shard_~{j}.non_trv_merged",
                    docker = utils_docker,
                    runtime_attr_override = runtime_attr_finalize_non_trv
            }

            call Helpers.ConcatVcfs as ConcatMergedTypesShard {
                input:
                    vcfs = [MergeTrvShard.merged_vcf, FinalizeNonTrvShard.merged_vcf],
                    vcf_idxs = [MergeTrvShard.merged_vcf_idx, FinalizeNonTrvShard.merged_vcf_idx],
                    allow_overlaps = true,
                    naive = false,
                    prefix = "~{prefix}.shard_~{j}.merged",
                    docker = utils_docker,
                    runtime_attr_override = runtime_attr_concat_merged
            }
        }

        call Helpers.ConcatVcfs as ConcatShards {
            input:
                vcfs = ConcatMergedTypesShard.concat_vcf,
                vcf_idxs = ConcatMergedTypesShard.concat_vcf_idx,
                allow_overlaps = false,
                naive = true,
                prefix = "~{prefix}.merged",
                docker = utils_docker,
                runtime_attr_override = runtime_attr_concat_shards
        }

        call Helpers.ConcatTsvs as ConcatShardSummaries {
            input:
                tsvs = FinalizeNonTrvShard.summary_tsv,
                sort_output = false,
                preserve_header = true,
                prefix = "~{prefix}.merge_summary",
                docker = utils_docker,
                runtime_attr_override = runtime_attr_concat_shard_summaries
        }
    }

    if (!defined(shard_bin_size)) {
        scatter (i in range(length(contig_vcfs))) {
            call Helpers.SubsetVcfByArgs as SubsetTrv {
                input:
                    vcf = contig_vcfs[i],
                    vcf_idx = contig_vcf_idxs[i],
                    include_args = 'INFO/allele_type=="trv"',
                    prefix = "~{prefix}.callset_~{i}.trv",
                    docker = utils_docker,
                    runtime_attr_override = runtime_attr_subset_by_type
            }

            call Helpers.SubsetVcfByArgs as SubsetNonTrv {
                input:
                    vcf = contig_vcfs[i],
                    vcf_idx = contig_vcf_idxs[i],
                    include_args = 'INFO/allele_type!="trv"',
                    prefix = "~{prefix}.callset_~{i}.non_trv",
                    docker = utils_docker,
                    runtime_attr_override = runtime_attr_subset_by_type
            }
        }

        call MergeTrvVcfs {
            input:
                vcfs = SubsetTrv.subset_vcf,
                vcf_idxs = SubsetTrv.subset_vcf_idx,
                ref_fa = ref_fa,
                ref_fai = ref_fai,
                prefix = "~{prefix}.trv_merged",
                docker = utils_docker,
                runtime_attr_override = runtime_attr_merge_trv
        }

        call MergeNonTrvVcfs {
            input:
                vcfs = SubsetNonTrv.subset_vcf,
                vcf_idxs = SubsetNonTrv.subset_vcf_idx,
                ref_fa = ref_fa,
                ref_fai = ref_fai,
                min_truvari_match = min_truvari_match,
                prefix = "~{prefix}.non_trv",
                docker = utils_docker,
                runtime_attr_override = runtime_attr_merge_non_trv
        }

        call Helpers.ConsolidateCollapsedSites as ConsolidateNonTrv {
            input:
                vcf = MergeNonTrvVcfs.unmatched_large_vcf,
                vcf_idx = MergeNonTrvVcfs.unmatched_large_vcf_idx,
                breakpoint_window =truvari_breakpoint_window,
                reciprocal_overlap = truvari_reciprocal_overlap,
                sample_similarity = truvari_sample_similarity,
                sequence_similarity = truvari_sequence_similarity,
                size_similarity = truvari_size_similarity,
                size_min = truvari_size_min,
                size_max = truvari_size_max,
                keep_strategy = "first",
                set_merge_annotations = true,
                strip_format_to_gt = true,
                prefix = "~{prefix}.non_trv.truvari",
                docker = utils_docker,
                runtime_attr_override = runtime_attr_consolidate_non_trv
        }

        call FinalizeNonTrvMerge {
            input:
                matched_vcf = MergeNonTrvVcfs.matched_vcf,
                matched_vcf_idx = MergeNonTrvVcfs.matched_vcf_idx,
                unmatched_small_vcf = MergeNonTrvVcfs.unmatched_small_vcf,
                unmatched_small_vcf_idx = MergeNonTrvVcfs.unmatched_small_vcf_idx,
                consolidated_large_vcf = ConsolidateNonTrv.consolidated_vcf,
                consolidated_large_vcf_idx = ConsolidateNonTrv.consolidated_vcf_idx,
                n_non_trv_input = MergeNonTrvVcfs.n_non_trv_input,
                n_truvari_input = MergeNonTrvVcfs.n_truvari_input,
                contig = contig,
                prefix = "~{prefix}.non_trv_merged",
                docker = utils_docker,
                runtime_attr_override = runtime_attr_finalize_non_trv
        }

        call Helpers.ConcatVcfs as ConcatMergedTypes {
            input:
                vcfs = [MergeTrvVcfs.merged_vcf, FinalizeNonTrvMerge.merged_vcf],
                vcf_idxs = [MergeTrvVcfs.merged_vcf_idx, FinalizeNonTrvMerge.merged_vcf_idx],
                allow_overlaps = true,
                naive = false,
                prefix = "~{prefix}.merged",
                docker = utils_docker,
                runtime_attr_override = runtime_attr_concat_merged
        }
    }

    output {
        File merged_vcf = select_first([ConcatShards.concat_vcf, ConcatMergedTypes.concat_vcf])
        File merged_vcf_idx = select_first([ConcatShards.concat_vcf_idx, ConcatMergedTypes.concat_vcf_idx])
        File merge_summary_tsv = select_first([ConcatShardSummaries.concatenated_tsv, FinalizeNonTrvMerge.summary_tsv])
    }
}

task MergeTrvVcfs {
    input {
        Array[File] vcfs
        Array[File] vcf_idxs
        File ref_fa
        File ref_fai
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        echo '##INFO=<ID=MERGE_COUNT,Number=1,Type=Integer,Description="Number of source VCFs containing this variant">' > mc_hdr.txt
        echo '##INFO=<ID=MERGE_TYPE,Number=1,Type=String,Description="Merge strategy: EXACT, TRV_EXACT, TRUVARI, or UNIQUE">' > mt_hdr.txt

        paste ~{write_lines(vcfs)} ~{write_lines(vcf_idxs)} > vcf_pairs.tsv

        i=0
        while IFS=$'\t' read -r vcf vcf_idx; do
            if [[ "$vcf_idx" != "${vcf}.tbi" ]]; then
                ln -sf "$vcf_idx" "${vcf}.tbi"
            fi

            bcftools annotate \
                -x FORMAT/AL \
                -Oz -o "stripped_${i}.vcf.gz" \
                "$vcf"

            tabix -f -p vcf "stripped_${i}.vcf.gz"

            bcftools norm \
                --check-ref s \
                --do-not-normalize \
                -f ~{ref_fa} \
                -Oz -o "cleaned_${i}.vcf.gz" \
                "stripped_${i}.vcf.gz"

            tabix -f -p vcf "cleaned_${i}.vcf.gz"

            # Tag with MERGE_COUNT=1
            bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t1\n' "cleaned_${i}.vcf.gz" | bgzip > "mc_${i}.tsv.gz"
            tabix -s1 -b2 -e2 "mc_${i}.tsv.gz"
            bcftools annotate \
                -a "mc_${i}.tsv.gz" -h mc_hdr.txt \
                -c CHROM,POS,REF,ALT,MERGE_COUNT \
                -Oz -o "tagged_${i}.vcf.gz" "cleaned_${i}.vcf.gz"
            tabix -f -p vcf "tagged_${i}.vcf.gz"

            rm -f \
                "stripped_${i}.vcf.gz" "stripped_${i}.vcf.gz.tbi" \
                "cleaned_${i}.vcf.gz" "cleaned_${i}.vcf.gz.tbi" \
                "mc_${i}.tsv.gz" "mc_${i}.tsv.gz.tbi"

            echo "tagged_${i}.vcf.gz" >> tagged_vcfs.list

            i=$((i + 1))
        done < vcf_pairs.tsv

        bcftools merge \
            -m all \
            -i MERGE_COUNT:sum \
            -Oz -o merged.raw.vcf.gz \
            -l tagged_vcfs.list

        while read -r tagged_vcf; do
            rm -f "$tagged_vcf" "$tagged_vcf.tbi"
        done < tagged_vcfs.list

        # replace ';' with '_' so each merged site keeps a parseable ID
        bcftools view merged.raw.vcf.gz \
            | awk 'BEGIN{OFS="\t"} /^#/{print; next} {gsub(/;/,"_",$3); print}' \
            | bgzip > merged.unsorted.vcf.gz

        rm -f merged.raw.vcf.gz

        bcftools sort \
            --max-mem ~{select_first([runtime_attr.mem_gb, default_attr.mem_gb]) - 1}G \
            -T . \
            -Oz -o sorted.vcf.gz \
            merged.unsorted.vcf.gz
        
        tabix -f -p vcf sorted.vcf.gz

        rm -f merged.unsorted.vcf.gz

        # Annotate MERGE_TYPE based on MERGE_COUNT
        bcftools query \
            -f '%CHROM\t%POS\t%REF\t%ALT\t%INFO/MERGE_COUNT\n' \
            sorted.vcf.gz | \
            awk -F'\t' -v OFS='\t' '{if($5>1) $5="TRV_EXACT"; else $5="UNIQUE"; print}' | \
            bgzip > mt_annot.tsv.gz

        tabix -s1 -b2 -e2 mt_annot.tsv.gz

        bcftools annotate \
            -a mt_annot.tsv.gz -h mt_hdr.txt \
            -c CHROM,POS,REF,ALT,MERGE_TYPE \
            -Oz -o ~{prefix}.vcf.gz sorted.vcf.gz

        rm -f sorted.vcf.gz sorted.vcf.gz.tbi mt_annot.tsv.gz mt_annot.tsv.gz.tbi

        tabix -f -p vcf ~{prefix}.vcf.gz
    >>>

    output {
        File merged_vcf = "~{prefix}.vcf.gz"
        File merged_vcf_idx = "~{prefix}.vcf.gz.tbi"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 8,
        disk_gb: 5 * ceil(size(vcfs, "GB")) + 25,
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

task MergeNonTrvVcfs {
    input {
        Array[File] vcfs
        Array[File] vcf_idxs
        File ref_fa
        File ref_fai
        Int min_truvari_match
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        echo '##INFO=<ID=MERGE_COUNT,Number=1,Type=Integer,Description="Number of source VCFs containing this variant">' > mc_hdr.txt
        echo '##INFO=<ID=MERGE_TYPE,Number=1,Type=String,Description="Merge strategy: EXACT, TRV_EXACT, TRUVARI, or UNIQUE">' > mt_hdr.txt

        paste ~{write_lines(vcfs)} ~{write_lines(vcf_idxs)} > vcf_pairs.tsv

        i=0
        while IFS=$'\t' read -r vcf vcf_idx; do
            if [[ "$vcf_idx" != "${vcf}.tbi" ]]; then
                ln -sf "$vcf_idx" "${vcf}.tbi"
            fi

            bcftools annotate \
                -x FORMAT/AL \
                -Oz -o "stripped_${i}.vcf.gz" \
                "$vcf"

            tabix -f -p vcf "stripped_${i}.vcf.gz"

            bcftools norm \
                --check-ref s \
                --do-not-normalize \
                -f ~{ref_fa} \
                -Oz -o "cleaned_${i}.vcf.gz" \
                "stripped_${i}.vcf.gz"

            tabix -f -p vcf "cleaned_${i}.vcf.gz"

            # Tag with MERGE_COUNT=1
            bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t1\n' "cleaned_${i}.vcf.gz" | bgzip > "mc_${i}.tsv.gz"
            tabix -s1 -b2 -e2 "mc_${i}.tsv.gz"
            bcftools annotate \
                -a "mc_${i}.tsv.gz" -h mc_hdr.txt \
                -c CHROM,POS,REF,ALT,MERGE_COUNT \
                -Oz -o "tagged_${i}.vcf.gz" "cleaned_${i}.vcf.gz"
            tabix -f -p vcf "tagged_${i}.vcf.gz"

            rm -f \
                "stripped_${i}.vcf.gz" "stripped_${i}.vcf.gz.tbi" \
                "cleaned_${i}.vcf.gz" "cleaned_${i}.vcf.gz.tbi" \
                "mc_${i}.tsv.gz" "mc_${i}.tsv.gz.tbi"

            echo "tagged_${i}.vcf.gz" >> tagged_vcfs.list

            i=$((i + 1))
        done < vcf_pairs.tsv

        bcftools merge \
            -m none \
            -i MERGE_COUNT:sum \
            -Oz -o exact_merged.raw.vcf.gz \
            -l tagged_vcfs.list

        while read -r tagged_vcf; do
            rm -f "$tagged_vcf" "$tagged_vcf.tbi"
        done < tagged_vcfs.list

        # Joined IDs with ';' are rewritten to '_' so the ID column stays parseable
        bcftools view exact_merged.raw.vcf.gz \
            | awk 'BEGIN{OFS="\t"} /^#/{print; next} {gsub(/;/,"_",$3); print}' \
            | bgzip > exact_merged.unsorted.vcf.gz

        rm -f exact_merged.raw.vcf.gz

        bcftools sort \
            --max-mem ~{select_first([runtime_attr.mem_gb, default_attr.mem_gb]) - 1}G \
            -T . \
            -Oz -o exact_merged.vcf.gz \
            exact_merged.unsorted.vcf.gz

        tabix -f -p vcf exact_merged.vcf.gz

        rm -f exact_merged.unsorted.vcf.gz

        bcftools view -H exact_merged.vcf.gz | wc -l | awk '{print $1}' > n_non_trv_input.txt

        # Split into matched (exact), unmatched large (truvari candidates), and unmatched small (unique)
        bcftools view \
            -i 'MERGE_COUNT>1' \
            -Oz -o matched.vcf.gz \
            exact_merged.vcf.gz

        bcftools view \
            -i 'MERGE_COUNT=1 && abs(INFO/allele_length) >= ~{min_truvari_match}' \
            -Oz -o ~{prefix}.unmatched_large.vcf.gz \
            exact_merged.vcf.gz

        bcftools view \
            -i 'MERGE_COUNT=1 && abs(INFO/allele_length) < ~{min_truvari_match}' \
            -Oz -o unmatched_small.vcf.gz \
            exact_merged.vcf.gz

        tabix -f -p vcf matched.vcf.gz
        tabix -f -p vcf ~{prefix}.unmatched_large.vcf.gz
        tabix -f -p vcf unmatched_small.vcf.gz

        rm -f exact_merged.vcf.gz exact_merged.vcf.gz.tbi

        bcftools view -H ~{prefix}.unmatched_large.vcf.gz | wc -l | awk '{print $1}' > n_truvari_input.txt

        # Annotate matched with MERGE_TYPE=EXACT
        bcftools query \
            -f '%CHROM\t%POS\t%REF\t%ALT\tEXACT\n' \
            matched.vcf.gz \
            | bgzip > matched_mt.tsv.gz

        tabix -s1 -b2 -e2 matched_mt.tsv.gz

        bcftools annotate \
            -a matched_mt.tsv.gz \
            -h mt_hdr.txt \
            -c CHROM,POS,REF,ALT,MERGE_TYPE \
            -Oz -o ~{prefix}.matched.vcf.gz \
            matched.vcf.gz

        tabix -f -p vcf ~{prefix}.matched.vcf.gz

        rm -f matched.vcf.gz matched.vcf.gz.tbi matched_mt.tsv.gz matched_mt.tsv.gz.tbi

        # Annotate unmatched_small with MERGE_TYPE=UNIQUE
        bcftools query \
            -f '%CHROM\t%POS\t%REF\t%ALT\tUNIQUE\n' \
            unmatched_small.vcf.gz \
            | bgzip > small_mt.tsv.gz

        tabix -s1 -b2 -e2 small_mt.tsv.gz

        bcftools annotate \
            -a small_mt.tsv.gz \
            -h mt_hdr.txt \
            -c CHROM,POS,REF,ALT,MERGE_TYPE \
            -Oz -o ~{prefix}.unmatched_small.vcf.gz \
            unmatched_small.vcf.gz

        tabix -f -p vcf ~{prefix}.unmatched_small.vcf.gz

        rm -f unmatched_small.vcf.gz unmatched_small.vcf.gz.tbi small_mt.tsv.gz small_mt.tsv.gz.tbi
    >>>

    output {
        File matched_vcf = "~{prefix}.matched.vcf.gz"
        File matched_vcf_idx = "~{prefix}.matched.vcf.gz.tbi"
        File unmatched_large_vcf = "~{prefix}.unmatched_large.vcf.gz"
        File unmatched_large_vcf_idx = "~{prefix}.unmatched_large.vcf.gz.tbi"
        File unmatched_small_vcf = "~{prefix}.unmatched_small.vcf.gz"
        File unmatched_small_vcf_idx = "~{prefix}.unmatched_small.vcf.gz.tbi"
        Int n_non_trv_input = read_int("n_non_trv_input.txt")
        Int n_truvari_input = read_int("n_truvari_input.txt")
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 2,
        mem_gb: 16,
        disk_gb: 5 * ceil(size(vcfs, "GB")) + 25,
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

task FinalizeNonTrvMerge {
    input {
        File matched_vcf
        File matched_vcf_idx
        File unmatched_small_vcf
        File unmatched_small_vcf_idx
        File consolidated_large_vcf
        File consolidated_large_vcf_idx
        Int n_non_trv_input
        Int n_truvari_input
        String contig
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        bcftools concat \
            -a \
            -Oz -o concat.unsorted.vcf.gz \
            ~{matched_vcf} ~{unmatched_small_vcf} ~{consolidated_large_vcf}

        bcftools sort \
            --max-mem ~{select_first([runtime_attr.mem_gb, default_attr.mem_gb]) - 1}G \
            -T . \
            -Oz -o sorted.vcf.gz \
            concat.unsorted.vcf.gz

        rm -f concat.unsorted.vcf.gz

        bcftools +fill-tags \
            sorted.vcf.gz \
            -Oz -o ~{prefix}.vcf.gz \
            -- -t AC,AN,AF

        rm -f sorted.vcf.gz

        tabix -f -p vcf ~{prefix}.vcf.gz

        n_output=$(bcftools view -H ~{prefix}.vcf.gz | wc -l | awk '{print $1}')
        n_truvari_output=$(bcftools view -H ~{consolidated_large_vcf} | wc -l | awk '{print $1}')
        n_truvari_collapsed=$((~{n_truvari_input} - n_truvari_output))

        cat > ~{prefix}.merge_summary.tsv <<EOF
contig	n_non_trv_input	n_non_trv_output	n_truvari_input	n_truvari_output	n_truvari_collapsed
~{contig}	~{n_non_trv_input}	${n_output}	~{n_truvari_input}	${n_truvari_output}	${n_truvari_collapsed}
EOF
    >>>

    output {
        File merged_vcf = "~{prefix}.vcf.gz"
        File merged_vcf_idx = "~{prefix}.vcf.gz.tbi"
        File summary_tsv = "~{prefix}.merge_summary.tsv"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 8,
        disk_gb: 3 * ceil(size(matched_vcf, "GB") + size(unmatched_small_vcf, "GB") + size(consolidated_large_vcf, "GB")) + 25,
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
