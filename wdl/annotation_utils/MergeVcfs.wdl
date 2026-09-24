version 1.0

import "../utils/Helpers.wdl"
import "../utils/Structs.wdl"

workflow MergeVcfs {
    meta {
        description: [
            "This utility merges per-contig VCFs that were called across distinct sample sets into a single callset for the contig, handling tandem-repeat and non-tandem-repeat variants separately.",
            "Tandem-repeat variants are merged on CHROM, POS and REF alone, so records that describe the same locus with different ALT alleles collapse into one multiallelic record whose ALT list is the union of the inputs and whose genotypes are remapped onto it.",
            "Non-tandem-repeat variants are first merged on an exact CHROM, POS, REF and ALT match. Records that stay unmatched and are at least `min_truvari_match` long are then collapsed with Truvari (https://github.com/ACEnglish/truvari) using a breakpoint distance, reciprocal overlap, and sequence, size and sample similarity; shorter unmatched records pass through untouched.",
            "The contig is split into bins of `shard_bin_size` and every merging step runs per shard, so Truvari never pairs records more than one bin apart. Both merged and unmerged records reach the output.",
            "Every output record carries `MERGE_COUNT`, the number of input records merged into it, `MERGE_TYPE`, one of EXACT, TRV_EXACT, TRUVARI or UNIQUE, and `MERGE_SOURCE`, which is MERGED for records drawn from more than one input VCF and otherwise the `vcf_names` entry of the single VCF that carried it.",
            "Where a merged record cannot hold both inputs' values, the ID and any INFO field other than `MERGE_COUNT` are taken from the first input VCF that carried the record, and AC, AN and AF are recomputed over the merged samples."
        ]
    }

    parameter_meta {
        contig_vcfs: "Per-callset VCFs for the contig being merged, each called across a distinct set of samples."
        contig_vcf_idxs: "Indexes for `contig_vcfs`."
        vcf_names: "Name of each entry of `contig_vcfs`, in the same order, used as the `MERGE_SOURCE` value for records that only one callset carried."
        contig: "Contig being merged."
        min_truvari_match: "Minimum variant length for Truvari matching."
        truvari_breakpoint_window: "Maximum breakpoint distance, in bp, for merging non-TR variants."
        truvari_reciprocal_overlap: "Minimum reciprocal overlap for merging non-TR variants."
        truvari_sample_similarity: "Minimum sample similarity for merging non-TR variants."
        truvari_sequence_similarity: "Minimum sequence similarity for merging non-TR variants."
        truvari_size_similarity: "Minimum size similarity for merging non-TR variants."
        truvari_size_max: "Maximum variant length Truvari will consider when collapsing."
        truvari_size_min: "Minimum variant length Truvari will consider when collapsing."
        ref_fa: "From references."
        ref_fai: "From references."
        shard_bin_size: "Region-bin size, in bp, used when sharding the contig."
        merged_vcf: "Merged VCF."
        merged_vcf_idx: "Index for the merged VCF."
        merge_summary_tsv: "TSV summarizing the merge."
    }

    input {
        Array[File] contig_vcfs
        Array[File] contig_vcf_idxs
        Array[String] vcf_names
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

        Int shard_bin_size = 10000000

        String utils_docker

        RuntimeAttr? runtime_attr_create_shards
        RuntimeAttr? runtime_attr_select_shards
        RuntimeAttr? runtime_attr_subset_shard
        RuntimeAttr? runtime_attr_merge_trv
        RuntimeAttr? runtime_attr_merge_non_trv
        RuntimeAttr? runtime_attr_consolidate_non_trv
        RuntimeAttr? runtime_attr_finalize_non_trv
        RuntimeAttr? runtime_attr_concat_merged
        RuntimeAttr? runtime_attr_concat_shards
        RuntimeAttr? runtime_attr_concat_shard_summaries
    }

    call Helpers.CreateContigShards {
        input:
            vcfs = contig_vcfs,
            vcf_idxs = contig_vcf_idxs,
            contig = contig,
            shard_bin_size = shard_bin_size,
            prefix = "~{prefix}.shards",
            docker = utils_docker,
            runtime_attr_override = runtime_attr_create_shards
    }

    call SelectPopulatedShards {
        input:
            vcfs = contig_vcfs,
            vcf_idxs = contig_vcf_idxs,
            regions = CreateContigShards.shard_regions,
            prefix = "~{prefix}.populated_shards",
            docker = utils_docker,
            runtime_attr_override = runtime_attr_select_shards
    }

    scatter (j in range(length(SelectPopulatedShards.populated_regions))) {
        scatter (i in range(length(contig_vcfs))) {
            call Helpers.SubsetVcfToRegion as SubsetCallsetToShard {
                input:
                    vcf = contig_vcfs[i],
                    vcf_idx = contig_vcf_idxs[i],
                    region = SelectPopulatedShards.populated_regions[j],
                    prefix = "~{prefix}.shard_~{j}.callset_~{i}",
                    docker = utils_docker,
                    runtime_attr_override = runtime_attr_subset_shard
            }
        }

        call MergeTrvVcfs {
            input:
                vcfs = SubsetCallsetToShard.subset_vcf,
                vcf_idxs = SubsetCallsetToShard.subset_vcf_idx,
                vcf_names = vcf_names,
                ref_fa = ref_fa,
                ref_fai = ref_fai,
                prefix = "~{prefix}.shard_~{j}.trv_merged",
                docker = utils_docker,
                runtime_attr_override = runtime_attr_merge_trv
        }

        call MergeNonTrvVcfs {
            input:
                vcfs = SubsetCallsetToShard.subset_vcf,
                vcf_idxs = SubsetCallsetToShard.subset_vcf_idx,
                vcf_names = vcf_names,
                ref_fa = ref_fa,
                ref_fai = ref_fai,
                min_truvari_match = min_truvari_match,
                prefix = "~{prefix}.shard_~{j}.non_trv",
                docker = utils_docker,
                runtime_attr_override = runtime_attr_merge_non_trv
        }

        call Helpers.ConsolidateCollapsedSites as ConsolidateNonTrv {
            input:
                vcf = MergeNonTrvVcfs.unmatched_large_vcf,
                vcf_idx = MergeNonTrvVcfs.unmatched_large_vcf_idx,
                breakpoint_window = truvari_breakpoint_window,
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
                prefix = "~{prefix}.shard_~{j}.non_trv_merged",
                docker = utils_docker,
                runtime_attr_override = runtime_attr_finalize_non_trv
        }

        call Helpers.ConcatVcfs as ConcatMergedTypes {
            input:
                vcfs = [MergeTrvVcfs.merged_vcf, FinalizeNonTrvMerge.merged_vcf],
                vcf_idxs = [MergeTrvVcfs.merged_vcf_idx, FinalizeNonTrvMerge.merged_vcf_idx],
                allow_overlaps = false,
                naive = false,
                sort_output = true,
                prefix = "~{prefix}.shard_~{j}.merged",
                docker = utils_docker,
                runtime_attr_override = runtime_attr_concat_merged
        }
    }

    call Helpers.ConcatVcfs as ConcatShards {
        input:
            vcfs = ConcatMergedTypes.concat_vcf,
            vcf_idxs = ConcatMergedTypes.concat_vcf_idx,
            allow_overlaps = false,
            naive = false,
            prefix = "~{prefix}.merged",
            docker = utils_docker,
            runtime_attr_override = runtime_attr_concat_shards
    }

    call Helpers.ConcatTsvs as ConcatShardSummaries {
        input:
            tsvs = FinalizeNonTrvMerge.summary_tsv,
            sort_output = false,
            preserve_header = true,
            prefix = "~{prefix}.merge_summary",
            docker = utils_docker,
            runtime_attr_override = runtime_attr_concat_shard_summaries
    }

    output {
        File merged_vcf = ConcatShards.concat_vcf
        File merged_vcf_idx = ConcatShards.concat_vcf_idx
        File merge_summary_tsv = ConcatShardSummaries.concatenated_tsv
    }
}

task SelectPopulatedShards {
    input {
        Array[File] vcfs
        Array[File] vcf_idxs
        Array[String] regions
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        paste ~{write_lines(vcfs)} ~{write_lines(vcf_idxs)} > vcf_pairs.tsv

        while IFS=$'\t' read -r vcf vcf_idx; do
            if [[ "$vcf_idx" != "${vcf}.tbi" ]]; then
                ln -sf "$vcf_idx" "${vcf}.tbi"
            fi
        done < vcf_pairs.tsv

        # Keep only the regions that at least one callset has a record in, so empty shards cost nothing
        : > ~{prefix}.txt
        while read -r region; do
            while IFS=$'\t' read -r vcf vcf_idx; do
                if [[ -n "$(bcftools view -H -r "$region" "$vcf" | head -n 1)" ]]; then
                    echo "$region" >> ~{prefix}.txt
                    break
                fi
            done < vcf_pairs.tsv
        done < ~{write_lines(regions)}

        # Fall back to the first region so that a contig with no records still yields one empty shard
        if [[ ! -s ~{prefix}.txt ]]; then
            head -n 1 ~{write_lines(regions)} > ~{prefix}.txt
        fi
    >>>

    output {
        Array[String] populated_regions = read_lines("~{prefix}.txt")
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 4,
        disk_gb: 2 * ceil(size(vcfs, "GB")) + 10,
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

task MergeTrvVcfs {
    input {
        Array[File] vcfs
        Array[File] vcf_idxs
        Array[String] vcf_names
        File ref_fa
        File ref_fai
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        echo '##INFO=<ID=ORIG_ID,Number=1,Type=String,Description="Variant ID in the input callset">' > merge_hdr.txt
        echo '##INFO=<ID=ORIG_REF,Number=1,Type=String,Description="REF allele as written in the input callset">' >> merge_hdr.txt
        echo '##INFO=<ID=MERGE_COUNT,Number=1,Type=Integer,Description="Number of input records merged into this site">' >> merge_hdr.txt
        echo '##INFO=<ID=MERGE_TYPE,Number=1,Type=String,Description="Merge strategy: EXACT, TRV_EXACT, TRUVARI, or UNIQUE">' >> merge_hdr.txt
        echo '##INFO=<ID=MERGE_SOURCE,Number=1,Type=String,Description="MERGED when the record came from more than one input VCF, otherwise the name of the input VCF that carried it">' >> merge_hdr.txt

        paste ~{write_lines(vcfs)} ~{write_lines(vcf_idxs)} ~{write_lines(vcf_names)} > vcf_triples.tsv

        i=0
        while IFS=$'\t' read -r vcf vcf_idx vcf_name; do
            if [[ "$vcf_idx" != "${vcf}.tbi" ]]; then
                ln -sf "$vcf_idx" "${vcf}.tbi"
            fi

            bcftools view -i 'INFO/allele_type=="trv"' -Oz -o "typed_${i}.vcf.gz" "$vcf"
            tabix -f -p vcf "typed_${i}.vcf.gz"

            # Drop annotations left by an earlier Truvari run, so this output's header matches the non-TR one
            STALE_TAGS=$(bcftools view -h "typed_${i}.vcf.gz" | grep '^##INFO' \
                | sed 's/.*ID=\([^,]*\).*/\1/' \
                | grep -Ex 'NumCollapsed|NumConsolidated|CollapseId|MatchId|TruScore' \
                | sed 's|^|INFO/|' | paste -sd',' - || true)

            if [[ -n "$STALE_TAGS" ]]; then
                bcftools annotate -x "$STALE_TAGS" -Oz -o "stripped_${i}.vcf.gz" "typed_${i}.vcf.gz"
            else
                mv "typed_${i}.vcf.gz" "stripped_${i}.vcf.gz"
            fi

            tabix -f -p vcf "stripped_${i}.vcf.gz"

            bcftools norm \
                --check-ref s \
                --do-not-normalize \
                -f ~{ref_fa} \
                -Oz -o "cleaned_${i}.vcf.gz" \
                "stripped_${i}.vcf.gz"

            tabix -f -p vcf "cleaned_${i}.vcf.gz"

            # Stamp provenance, and stash the ID and REF that bcftools merge would otherwise rewrite
            bcftools query \
                -f "%CHROM\t%POS\t%REF\t%ALT\t%ID\t%REF\t1\t${vcf_name}\n" \
                "cleaned_${i}.vcf.gz" \
                | bgzip > "annot_${i}.tsv.gz"

            tabix -s1 -b2 -e2 "annot_${i}.tsv.gz"

            bcftools annotate \
                -a "annot_${i}.tsv.gz" -h merge_hdr.txt \
                -c CHROM,POS,REF,ALT,ORIG_ID,ORIG_REF,MERGE_COUNT,MERGE_SOURCE \
                -Oz -o "stamped_${i}.vcf.gz" \
                "cleaned_${i}.vcf.gz"

            tabix -f -p vcf "stamped_${i}.vcf.gz"

            # Set the ID to the REF so that 'bcftools merge -m id' keys on CHROM, POS and REF alone
            bcftools annotate --set-id '%REF' -Oz -o "tagged_${i}.vcf.gz" "stamped_${i}.vcf.gz"

            tabix -f -p vcf "tagged_${i}.vcf.gz"

            rm -f \
                "typed_${i}.vcf.gz" "typed_${i}.vcf.gz.tbi" \
                "stripped_${i}.vcf.gz" "stripped_${i}.vcf.gz.tbi" \
                "cleaned_${i}.vcf.gz" "cleaned_${i}.vcf.gz.tbi" \
                "stamped_${i}.vcf.gz" "stamped_${i}.vcf.gz.tbi" \
                "annot_${i}.tsv.gz" "annot_${i}.tsv.gz.tbi"

            echo "tagged_${i}.vcf.gz" >> tagged_vcfs.list

            i=$((i + 1))
        done < vcf_triples.tsv

        n_records=0
        while read -r tagged_vcf; do
            n_records=$((n_records + $(bcftools view -H "$tagged_vcf" | wc -l)))
        done < tagged_vcfs.list

        # bcftools merge crashes when every input is empty, so build the header on its own and stop
        if (( n_records == 0 )); then
            bcftools merge --print-header -m id -l tagged_vcfs.list \
                | bcftools annotate -x INFO/ORIG_ID,INFO/ORIG_REF -Oz -o ~{prefix}.vcf.gz
            tabix -f -p vcf ~{prefix}.vcf.gz
            exit 0
        fi

        bcftools merge \
            -m id \
            -i MERGE_COUNT:sum \
            -Ov -o merged.raw.vcf \
            -l tagged_vcfs.list

        while read -r tagged_vcf; do
            rm -f "$tagged_vcf" "$tagged_vcf.tbi"
        done < tagged_vcfs.list

        # bcftools merge trims the bases REF and ALT share, so restore the input representation from ORIG_REF
        awk -F'\t' -v OFS='\t' '
            /^#/ { print; next }
            {
                orig = ""
                if (match($8, /(^|;)ORIG_REF=[^;]*/)) {
                    orig = substr($8, RSTART, RLENGTH)
                    sub(/^;/, "", orig)
                    sub(/^ORIG_REF=/, "", orig)
                }
                if (orig != "" && length($4) < length(orig)) {
                    suffix = substr(orig, length($4) + 1)
                    $4 = orig
                    n = split($5, alts, ",")
                    rebuilt = alts[1] suffix
                    for (k = 2; k <= n; k++) {
                        rebuilt = rebuilt "," alts[k] suffix
                    }
                    $5 = rebuilt
                }
                print
            }
        ' merged.raw.vcf | bgzip > untrimmed.vcf.gz

        tabix -f -p vcf untrimmed.vcf.gz

        rm -f merged.raw.vcf

        # A site seen in more than one callset is a tandem-repeat match on CHROM, POS and REF
        bcftools query \
            -f '%CHROM\t%POS\t%REF\t%ALT\t%INFO/MERGE_COUNT\n' \
            untrimmed.vcf.gz \
            | awk -F'\t' -v OFS='\t' '{ if ($5 > 1) { $5 = "TRV_EXACT" } else { $5 = "UNIQUE" } ; print }' \
            | bgzip > merge_type.tsv.gz

        tabix -s1 -b2 -e2 merge_type.tsv.gz

        bcftools query \
            -i 'MERGE_COUNT>1' \
            -f '%CHROM\t%POS\t%REF\t%ALT\tMERGED\n' \
            untrimmed.vcf.gz \
            | bgzip > merge_source.tsv.gz

        tabix -s1 -b2 -e2 merge_source.tsv.gz

        bcftools annotate \
            -a merge_type.tsv.gz \
            -c CHROM,POS,REF,ALT,MERGE_TYPE \
            -Oz -o type_filled.vcf.gz \
            untrimmed.vcf.gz

        tabix -f -p vcf type_filled.vcf.gz

        bcftools annotate \
            -a merge_source.tsv.gz \
            -c CHROM,POS,REF,ALT,MERGE_SOURCE \
            -Oz -o source_filled.vcf.gz \
            type_filled.vcf.gz

        tabix -f -p vcf source_filled.vcf.gz

        bcftools annotate \
            --set-id '%INFO/ORIG_ID' \
            -Oz -o id_restored.vcf.gz \
            source_filled.vcf.gz

        tabix -f -p vcf id_restored.vcf.gz

        bcftools annotate \
            -x INFO/ORIG_ID,INFO/ORIG_REF \
            -Oz -o ~{prefix}.vcf.gz \
            id_restored.vcf.gz

        tabix -f -p vcf ~{prefix}.vcf.gz

        rm -f untrimmed.vcf.gz untrimmed.vcf.gz.tbi merge_type.tsv.gz* merge_source.tsv.gz* \
            type_filled.vcf.gz* source_filled.vcf.gz* id_restored.vcf.gz*
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
        Array[String] vcf_names
        File ref_fa
        File ref_fai
        Int min_truvari_match
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        echo '##INFO=<ID=ORIG_ID,Number=1,Type=String,Description="Variant ID in the input callset">' > merge_hdr.txt
        echo '##INFO=<ID=ORIG_REF,Number=1,Type=String,Description="REF allele as written in the input callset">' >> merge_hdr.txt
        echo '##INFO=<ID=MERGE_COUNT,Number=1,Type=Integer,Description="Number of input records merged into this site">' >> merge_hdr.txt
        echo '##INFO=<ID=MERGE_TYPE,Number=1,Type=String,Description="Merge strategy: EXACT, TRV_EXACT, TRUVARI, or UNIQUE">' >> merge_hdr.txt
        echo '##INFO=<ID=MERGE_SOURCE,Number=1,Type=String,Description="MERGED when the record came from more than one input VCF, otherwise the name of the input VCF that carried it">' >> merge_hdr.txt

        paste ~{write_lines(vcfs)} ~{write_lines(vcf_idxs)} ~{write_lines(vcf_names)} > vcf_triples.tsv

        i=0
        while IFS=$'\t' read -r vcf vcf_idx vcf_name; do
            if [[ "$vcf_idx" != "${vcf}.tbi" ]]; then
                ln -sf "$vcf_idx" "${vcf}.tbi"
            fi

            bcftools view -e 'INFO/allele_type=="trv"' -Oz -o "typed_${i}.vcf.gz" "$vcf"
            tabix -f -p vcf "typed_${i}.vcf.gz"

            # Drop annotations left by an earlier Truvari run, whose ids would collide with this one
            STALE_TAGS=$(bcftools view -h "typed_${i}.vcf.gz" | grep '^##INFO' \
                | sed 's/.*ID=\([^,]*\).*/\1/' \
                | grep -Ex 'NumCollapsed|NumConsolidated|CollapseId|MatchId|TruScore' \
                | sed 's|^|INFO/|' | paste -sd',' - || true)

            if [[ -n "$STALE_TAGS" ]]; then
                bcftools annotate -x "$STALE_TAGS" -Oz -o "stripped_${i}.vcf.gz" "typed_${i}.vcf.gz"
            else
                mv "typed_${i}.vcf.gz" "stripped_${i}.vcf.gz"
            fi

            tabix -f -p vcf "stripped_${i}.vcf.gz"

            bcftools norm \
                --check-ref s \
                --do-not-normalize \
                -f ~{ref_fa} \
                -Oz -o "cleaned_${i}.vcf.gz" \
                "stripped_${i}.vcf.gz"

            tabix -f -p vcf "cleaned_${i}.vcf.gz"

            # Stamp provenance, and stash the ID and REF that bcftools merge would otherwise rewrite
            bcftools query \
                -f "%CHROM\t%POS\t%REF\t%ALT\t%ID\t%REF\t1\t${vcf_name}\n" \
                "cleaned_${i}.vcf.gz" \
                | bgzip > "annot_${i}.tsv.gz"

            tabix -s1 -b2 -e2 "annot_${i}.tsv.gz"

            bcftools annotate \
                -a "annot_${i}.tsv.gz" -h merge_hdr.txt \
                -c CHROM,POS,REF,ALT,ORIG_ID,ORIG_REF,MERGE_COUNT,MERGE_SOURCE \
                -Oz -o "stamped_${i}.vcf.gz" \
                "cleaned_${i}.vcf.gz"

            tabix -f -p vcf "stamped_${i}.vcf.gz"

            # Qualify the ID with the callset name, since Truvari needs the IDs it collapses to be unique
            bcftools annotate --set-id "${vcf_name}:%ID" -Oz -o "tagged_${i}.vcf.gz" "stamped_${i}.vcf.gz"

            tabix -f -p vcf "tagged_${i}.vcf.gz"

            rm -f \
                "typed_${i}.vcf.gz" "typed_${i}.vcf.gz.tbi" \
                "stripped_${i}.vcf.gz" "stripped_${i}.vcf.gz.tbi" \
                "cleaned_${i}.vcf.gz" "cleaned_${i}.vcf.gz.tbi" \
                "stamped_${i}.vcf.gz" "stamped_${i}.vcf.gz.tbi" \
                "annot_${i}.tsv.gz" "annot_${i}.tsv.gz.tbi"

            echo "tagged_${i}.vcf.gz" >> tagged_vcfs.list

            i=$((i + 1))
        done < vcf_triples.tsv

        n_records=0
        while read -r tagged_vcf; do
            n_records=$((n_records + $(bcftools view -H "$tagged_vcf" | wc -l)))
        done < tagged_vcfs.list

        # bcftools merge crashes when every input is empty, so build the header on its own and stop
        if (( n_records == 0 )); then
            bcftools merge --print-header -m none -l tagged_vcfs.list \
                | bcftools annotate -x INFO/ORIG_REF -Oz -o empty.vcf.gz

            for split in matched unmatched_large unmatched_small; do
                cp empty.vcf.gz "~{prefix}.${split}.vcf.gz"
                tabix -f -p vcf "~{prefix}.${split}.vcf.gz"
            done

            echo 0 > n_non_trv_input.txt
            echo 0 > n_truvari_input.txt
            exit 0
        fi

        bcftools merge \
            -m none \
            -i MERGE_COUNT:sum \
            -Ov -o exact_merged.raw.vcf \
            -l tagged_vcfs.list

        while read -r tagged_vcf; do
            rm -f "$tagged_vcf" "$tagged_vcf.tbi"
        done < tagged_vcfs.list

        # bcftools merge trims the bases REF and ALT share, so restore the input representation from ORIG_REF
        awk -F'\t' -v OFS='\t' '
            /^#/ { print; next }
            {
                orig = ""
                if (match($8, /(^|;)ORIG_REF=[^;]*/)) {
                    orig = substr($8, RSTART, RLENGTH)
                    sub(/^;/, "", orig)
                    sub(/^ORIG_REF=/, "", orig)
                }
                if (orig != "" && length($4) < length(orig)) {
                    suffix = substr(orig, length($4) + 1)
                    $4 = orig
                    n = split($5, alts, ",")
                    rebuilt = alts[1] suffix
                    for (k = 2; k <= n; k++) {
                        rebuilt = rebuilt "," alts[k] suffix
                    }
                    $5 = rebuilt
                }
                print
            }
        ' exact_merged.raw.vcf \
            | bcftools annotate -x INFO/ORIG_REF -Oz -o exact_merged.vcf.gz

        tabix -f -p vcf exact_merged.vcf.gz

        rm -f exact_merged.raw.vcf

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

        # A plain concat followed by the sort below, since 'bcftools concat -a' crashes when every input is empty
        bcftools concat \
            -Oz -o concat.unsorted.vcf.gz \
            ~{matched_vcf} ~{unmatched_small_vcf} ~{consolidated_large_vcf}

        bcftools sort \
            --max-mem ~{select_first([runtime_attr.mem_gb, default_attr.mem_gb]) - 1}G \
            -T . \
            -Oz -o sorted.vcf.gz \
            concat.unsorted.vcf.gz

        tabix -f -p vcf sorted.vcf.gz

        rm -f concat.unsorted.vcf.gz

        # Truvari re-emits records it declined to collapse on sample similarity without a MERGE_TYPE
        bcftools query \
            -e 'INFO/MERGE_TYPE!="."' \
            -f '%CHROM\t%POS\t%REF\t%ALT\tUNIQUE\n' \
            sorted.vcf.gz \
            | bgzip > missing_mt.tsv.gz

        tabix -s1 -b2 -e2 missing_mt.tsv.gz

        bcftools query \
            -i 'MERGE_COUNT>1' \
            -f '%CHROM\t%POS\t%REF\t%ALT\tMERGED\n' \
            sorted.vcf.gz \
            | bgzip > merge_source.tsv.gz

        tabix -s1 -b2 -e2 merge_source.tsv.gz

        bcftools annotate \
            -a missing_mt.tsv.gz \
            -c CHROM,POS,REF,ALT,MERGE_TYPE \
            -Oz -o type_filled.vcf.gz \
            sorted.vcf.gz

        tabix -f -p vcf type_filled.vcf.gz

        bcftools annotate \
            -a merge_source.tsv.gz \
            -c CHROM,POS,REF,ALT,MERGE_SOURCE \
            -Oz -o source_filled.vcf.gz \
            type_filled.vcf.gz

        tabix -f -p vcf source_filled.vcf.gz

        bcftools annotate \
            --set-id '%INFO/ORIG_ID' \
            -Oz -o id_restored.vcf.gz \
            source_filled.vcf.gz

        tabix -f -p vcf id_restored.vcf.gz

        bcftools annotate \
            -x INFO/ORIG_ID \
            -Oz -o stripped.vcf.gz \
            id_restored.vcf.gz

        tabix -f -p vcf stripped.vcf.gz

        bcftools +fill-tags stripped.vcf.gz -Oz -o ~{prefix}.vcf.gz -- -t AC,AN,AF

        rm -f sorted.vcf.gz sorted.vcf.gz.tbi missing_mt.tsv.gz* merge_source.tsv.gz* \
            type_filled.vcf.gz* source_filled.vcf.gz* id_restored.vcf.gz* stripped.vcf.gz*

        tabix -f -p vcf ~{prefix}.vcf.gz

        n_output=$(bcftools view -H ~{prefix}.vcf.gz | wc -l | awk '{print $1}')
        n_truvari_output=$(bcftools view -H ~{consolidated_large_vcf} | wc -l | awk '{print $1}')
        n_truvari_collapsed=$((~{n_truvari_input} - n_truvari_output))

        printf '%s\t%s\t%s\t%s\t%s\t%s\n' \
            contig n_non_trv_input n_non_trv_output n_truvari_input n_truvari_output n_truvari_collapsed \
            > ~{prefix}.merge_summary.tsv

        printf '%s\t%s\t%s\t%s\t%s\t%s\n' \
            "~{contig}" "~{n_non_trv_input}" "${n_output}" "~{n_truvari_input}" "${n_truvari_output}" "${n_truvari_collapsed}" \
            >> ~{prefix}.merge_summary.tsv
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
