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
            "Every output record carries `MERGE_COUNT`, the number of input records merged into it, and `MERGE_TYPE`, one of EXACT, TRV_EXACT, TRUVARI or UNIQUE.",
            "Provenance is recorded in four parallel lists with one entry per merged input record: `SOURCE_NAMES`, the `vcf_names` entry of the callset that carried it, `SOURCE_IDS`, its ID there, and `SOURCE_REFS` and `SOURCE_ALTS`, its REF and ALT as that callset wrote them. The ALT alleles of a single record are separated by a pipe, so that a multiallelic record stays one entry. A name repeats when Truvari collapses records that came from the same callset.",
            "Where a merged record cannot hold both inputs' values, the ID and any INFO field other than `MERGE_COUNT` and the `SOURCE_` lists are taken from the first input VCF that carried the record, and AC, AN and AF are recomputed over the merged samples."
        ]
    }

    parameter_meta {
        vcfs: "Per-callset VCFs for the contig being merged, each called across a distinct set of samples."
        vcf_idxs: "Indexes for `vcfs`."
        vcf_names: "Name of each entry of `vcfs`, in the same order, recorded in `SOURCE_NAMES`. A name must not contain a comma, semicolon, pipe or whitespace."
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
        merge_summary_tsv: "TSV with one row per allele type and size bin, tandem repeats sharing a single row without a bin, and three columns per entry of `vcf_names` holding how many of that callset's records went in, how many came out merged with another record, and how many came out on their own. The merged and unmerged counts sum to the input count."
    }

    input {
        Array[File] vcfs
        Array[File] vcf_idxs
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

        Int shard_bin_size

        String utils_docker

        RuntimeAttr? runtime_attr_create_shards
        RuntimeAttr? runtime_attr_subset_shard
        RuntimeAttr? runtime_attr_merge_trv
        RuntimeAttr? runtime_attr_merge_non_trv
        RuntimeAttr? runtime_attr_consolidate_non_trv
        RuntimeAttr? runtime_attr_finalize_non_trv
        RuntimeAttr? runtime_attr_concat_merged
        RuntimeAttr? runtime_attr_concat_shards
        RuntimeAttr? runtime_attr_summarize_merge
    }

    call Helpers.CreateContigShards {
        input:
            vcfs = vcfs,
            vcf_idxs = vcf_idxs,
            contig = contig,
            shard_bin_size = shard_bin_size,
            prefix = "~{prefix}.shards",
            docker = utils_docker,
            runtime_attr_override = runtime_attr_create_shards
    }

    scatter (j in range(length(CreateContigShards.shard_regions))) {
        scatter (i in range(length(vcfs))) {
            call Helpers.SubsetVcfToRegion as SubsetCallsetToShard {
                input:
                    vcf = vcfs[i],
                    vcf_idx = vcf_idxs[i],
                    region = CreateContigShards.shard_regions[j],
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
                join_info_fields = ["SOURCE_IDS", "SOURCE_NAMES", "SOURCE_REFS", "SOURCE_ALTS"],
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

    call SummarizeMerge {
        input:
            vcfs = vcfs,
            vcf_idxs = vcf_idxs,
            vcf_names = vcf_names,
            merged_vcf = ConcatShards.concat_vcf,
            merged_vcf_idx = ConcatShards.concat_vcf_idx,
            prefix = "~{prefix}.merge_summary",
            docker = utils_docker,
            runtime_attr_override = runtime_attr_summarize_merge
    }

    output {
        File merged_vcf = ConcatShards.concat_vcf
        File merged_vcf_idx = ConcatShards.concat_vcf_idx
        File merge_summary_tsv = SummarizeMerge.summary_tsv
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

        echo '##INFO=<ID=SOURCE_IDS,Number=.,Type=String,Description="ID of each input record merged into this site, parallel to SOURCE_NAMES">' > merge_hdr.txt
        echo '##INFO=<ID=SOURCE_NAMES,Number=.,Type=String,Description="Name of the input VCF that carried each merged record; a name repeats when Truvari collapsed records from one callset">' >> merge_hdr.txt
        echo '##INFO=<ID=SOURCE_REFS,Number=.,Type=String,Description="REF allele of each merged input record as written in its callset, parallel to SOURCE_NAMES">' >> merge_hdr.txt
        echo '##INFO=<ID=SOURCE_ALTS,Number=.,Type=String,Description="ALT alleles of each merged input record as written in its callset, parallel to SOURCE_NAMES, with the ALTs of one record separated by |">' >> merge_hdr.txt
        echo '##INFO=<ID=MERGE_COUNT,Number=1,Type=Integer,Description="Number of input records merged into this site">' >> merge_hdr.txt
        echo '##INFO=<ID=MERGE_TYPE,Number=1,Type=String,Description="Merge strategy: EXACT, TRV_EXACT, TRUVARI, or UNIQUE">' >> merge_hdr.txt

        paste ~{write_lines(vcfs)} ~{write_lines(vcf_idxs)} ~{write_lines(vcf_names)} > vcf_triples.tsv

        i=0
        while IFS=$'\t' read -r vcf vcf_idx vcf_name; do
            if [[ "$vcf_idx" != "${vcf}.tbi" ]]; then
                ln -sf "$vcf_idx" "${vcf}.tbi"
            fi

            # Callsets type AL differently, and bcftools merge then writes raw bytes into the sample columns
            bcftools view -h "$vcf" \
                | sed 's|^##FORMAT=<ID=AL,.*$|##FORMAT=<ID=AL,Number=.,Type=Integer,Description="Length of each allele">|' \
                > "al_header_${i}.txt"

            bcftools reheader -h "al_header_${i}.txt" -o "retyped_${i}.vcf.gz" "$vcf"
            tabix -f -p vcf "retyped_${i}.vcf.gz"

            bcftools view -i 'INFO/allele_type=="trv"' -Oz -o "typed_${i}.vcf.gz" "retyped_${i}.vcf.gz"
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

            # Stamp provenance, separating one record's ALTs with '|' so the merged lists stay parseable
            bcftools query \
                -f "%CHROM\t%POS\t%REF\t%ALT\t%ID\t${vcf_name}\t%REF\t%ALT\t1\n" \
                "cleaned_${i}.vcf.gz" \
                | awk -F'\t' -v OFS='\t' '{ gsub(",", "|", $8); print }' \
                | bgzip > "annot_${i}.tsv.gz"

            tabix -s1 -b2 -e2 "annot_${i}.tsv.gz"

            bcftools annotate \
                -a "annot_${i}.tsv.gz" -h merge_hdr.txt \
                -c CHROM,POS,REF,ALT,SOURCE_IDS,SOURCE_NAMES,SOURCE_REFS,SOURCE_ALTS,MERGE_COUNT \
                -Oz -o "stamped_${i}.vcf.gz" \
                "cleaned_${i}.vcf.gz"

            tabix -f -p vcf "stamped_${i}.vcf.gz"

            # Set the ID to the REF so that 'bcftools merge -m id' keys on CHROM, POS and REF alone
            bcftools annotate --set-id '%REF' -Oz -o "tagged_${i}.vcf.gz" "stamped_${i}.vcf.gz"

            tabix -f -p vcf "tagged_${i}.vcf.gz"

            rm -f \
                "al_header_${i}.txt" \
                "retyped_${i}.vcf.gz" "retyped_${i}.vcf.gz.tbi" \
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
            bcftools merge --print-header -m id -l tagged_vcfs.list | bgzip > header_only.vcf.gz
            bcftools +fill-tags header_only.vcf.gz -Oz -o ~{prefix}.vcf.gz -- -t AC,AN,AF
            tabix -f -p vcf ~{prefix}.vcf.gz
            rm -f header_only.vcf.gz
            exit 0
        fi

        bcftools merge \
            -m id \
            -i MERGE_COUNT:sum,SOURCE_IDS:join,SOURCE_NAMES:join,SOURCE_REFS:join,SOURCE_ALTS:join \
            -Ov -o merged.raw.vcf \
            -l tagged_vcfs.list

        while read -r tagged_vcf; do
            rm -f "$tagged_vcf" "$tagged_vcf.tbi"
        done < tagged_vcfs.list

        # Undo the REF/ALT trimming and the ID set to REF; every record merged here shares one REF
        awk -F'\t' -v OFS='\t' '
            /^#/ { print; next }
            {
                orig = ""
                if (match($8, /(^|;)SOURCE_REFS=[^;]*/)) {
                    orig = substr($8, RSTART, RLENGTH)
                    sub(/^;/, "", orig)
                    sub(/^SOURCE_REFS=/, "", orig)
                    split(orig, refs, ",")
                    orig = refs[1]
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
                if (match($8, /(^|;)SOURCE_IDS=[^;]*/)) {
                    ids = substr($8, RSTART, RLENGTH)
                    sub(/^;/, "", ids)
                    sub(/^SOURCE_IDS=/, "", ids)
                    split(ids, id_list, ",")
                    $3 = id_list[1]
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

        bcftools annotate \
            -a merge_type.tsv.gz \
            -c CHROM,POS,REF,ALT,MERGE_TYPE \
            -Oz -o type_filled.vcf.gz \
            untrimmed.vcf.gz

        # Recompute the allele counts over the merged sample set, as the non-TR branch also does
        bcftools +fill-tags type_filled.vcf.gz -Oz -o ~{prefix}.vcf.gz -- -t AC,AN,AF

        tabix -f -p vcf ~{prefix}.vcf.gz

        rm -f untrimmed.vcf.gz untrimmed.vcf.gz.tbi type_filled.vcf.gz merge_type.tsv.gz*
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

        echo '##INFO=<ID=SOURCE_IDS,Number=.,Type=String,Description="ID of each input record merged into this site, parallel to SOURCE_NAMES">' > merge_hdr.txt
        echo '##INFO=<ID=SOURCE_NAMES,Number=.,Type=String,Description="Name of the input VCF that carried each merged record; a name repeats when Truvari collapsed records from one callset">' >> merge_hdr.txt
        echo '##INFO=<ID=SOURCE_REFS,Number=.,Type=String,Description="REF allele of each merged input record as written in its callset, parallel to SOURCE_NAMES">' >> merge_hdr.txt
        echo '##INFO=<ID=SOURCE_ALTS,Number=.,Type=String,Description="ALT alleles of each merged input record as written in its callset, parallel to SOURCE_NAMES, with the ALTs of one record separated by |">' >> merge_hdr.txt
        echo '##INFO=<ID=MERGE_COUNT,Number=1,Type=Integer,Description="Number of input records merged into this site">' >> merge_hdr.txt
        echo '##INFO=<ID=MERGE_TYPE,Number=1,Type=String,Description="Merge strategy: EXACT, TRV_EXACT, TRUVARI, or UNIQUE">' >> merge_hdr.txt

        paste ~{write_lines(vcfs)} ~{write_lines(vcf_idxs)} ~{write_lines(vcf_names)} > vcf_triples.tsv

        i=0
        while IFS=$'\t' read -r vcf vcf_idx vcf_name; do
            if [[ "$vcf_idx" != "${vcf}.tbi" ]]; then
                ln -sf "$vcf_idx" "${vcf}.tbi"
            fi

            # Callsets type AL differently, and bcftools merge then writes raw bytes into the sample columns
            bcftools view -h "$vcf" \
                | sed 's|^##FORMAT=<ID=AL,.*$|##FORMAT=<ID=AL,Number=.,Type=Integer,Description="Length of each allele">|' \
                > "al_header_${i}.txt"

            bcftools reheader -h "al_header_${i}.txt" -o "retyped_${i}.vcf.gz" "$vcf"
            tabix -f -p vcf "retyped_${i}.vcf.gz"

            bcftools view -e 'INFO/allele_type=="trv"' -Oz -o "typed_${i}.vcf.gz" "retyped_${i}.vcf.gz"
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

            # Stamp provenance, separating one record's ALTs with '|' so the merged lists stay parseable
            bcftools query \
                -f "%CHROM\t%POS\t%REF\t%ALT\t%ID\t${vcf_name}\t%REF\t%ALT\t1\n" \
                "cleaned_${i}.vcf.gz" \
                | awk -F'\t' -v OFS='\t' '{ gsub(",", "|", $8); print }' \
                | bgzip > "annot_${i}.tsv.gz"

            tabix -s1 -b2 -e2 "annot_${i}.tsv.gz"

            bcftools annotate \
                -a "annot_${i}.tsv.gz" -h merge_hdr.txt \
                -c CHROM,POS,REF,ALT,SOURCE_IDS,SOURCE_NAMES,SOURCE_REFS,SOURCE_ALTS,MERGE_COUNT \
                -Oz -o "stamped_${i}.vcf.gz" \
                "cleaned_${i}.vcf.gz"

            tabix -f -p vcf "stamped_${i}.vcf.gz"

            # Qualify the ID with the callset name, since Truvari needs the IDs it collapses to be unique
            bcftools annotate --set-id "${vcf_name}:%ID" -Oz -o "tagged_${i}.vcf.gz" "stamped_${i}.vcf.gz"

            tabix -f -p vcf "tagged_${i}.vcf.gz"

            rm -f \
                "al_header_${i}.txt" \
                "retyped_${i}.vcf.gz" "retyped_${i}.vcf.gz.tbi" \
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
            bcftools merge --print-header -m none -l tagged_vcfs.list | bgzip > empty.vcf.gz

            for split in matched unmatched_large unmatched_small; do
                cp empty.vcf.gz "~{prefix}.${split}.vcf.gz"
                tabix -f -p vcf "~{prefix}.${split}.vcf.gz"
            done

            exit 0
        fi

        bcftools merge \
            -m none \
            -i MERGE_COUNT:sum,SOURCE_IDS:join,SOURCE_NAMES:join,SOURCE_REFS:join,SOURCE_ALTS:join \
            -Ov -o exact_merged.raw.vcf \
            -l tagged_vcfs.list

        while read -r tagged_vcf; do
            rm -f "$tagged_vcf" "$tagged_vcf.tbi"
        done < tagged_vcfs.list

        # Undo the REF/ALT trimming; the ID keeps its callset prefix, which Truvari needs to stay unique
        awk -F'\t' -v OFS='\t' '
            /^#/ { print; next }
            {
                orig = ""
                if (match($8, /(^|;)SOURCE_REFS=[^;]*/)) {
                    orig = substr($8, RSTART, RLENGTH)
                    sub(/^;/, "", orig)
                    sub(/^SOURCE_REFS=/, "", orig)
                    split(orig, refs, ",")
                    orig = refs[1]
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
        ' exact_merged.raw.vcf | bgzip > exact_merged.vcf.gz

        tabix -f -p vcf exact_merged.vcf.gz

        rm -f exact_merged.raw.vcf

        # Split into matched (exact), unmatched large (truvari candidates), and unmatched small (unique)
        bcftools view \
            -i 'MERGE_COUNT>1' \
            -Oz -o matched.vcf.gz \
            exact_merged.vcf.gz

        bcftools view \
            -i 'MERGE_COUNT=1 && abs(INFO/allele_length) >= ~{min_truvari_match}' \
            -Oz -o ~{prefix}.unmatched_large.vcf.gz \
            exact_merged.vcf.gz

        # The complement of the split above, so that a record with no allele_length is kept rather than dropped
        bcftools view \
            -e 'MERGE_COUNT!=1 || abs(INFO/allele_length) >= ~{min_truvari_match}' \
            -Oz -o unmatched_small.vcf.gz \
            exact_merged.vcf.gz

        tabix -f -p vcf matched.vcf.gz
        tabix -f -p vcf ~{prefix}.unmatched_large.vcf.gz
        tabix -f -p vcf unmatched_small.vcf.gz

        rm -f exact_merged.vcf.gz exact_merged.vcf.gz.tbi

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

        bcftools annotate \
            -a missing_mt.tsv.gz \
            -c CHROM,POS,REF,ALT,MERGE_TYPE \
            -Oz -o type_filled.vcf.gz \
            sorted.vcf.gz

        tabix -f -p vcf type_filled.vcf.gz

        # Take back the ID that the callset prefix overwrote, now that Truvari no longer needs it
        bcftools view type_filled.vcf.gz \
            | awk -F'\t' -v OFS='\t' '
                /^#/ { print; next }
                {
                    if (match($8, /(^|;)SOURCE_IDS=[^;]*/)) {
                        ids = substr($8, RSTART, RLENGTH)
                        sub(/^;/, "", ids)
                        sub(/^SOURCE_IDS=/, "", ids)
                        split(ids, id_list, ",")
                        $3 = id_list[1]
                    }
                    print
                }
            ' | bgzip > id_restored.vcf.gz

        tabix -f -p vcf id_restored.vcf.gz

        bcftools +fill-tags id_restored.vcf.gz -Oz -o ~{prefix}.vcf.gz -- -t AC,AN,AF

        rm -f sorted.vcf.gz sorted.vcf.gz.tbi missing_mt.tsv.gz* \
            type_filled.vcf.gz* id_restored.vcf.gz*

        tabix -f -p vcf ~{prefix}.vcf.gz
    >>>

    output {
        File merged_vcf = "~{prefix}.vcf.gz"
        File merged_vcf_idx = "~{prefix}.vcf.gz.tbi"
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

task SummarizeMerge {
    input {
        Array[File] vcfs
        Array[File] vcf_idxs
        Array[String] vcf_names
        File merged_vcf
        File merged_vcf_idx
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        bcftools query \
            -f '%INFO/MERGE_COUNT\t%INFO/SOURCE_NAMES\t%INFO/SOURCE_IDS\n' \
            ~{merged_vcf} \
            > output_slots.tsv

        paste ~{write_lines(vcfs)} ~{write_lines(vcf_idxs)} ~{write_lines(vcf_names)} > vcf_triples.tsv

        : > input_records.tsv
        while IFS=$'\t' read -r vcf vcf_idx vcf_name; do
            if [[ "$vcf_idx" != "${vcf}.tbi" ]]; then
                ln -sf "$vcf_idx" "${vcf}.tbi"
            fi

            bcftools query \
                -f "${vcf_name}\t%ID\t%INFO/allele_type\t%INFO/allele_length\n" \
                "$vcf" \
                >> input_records.tsv
        done < vcf_triples.tsv

        python3 <<CODE
import sys

names = "~{sep='\t' vcf_names}".split("\t")

BINS = [
    ("0", lambda v: v == 0),
    ("1-19", lambda v: 1 <= v <= 19),
    ("20-49", lambda v: 20 <= v <= 49),
    ("50-99", lambda v: 50 <= v <= 99),
    ("100-499", lambda v: 100 <= v <= 499),
    ("500-4999", lambda v: 500 <= v <= 4999),
    ("5000+", lambda v: v >= 5000),
]
BIN_ORDER = {label: index for index, (label, _) in enumerate(BINS)}
BIN_ORDER["NA"] = len(BINS)


def size_bin(length):
    if length in (".", ""):
        return "NA"
    try:
        value = abs(int(float(length)))
    except ValueError:
        return "NA"
    for label, matches in BINS:
        if matches(value):
            return label
    return "NA"


# Every input record reaches the output inside exactly one slot of one merged record
merged_pairs = set()
unmerged_pairs = set()
total_slots = 0
with open("output_slots.tsv") as handle:
    for line in handle:
        count, source_names, source_ids = line.rstrip("\n").split("\t")
        count = int(count)
        record_names = source_names.split(",")
        record_ids = source_ids.split(",")
        if len(record_names) != count or len(record_ids) != count:
            sys.exit(
                "SOURCE_NAMES/SOURCE_IDS length does not match MERGE_COUNT: " + line
            )
        target = merged_pairs if count > 1 else unmerged_pairs
        for pair in zip(record_names, record_ids):
            if pair in merged_pairs or pair in unmerged_pairs:
                sys.exit("input record appears in more than one output record: %s" % (pair,))
            target.add(pair)
        total_slots += count

rows = {}
n_inputs = 0
with open("input_records.tsv") as handle:
    for line in handle:
        vcf_name, variant_id, allele_type, allele_length = line.rstrip("\n").split("\t")
        n_inputs += 1
        if allele_type == "trv":
            key = ("trv", "NA")
        else:
            key = (allele_type, size_bin(allele_length))
        row = rows.setdefault(key, {name: [0, 0, 0] for name in names})
        counts = row[vcf_name]
        counts[0] += 1
        pair = (vcf_name, variant_id)
        if pair in merged_pairs:
            counts[1] += 1
        elif pair in unmerged_pairs:
            counts[2] += 1
        else:
            sys.exit("input record is missing from the merged VCF: %s" % (pair,))

if total_slots != n_inputs:
    sys.exit("sum of MERGE_COUNT (%d) does not match the input record count (%d)" % (total_slots, n_inputs))

header = ["type", "size"]
for name in names:
    header += ["%s_input" % name, "%s_merged" % name, "%s_unmerged" % name]

with open("~{prefix}.tsv", "w") as out:
    out.write("\t".join(header) + "\n")
    for key in sorted(rows, key=lambda k: (k[0], BIN_ORDER.get(k[1], len(BINS)))):
        fields = list(key)
        for name in names:
            fields += [str(value) for value in rows[key][name]]
        out.write("\t".join(fields) + "\n")
CODE
    >>>

    output {
        File summary_tsv = "~{prefix}.tsv"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 16,
        disk_gb: 2 * ceil(size(vcfs, "GB") + size(merged_vcf, "GB")) + 25,
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
