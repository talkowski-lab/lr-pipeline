#!/bin/bash
# One-pass cumulative-union counting for one chromosome VCF.
#
# For the rank-ordered sample list, the number of distinct FILTER=PASS sites where
# at least one of the top-i samples is non-ref (the "union" for i) equals the
# number of sites whose best-ranked non-ref carrier has rank <= i. So instead of
# one bcftools pass per i, this makes a single pass per VCF:
#   1. bcftools view -f PASS -S <ranked samples present in this VCF> -c 1
#      (subset to the ranked samples; keep PASS sites non-ref in that subset, the
#      same site set `-S <top-N> -c 1` keeps),
#   2. for each site, find the smallest rank among samples with a non-ref GT,
#   3. histogram sites by that rank, split by the same 7 type/size categories
#      (SNV, del/ins x 1-49 / 50-499 / >499 bp) plus the total.
# Summing these histograms over contigs and taking the cumulative sum over rank
# gives the cumulative union table (cumulative_union_from_min_rank.py).
#
# Usage: union_min_rank_histogram.sh <in.vcf.gz> <ordered_samples.txt> <out_prefix>
#   ordered_samples.txt: one sample ID per line, rank 1 first (RankSamples output)
# Writes:
#   <out_prefix>.min_rank_histogram.tsv: rank, n_total_nonref, n_snv, n_del_1_49,
#       n_ins_1_49, n_del_50_499, n_ins_50_499, n_del_gt499, n_ins_gt499
#       (sites whose best-ranked carrier has that rank; ranks with no sites omitted)
set -euo pipefail

VCF=$1
ORDERED_SAMPLES=$2
OUT_PREFIX=$3

WORKDIR=$(mktemp -d)
trap 'rm -rf "$WORKDIR"' EXIT

# Ranked samples present in this VCF, kept in rank order.
bcftools query -l "$VCF" > "$WORKDIR/vcf_samples.txt"
awk 'NR==FNR { in_vcf[$1] = 1; next } ($1 in in_vcf)' "$WORKDIR/vcf_samples.txt" "$ORDERED_SAMPLES" \
    > "$WORKDIR/present.txt"
n_ranked=$(wc -l < "$ORDERED_SAMPLES" | tr -d ' ')
n_present=$(wc -l < "$WORKDIR/present.txt" | tr -d ' ')
if [ "$n_present" -eq 0 ]; then
    echo "ERROR: none of the $n_ranked ranked samples are in $VCF" >&2
    exit 1
fi
if [ "$n_present" -lt "$n_ranked" ]; then
    echo "WARNING: $((n_ranked - n_present)) ranked samples absent from $VCF; they contribute no sites here" >&2
fi

bcftools view -f PASS -S "$WORKDIR/present.txt" -c 1 -Ou "$VCF" > "$WORKDIR/subset.bcf"
# Global rank of each output genotype column, in the subset file's own sample order.
bcftools query -l "$WORKDIR/subset.bcf" \
    | awk 'NR==FNR { rank[$1] = FNR; next } { print rank[$1] }' "$ORDERED_SAMPLES" - \
    > "$WORKDIR/column_ranks.txt"

bcftools query -f '%INFO/allele_type\t%INFO/allele_length[\t%GT]\n' "$WORKDIR/subset.bcf" \
    | awk -F'\t' -v rank_file="$WORKDIR/column_ranks.txt" '
        BEGIN {
            # Visit genotype columns in rank order, so the first non-ref one found
            # is the best-ranked carrier and the scan can stop there.
            n = 0
            while ((getline r < rank_file) > 0) { n++; col_of_rank[r] = n + 2 }
            m = 0
            for (r in col_of_rank) ranks[++m] = r + 0
            asort_n = m
            for (a = 2; a <= m; a++) {           # insertion sort (n <= a few thousand)
                v = ranks[a]; b = a - 1
                while (b >= 1 && ranks[b] > v) { ranks[b + 1] = ranks[b]; b-- }
                ranks[b + 1] = v
            }
            for (a = 1; a <= m; a++) order[a] = col_of_rank[ranks[a]]
        }
        {
            best = 0
            for (a = 1; a <= asort_n; a++) {
                if ($(order[a]) ~ /[1-9]/) { best = ranks[a]; break }
            }
            if (!best) { n_no_carrier++; next }   # cannot happen after -c 1; counted for safety
            total[best]++
            len = $2; if (len < 0) len = -len
            if ($1 == "snv")                                   snv[best]++
            else if ($1 == "del" && len >= 1  && len <= 49)    d1[best]++
            else if ($1 == "ins" && len >= 1  && len <= 49)    i1[best]++
            else if ($1 == "del" && len >= 50 && len <= 499)   d2[best]++
            else if ($1 == "ins" && len >= 50 && len <= 499)   i2[best]++
            else if ($1 == "del" && len > 499)                 d3[best]++
            else if ($1 == "ins" && len > 499)                 i3[best]++
        }
        END {
            print "rank\tn_total_nonref\tn_snv\tn_del_1_49\tn_ins_1_49\tn_del_50_499\tn_ins_50_499\tn_del_gt499\tn_ins_gt499"
            for (a = 1; a <= asort_n; a++) {
                r = ranks[a]
                if (total[r])
                    printf "%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\n", r, total[r], snv[r]+0, d1[r]+0, i1[r]+0,
                        d2[r]+0, i2[r]+0, d3[r]+0, i3[r]+0
            }
            if (n_no_carrier) printf "WARNING: %d sites with no non-ref carrier skipped\n", n_no_carrier > "/dev/stderr"
        }
    ' > "${OUT_PREFIX}.min_rank_histogram.tsv"

echo "Wrote ${OUT_PREFIX}.min_rank_histogram.tsv ($n_present of $n_ranked ranked samples in $VCF)"
