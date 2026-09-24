#!/bin/bash
# For the top-i samples in a rank-ordered sample list, count the number of
# distinct FILTER=PASS variant SITES where at least one of those i samples
# carries a non-ref genotype (i.e. the union of each sample's non-ref sites,
# not the sum -- a site shared by 2 of the i samples is counted once), split
# by the same 7 type/size categories used elsewhere (SNV, del/ins x
# 1-49bp/50-499bp/>499bp). Implements this exactly the way bcftools does it
# natively: `bcftools view -S <sample_list> -c 1` re-computes INFO/AC from
# only the listed samples' genotypes and keeps sites where that recomputed
# AC is >=1, i.e. sites non-ref in the subset.
#
# Usage: cumulative_union_counts.sh <ordered_samples.txt> <i> <out_prefix> <vcf1> [<vcf2> ...]
# Writes:
#   <out_prefix>.tsv -- one row: sample_id(the i-th/newest-ranked sample),
#                        i, n_total_nonref, n_snv, n_del_1_49, n_ins_1_49,
#                        n_del_50_499, n_ins_50_499, n_del_gt499, n_ins_gt499
#                        (all counts = union over the top-i samples)
set -euo pipefail

ORDERED_SAMPLES=$1
I=$2
OUT_PREFIX=$3
shift 3
VCFS=("$@")

WORKDIR=$(mktemp -d)
trap 'rm -rf "$WORKDIR"' EXIT

head -n "$I" "$ORDERED_SAMPLES" > "$WORKDIR/sample_list.txt"
SAMPLE_ID=$(tail -n 1 "$WORKDIR/sample_list.txt")

: > "$WORKDIR/per_vcf_counts.tsv"
for vcf in "${VCFS[@]}"; do
    bcftools view -f PASS -S "$WORKDIR/sample_list.txt" -c 1 "$vcf" \
        | bcftools query -f '%INFO/allele_type\t%INFO/allele_length\n' \
        | awk -F'\t' '
            {
                total++
                len = $2; if (len < 0) len = -len
                if ($1 == "snv")                              snv++
                else if ($1 == "del" && len >= 1  && len <= 49)  del_1_49++
                else if ($1 == "ins" && len >= 1  && len <= 49)  ins_1_49++
                else if ($1 == "del" && len >= 50 && len <= 499) del_50_499++
                else if ($1 == "ins" && len >= 50 && len <= 499) ins_50_499++
                else if ($1 == "del" && len > 499)               del_gt499++
                else if ($1 == "ins" && len > 499)               ins_gt499++
            }
            END {
                printf "%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\n",
                    total+0, snv+0, del_1_49+0, ins_1_49+0,
                    del_50_499+0, ins_50_499+0, del_gt499+0, ins_gt499+0
            }
        ' >> "$WORKDIR/per_vcf_counts.tsv"
done

COUNTS=$(awk -F'\t' '
    { for (k = 1; k <= 8; k++) sum[k] += $k }
    END { for (k = 1; k <= 8; k++) printf "%s%s", sum[k]+0, (k < 8 ? "\t" : "\n") }
' "$WORKDIR/per_vcf_counts.tsv")

printf '%s\t%s\t%s\n' "$SAMPLE_ID" "$I" "$COUNTS" > "${OUT_PREFIX}.tsv"

echo "Wrote ${OUT_PREFIX}.tsv (i=${I}, sample=${SAMPLE_ID})"
