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
# The VCFs (one per chromosome) are processed in parallel, N_PARALLEL at a
# time (default: number of CPUs), and their counts summed.
#
# BY_ALLELE_TYPE (default true): split counts by type/size, which needs
# INFO/allele_type and INFO/allele_length (checked in every VCF before counting;
# fails with a clear message if missing). false: only the total union count.
# ALLELE_TYPE_TAG (default allele_type) / ALLELE_LENGTH_TAG (default allele_length):
# names of the INFO tags holding the variant type and length (e.g. SVTYPE / SVLEN).
# Type values are matched as snv/del/ins in either case (snv or SNV, del or DEL, ...).
#
# Usage: [N_PARALLEL=n] [BY_ALLELE_TYPE=true|false] [ALLELE_TYPE_TAG=allele_type] [ALLELE_LENGTH_TAG=allele_length] cumulative_union_counts.sh <ordered_samples.txt> <i> <out_prefix> <vcf1> [<vcf2> ...]
# Writes:
#   <out_prefix>.tsv -- one row: sample_id(the i-th/newest-ranked sample),
#                        i, n_total_nonref, then (BY_ALLELE_TYPE=true only)
#                        n_snv, n_del_1_49, n_ins_1_49, n_del_50_499,
#                        n_ins_50_499, n_del_gt499, n_ins_gt499
#                        (all counts = union over the top-i samples)
set -euo pipefail

ORDERED_SAMPLES=$1
I=$2
OUT_PREFIX=$3
shift 3
VCFS=("$@")
N_PARALLEL=${N_PARALLEL:-$(nproc 2>/dev/null || echo 1)}
BY_ALLELE_TYPE=${BY_ALLELE_TYPE:-true}
ALLELE_TYPE_TAG=${ALLELE_TYPE_TAG:-allele_type}
ALLELE_LENGTH_TAG=${ALLELE_LENGTH_TAG:-allele_length}
if [ "$BY_ALLELE_TYPE" != "true" ] && [ "$BY_ALLELE_TYPE" != "false" ]; then
    echo "ERROR: BY_ALLELE_TYPE must be true or false, got '$BY_ALLELE_TYPE'" >&2
    exit 1
fi
if [ "$BY_ALLELE_TYPE" = "true" ]; then
    for vcf in "${VCFS[@]}"; do
        for tag in "$ALLELE_TYPE_TAG" "$ALLELE_LENGTH_TAG"; do
            if ! bcftools view -h "$vcf" | grep -q "^##INFO=<ID=${tag},"; then
                echo "ERROR: INFO/$tag is not defined in $vcf, so sites cannot be split by type/size." >&2
                echo "       Set split_by_allele_type = false (BY_ALLELE_TYPE=false) to count total sites only." >&2
                exit 1
            fi
        done
    done
fi
N_COLS=8
[ "$BY_ALLELE_TYPE" = "true" ] || N_COLS=1

WORKDIR=$(mktemp -d)
trap 'rm -rf "$WORKDIR"' EXIT

head -n "$I" "$ORDERED_SAMPLES" > "$WORKDIR/sample_list.txt"
SAMPLE_ID=$(tail -n 1 "$WORKDIR/sample_list.txt")
if [ "$(wc -l < "$WORKDIR/sample_list.txt" | tr -d ' ')" -ne "$I" ]; then
    echo "ERROR: i=$I but $ORDERED_SAMPLES has fewer samples" >&2
    exit 1
fi

# Counts for one VCF -> $WORKDIR/counts.<k>.tsv (one line: 8 columns, or 1 if
# BY_ALLELE_TYPE=false).
count_one() {
    set -euo pipefail
    local k=$1 vcf=$2
    if [ "$BY_ALLELE_TYPE" = "false" ]; then
        bcftools view -H -f PASS -S "$WORKDIR/sample_list.txt" -c 1 "$vcf" | wc -l | tr -d ' ' \
            > "$WORKDIR/counts.$k.tsv"
        return
    fi
    bcftools view -f PASS -S "$WORKDIR/sample_list.txt" -c 1 "$vcf" \
        | bcftools query -f "%INFO/${ALLELE_TYPE_TAG}\t%INFO/${ALLELE_LENGTH_TAG}\n" \
        | awk -F'\t' '
            {
                total++
                $1 = tolower($1)
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
        ' > "$WORKDIR/counts.$k.tsv"
}
export -f count_one
export WORKDIR BY_ALLELE_TYPE ALLELE_TYPE_TAG ALLELE_LENGTH_TAG

# Pairs "k vcf" -> xargs runs count_one on N_PARALLEL of them at a time; a failure
# in any VCF makes xargs (and so this script) exit non-zero.
for k in "${!VCFS[@]}"; do printf '%s\n%s\n' "$k" "${VCFS[$k]}"; done \
    | xargs -n 2 -P "$N_PARALLEL" bash -c 'count_one "$1" "$2"' _

n_done=$(ls "$WORKDIR"/counts.*.tsv | wc -l | tr -d ' ')
if [ "$n_done" -ne "${#VCFS[@]}" ]; then
    echo "ERROR: counted $n_done of ${#VCFS[@]} VCFs" >&2
    exit 1
fi

COUNTS=$(cat "$WORKDIR"/counts.*.tsv | awk -F'\t' -v n="$N_COLS" '
    { for (k = 1; k <= n; k++) sum[k] += $k }
    END { for (k = 1; k <= n; k++) printf "%s%s", sum[k]+0, (k < n ? "\t" : "\n") }
')

printf '%s\t%s\t%s\n' "$SAMPLE_ID" "$I" "$COUNTS" > "${OUT_PREFIX}.tsv"

echo "Wrote ${OUT_PREFIX}.tsv (i=${I}, sample=${SAMPLE_ID}, ${#VCFS[@]} VCFs, ${N_PARALLEL} in parallel)"
