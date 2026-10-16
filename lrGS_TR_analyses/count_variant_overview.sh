#!/usr/bin/env bash
# Variant overview of a gnomAD_LR vep_parsed annotated bed: counts by variant class, all sites and PASS.
# Class: snv; del / ins split into indel (|allele_length| < 50) and SV (>= 50); trv; every other allele_type as is.
# Columns used: 8 FILTER, 9 allele_type, 10 allele_length.
# Usage: count_variant_overview.sh <annotated.bed.gz> <label> <out.tsv>
set -euo pipefail

if [ $# -ne 3 ]; then
    echo "Usage: $0 <annotated.bed.gz> <label> <out.tsv>" >&2
    exit 1
fi
bed=$1
label=$2
out=$3

gzip -cd "$bed" \
    | awk -F'\t' -v label="$label" 'BEGIN{OFS="\t"} NR>1 {
        t = $9; L = $10 + 0; if (L < 0) L = -L
        if (t == "del" || t == "ins") c = t "_" (L >= 50 ? "SV" : "indel")
        else c = t
        n[c]++; if ($8 == "PASS") p[c]++
      }
      END { for (c in n) print label, c, n[c], p[c] + 0 }' \
    | sort -t$'\t' -k3,3nr \
    | { printf "cohort\tvariant_class\tn_all\tn_PASS\n"; cat; } > "$out"
cat "$out"
