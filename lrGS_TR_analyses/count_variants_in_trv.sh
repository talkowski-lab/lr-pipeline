#!/usr/bin/env bash
# For every TRV, count the routine (non-TRV) variants that fall inside it, i.e. rows whose TRID column (col 13) holds
# that TRV's ID, and how many of them are lrGS-unique: "." in both dbGaP_ID (col 14) and gnomAD_V4_match_ID (col 16).
# Only FILTER PASS rows (col 8) with summed AC (col 19) > 0 are counted.
# Output TSV (no header order guarantee besides the header line): TRV_ID n_variants n_novel n_snv n_novel_snv
# Usage: count_variants_in_trv.sh <annotated.bed.gz> <out.tsv.gz>
set -euo pipefail

if [ $# -ne 2 ]; then
    echo "Usage: $0 <annotated.bed.gz> <out.tsv.gz>" >&2
    exit 1
fi
bed=$1
out=$2

gzip -cd "$bed" \
    | awk -F'\t' 'BEGIN{OFS="\t"} NR>1 && $9 != "trv" && $13 != "." && $8 == "PASS" {
        k = split($19, a, ","); s = 0; for (i = 1; i <= k; i++) s += a[i]
        if (s <= 0) next
        nv = ($14 == "." && $16 == ".")
        n[$13]++; nn[$13] += nv
        if ($9 == "snv") { ns[$13]++; nns[$13] += nv }
      }
      END {
        print "TRV_ID", "n_variants", "n_novel", "n_snv", "n_novel_snv"
        for (t in n) print t, n[t], nn[t], ns[t] + 0, nns[t] + 0
      }' \
    | gzip -c > "$out"
gzip -cd "$out" | awk -F'\t' 'NR>1{t++; v+=$2; nv+=$3; if($3==$2) all++} END{printf "TRVs with >=1 variant: %d; variants: %d; novel: %d (%.1f%%); TRVs with all variants novel: %d\n", t, v, nv, 100*nv/v, all}'
