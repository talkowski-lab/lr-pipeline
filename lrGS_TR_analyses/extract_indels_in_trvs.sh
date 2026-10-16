#!/usr/bin/env bash
# Extract routine indels (allele_type del / ins) that fall inside a given set of TRVs (TRID column = TRV ID) from a
# gnomAD-LR annotated bed: FILTER PASS, summed AC > 0. Output TSV: TRV_ID variant_ID allele_type abs_length dbGaP_ID
# gnomAD_V4_match_ID srGS_captured (0 if "." in both match columns, else 1).
# Usage: extract_indels_in_trvs.sh <annotated.bed.gz> <trv_ids.txt> <out.tsv.gz>
set -euo pipefail

if [ $# -ne 3 ]; then
    echo "Usage: $0 <annotated.bed.gz> <trv_ids.txt> <out.tsv.gz>" >&2
    exit 1
fi
bed=$1
ids=$2
out=$3

gzip -cd "$bed" \
    | awk -F'\t' -v idf="$ids" 'BEGIN{OFS="\t"; while ((getline l < idf) > 0) keep[l] = 1
        print "TRV_ID", "variant_ID", "allele_type", "abs_length", "dbGaP_ID", "gnomAD_V4_match_ID", "srGS_captured"}
      NR>1 && ($9 == "del" || $9 == "ins") && ($13 in keep) && $8 == "PASS" {
        k = split($19, a, ","); s = 0; for (i = 1; i <= k; i++) s += a[i]
        if (s <= 0) next
        L = $10 + 0; if (L < 0) L = -L
        print $13, $4, $9, L, $14, $16, ($14 == "." && $16 == "." ? 0 : 1)
      }' \
    | gzip -c > "$out"
