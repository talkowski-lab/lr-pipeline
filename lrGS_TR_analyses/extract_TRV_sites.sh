#!/usr/bin/env bash
# Extract tandem-repeat variant (TRV) sites from a gnomAD_LR vep_parsed annotated bed.gz.
# A site is kept if its ID (col 4) contains "-TRV-" (e.g. chr1-10616-TRV-381). Header line is retained.
# Output is bgzipped + tabix-indexed, plus a FILTER count summary.
# Usage: extract_TRV_sites.sh <input.bed.gz> <output.bed.gz>
set -euo pipefail

if [ $# -ne 2 ]; then
    echo "Usage: $0 <input.bed.gz> <output.bed.gz>" >&2
    exit 1
fi
in_bed=$1
out_bed=$2

gzip -cd "$in_bed" \
    | awk -F'\t' 'NR==1 && /^#/ {print; next} $4 ~ /-TRV-/' \
    | bgzip -c > "$out_bed"
tabix -f -p bed -S 1 "$out_bed"

summary=${out_bed%.bed.gz}.FILTER_counts.tsv
gzip -cd "$out_bed" \
    | awk -F'\t' 'NR>1 {n[$8]++; tot++} END {print "FILTER\tn_sites"; for (f in n) print f"\t"n[f]; print "TOTAL\t"tot}' \
    > "$summary"
cat "$summary"
