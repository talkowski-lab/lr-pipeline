#!/usr/bin/env bash
# Merge the annotated TRV site BED (AnalyzeTRVariants + annotate_trv_sites.py) with the columns that exist only in the
# VEP-parsed TRV bed (extract_TRV_sites.sh output): full REF / ALT sequences, SOURCE (TR catalog), REGION, nhomref / nhet /
# nhomalt (bi-allelic sites only) and the region flags in_telomere, in_centromere, in_nmasked, in_srGS_lowcov,
# in_lrGS_lowcov, in_centromere_blacklist. The two inputs must list the same TRVs in the same order (checked on ID).
# Usage: merge_trv_tables.sh <sites.annotated.bed.gz> <vep_parsed.annotated.TRV.bed.gz> <out.bed.gz>
set -euo pipefail

if [ $# -ne 3 ]; then
    echo "Usage: $0 <sites.annotated.bed.gz> <vep_parsed.annotated.TRV.bed.gz> <out.bed.gz>" >&2
    exit 1
fi
sites=$1
vep=$2
out=$3

if ! cmp -s <(gzip -cd "$sites" | cut -f4) <(gzip -cd "$vep" | cut -f4); then
    echo "Error: the two inputs do not list the same IDs in the same order" >&2
    exit 1
fi

# VEP-parsed columns kept: 5 REF, 6 ALT, 11 SOURCE, 12 REGION, 41-43 nhomref nhet nhomalt, 44-49 region flags
paste <(gzip -cd "$sites") <(gzip -cd "$vep" | cut -f5,6,11,12,41-49) | bgzip -c > "$out"
tabix -f -p bed -S 1 "$out"
echo "$out: $(gzip -cd "$out" | tail -n +2 | wc -l | tr -d ' ') rows, $(gzip -cd "$out" | sed -n '1p' | awk -F'\t' '{print NF}') columns"
