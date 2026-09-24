#!/bin/bash
# Per-sample non-ref genotype counts from one chromosome VCF, overall, by
# variant type/size (using the existing INFO/allele_type + INFO/allele_length
# annotations -- no need to re-derive type/size from REF/ALT), and by each
# type/size category further split by INFO/REGION (genomic context: US, RM,
# SD, SR -- see INFO/REGION header description; single-valued per site, so
# the split is mutually exclusive).
#
# For a filtered (or unfiltered) subset of sites, per-sample non-ref count is
# extracted from `bcftools stats -s -`'s PSC (SNV het/hom) + PSI (indel
# het/hom) sections combined:
#   count = PSC.nNonRefHom + PSC.nHets + PSI.nInsHets + PSI.nDelHets
#           + PSI.nInsAltHoms + PSI.nDelAltHoms
# (safe to always sum all six: after filtering to a single type, the
# irrelevant terms are 0 anyway; for the unfiltered "total" pass this
# correctly combines SNV + indel genotypes). Caveat inherited from bcftools
# itself: a het genotype with one ins allele and one del allele is counted
# in both nInsHets and nDelHets (documented bcftools behavior).
#
# Usage: per_sample_type_counts.sh <in.vcf.gz> <out_prefix>
# Writes:
#   <out_prefix>.n_samples.txt         -- single integer
#   <out_prefix>.per_sample_counts.tsv -- sample_id, n_total_nonref, then the
#                                          7 type/size categories (n_snv,
#                                          n_del_1_49, n_ins_1_49,
#                                          n_del_50_499, n_ins_50_499,
#                                          n_del_gt499, n_ins_gt499), then
#                                          each of those 7 further split by
#                                          REGION, e.g. n_snv_US, n_snv_RM,
#                                          n_snv_SD, n_snv_SR, ...
set -euo pipefail

VCF=$1
OUT_PREFIX=$2

REGIONS=(US RM SD SR)

bcftools query -l "$VCF" | wc -l | tr -d ' ' > "${OUT_PREFIX}.n_samples.txt"

extract_counts() {
    # $1 = bcftools -i filter expression ("" = no filter)
    local filt="$1"
    if [ -n "$filt" ]; then
        bcftools view -i "$filt" "$VCF"
    else
        bcftools view "$VCF"
    fi | bcftools stats -s - - | awk -F'\t' '
        $1=="PSC" { nonref[$3] += $5 + $6 }
        $1=="PSI" { nonref[$3] += $8 + $9 + $10 + $11 }
        END { for (s in nonref) print s"\t"nonref[s] }
    '
}

WORKDIR=$(mktemp -d)
trap 'rm -rf "$WORKDIR"' EXIT

# name<TAB>filter manifest, built up below, then run through extract_counts
# once per line (one bcftools view+stats pass per category/region combo).
COL_NAMES=()
COL_FILTERS=()

COL_NAMES+=("n_total_nonref");  COL_FILTERS+=("")
COL_NAMES+=("n_snv");           COL_FILTERS+=('INFO/allele_type="snv"')
COL_NAMES+=("n_del_1_49");      COL_FILTERS+=('INFO/allele_type="del" && abs(INFO/allele_length)>=1   && abs(INFO/allele_length)<=49')
COL_NAMES+=("n_ins_1_49");      COL_FILTERS+=('INFO/allele_type="ins" && abs(INFO/allele_length)>=1   && abs(INFO/allele_length)<=49')
COL_NAMES+=("n_del_50_499");    COL_FILTERS+=('INFO/allele_type="del" && abs(INFO/allele_length)>=50  && abs(INFO/allele_length)<=499')
COL_NAMES+=("n_ins_50_499");    COL_FILTERS+=('INFO/allele_type="ins" && abs(INFO/allele_length)>=50  && abs(INFO/allele_length)<=499')
COL_NAMES+=("n_del_gt499");     COL_FILTERS+=('INFO/allele_type="del" && abs(INFO/allele_length)>499')
COL_NAMES+=("n_ins_gt499");     COL_FILTERS+=('INFO/allele_type="ins" && abs(INFO/allele_length)>499')

# The 7 type/size categories above, each further split by REGION.
n_base_cats=${#COL_NAMES[@]}
for (( i=1; i<n_base_cats; i++ )); do
    base_name="${COL_NAMES[$i]}"
    base_filt="${COL_FILTERS[$i]}"
    for region in "${REGIONS[@]}"; do
        COL_NAMES+=("${base_name}_${region}")
        COL_FILTERS+=("${base_filt} && INFO/REGION=\"${region}\"")
    done
done

for (( i=0; i<${#COL_NAMES[@]}; i++ )); do
    extract_counts "${COL_FILTERS[$i]}" > "$WORKDIR/${i}.tsv"
done

# Join all per-category (sample -> count) files into one wide table. They
# all share the same sample set (bcftools stats -s - always emits every
# sample from the VCF header, count 0 if none matched the filter), so a
# single-pass awk merge keyed by sample is safe and avoids needing python3
# in this bcftools-only container. Row order doesn't matter here -- any
# downstream aggregation re-sorts/re-sums across chromosome shards.
{
    printf 'sample'
    for name in "${COL_NAMES[@]}"; do printf '\t%s' "$name"; done
    printf '\n'
} > "${OUT_PREFIX}.per_sample_counts.tsv"

awk -F'\t' -v workdir="$WORKDIR" -v n_cols="${#COL_NAMES[@]}" '
    BEGIN {
        for (k = 0; k < n_cols; k++) {
            f = workdir "/" k ".tsv"
            while ((getline line < f) > 0) {
                split(line, a, "\t")
                cnt[a[1], k] = a[2]
                samples[a[1]] = 1
            }
            close(f)
        }
        for (s in samples) {
            out = s
            for (k = 0; k < n_cols; k++) {
                key = s SUBSEP k
                out = out "\t" ((key in cnt) ? cnt[key] : 0)
            }
            print out
        }
    }
' >> "${OUT_PREFIX}.per_sample_counts.tsv"

echo "Wrote ${OUT_PREFIX}.per_sample_counts.tsv and ${OUT_PREFIX}.n_samples.txt"
