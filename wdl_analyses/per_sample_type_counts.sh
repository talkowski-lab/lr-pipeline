#!/bin/bash
# Per-sample non-ref genotype counts from one chromosome VCF, overall and by
# variant type/size (using the existing INFO/allele_type + INFO/allele_length
# annotations -- no need to re-derive type/size from REF/ALT).
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
#   <out_prefix>.n_samples.txt        -- single integer
#   <out_prefix>.per_sample_counts.tsv -- sample_id, n_total_nonref, n_snv,
#                                          n_del_1_49, n_ins_1_49,
#                                          n_del_50_499, n_ins_50_499,
#                                          n_del_gt499, n_ins_gt499
set -euo pipefail

VCF=$1
OUT_PREFIX=$2

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

extract_counts ''                                                                              > "$WORKDIR/1.total.tsv"
extract_counts 'INFO/allele_type="snv"'                                                        > "$WORKDIR/2.snv.tsv"
extract_counts 'INFO/allele_type="del" && abs(INFO/allele_length)>=1   && abs(INFO/allele_length)<=49'  > "$WORKDIR/3.del_1_49.tsv"
extract_counts 'INFO/allele_type="ins" && abs(INFO/allele_length)>=1   && abs(INFO/allele_length)<=49'  > "$WORKDIR/4.ins_1_49.tsv"
extract_counts 'INFO/allele_type="del" && abs(INFO/allele_length)>=50  && abs(INFO/allele_length)<=499' > "$WORKDIR/5.del_50_499.tsv"
extract_counts 'INFO/allele_type="ins" && abs(INFO/allele_length)>=50  && abs(INFO/allele_length)<=499' > "$WORKDIR/6.ins_50_499.tsv"
extract_counts 'INFO/allele_type="del" && abs(INFO/allele_length)>499'                          > "$WORKDIR/7.del_gt499.tsv"
extract_counts 'INFO/allele_type="ins" && abs(INFO/allele_length)>499'                          > "$WORKDIR/8.ins_gt499.tsv"

# Join the 8 per-category (sample -> count) files into one wide table. All 8
# share the same sample set (bcftools stats -s - always emits every sample
# from the VCF header, count 0 if none matched the filter), so a single-pass
# awk merge keyed by sample is safe and avoids needing python3 in this
# bcftools-only container. Row order doesn't matter here -- the downstream
# combiner re-sorts by total non-ref count across all chromosome shards.
awk -F'\t' -v f1="$WORKDIR/1.total.tsv" -v f2="$WORKDIR/2.snv.tsv" -v f3="$WORKDIR/3.del_1_49.tsv" \
    -v f4="$WORKDIR/4.ins_1_49.tsv" -v f5="$WORKDIR/5.del_50_499.tsv" -v f6="$WORKDIR/6.ins_50_499.tsv" \
    -v f7="$WORKDIR/7.del_gt499.tsv" -v f8="$WORKDIR/8.ins_gt499.tsv" '
    BEGIN {
        files[1]=f1; files[2]=f2; files[3]=f3; files[4]=f4
        files[5]=f5; files[6]=f6; files[7]=f7; files[8]=f8
        for (k = 1; k <= 8; k++) {
            while ((getline line < files[k]) > 0) {
                split(line, a, "\t")
                cnt[a[1], k] = a[2]
                samples[a[1]] = 1
            }
            close(files[k])
        }
        print "sample\tn_total_nonref\tn_snv\tn_del_1_49\tn_ins_1_49\tn_del_50_499\tn_ins_50_499\tn_del_gt499\tn_ins_gt499"
        for (s in samples) {
            out = s
            for (k = 1; k <= 8; k++) {
                key = s SUBSEP k
                out = out "\t" ((key in cnt) ? cnt[key] : 0)
            }
            print out
        }
    }
' > "${OUT_PREFIX}.per_sample_counts.tsv"

echo "Wrote ${OUT_PREFIX}.per_sample_counts.tsv and ${OUT_PREFIX}.n_samples.txt"
