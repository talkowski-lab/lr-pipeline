#!/usr/bin/env bash
# Assign each TRV site to a mutually exclusive genic context using a GENCODE GTF:
#   coding (overlaps CDS) > UTR (overlaps exon, no CDS) > intronic (inside transcript span) > intergenic.
# Overlap = >=1 bp between the TR locus [START, END) and the feature.
# Also emits per-site size fields: locus length, min motif length across TRID components,
# and the per-ALT-allele length change (len(ALT) - len(REF)).
# Usage: annotate_TRV_genic_context.sh <TRV.bed.gz> <gencode.gtf.gz> <out.tsv.gz>
set -euo pipefail

if [ $# -ne 3 ]; then
    echo "Usage: $0 <TRV.bed.gz> <gencode.gtf.gz> <out.tsv.gz>" >&2
    exit 1
fi
trv_bed=$1
gtf=$2
out_tsv=$3

tmp=$(mktemp -d)
trap 'rm -rf "$tmp"' EXIT

# Feature beds (0-based half-open) labelled with gene_name, sorted and merged per gene name
for feat in CDS exon transcript; do
    gzip -cd "$gtf" \
        | awk -F'\t' -v f="$feat" 'BEGIN{OFS="\t"} $3==f {
            match($9, /gene_name "[^"]+"/); g=substr($9, RSTART+11, RLENGTH-12);
            print $1, $4-1, $5, g}' \
        | sort -k1,1 -k2,2n > "$tmp/$feat.bed"
done

# Slim site table: chrom start end ID FILTER SOURCE REGION n_TRID min_motif_len locus_len n_alt alt_len_diffs
gzip -cd "$trv_bed" \
    | awk -F'\t' 'BEGIN{OFS="\t"} NR>1 {
        k=split($13, t, ","); mm=-1;
        for (i=1; i<=k; i++) { m=split(t[i], p, "-"); ml=length(p[m]); if (mm<0 || ml<mm) mm=ml }
        n=split($6, a, ","); d="";
        for (i=1; i<=n; i++) d=d (i>1 ? "," : "") (length(a[i]) - length($5));
        print $1, $2, $3, $4, $8, $11, $12, k, mm, $3-$2, n, d}' \
    | sort -k1,1 -k2,2n > "$tmp/sites.bed"

bedtools map -a "$tmp/sites.bed" -b "$tmp/CDS.bed" -c 4 -o distinct \
    | bedtools map -a - -b "$tmp/exon.bed" -c 4 -o distinct \
    | bedtools map -a - -b "$tmp/transcript.bed" -c 4 -o distinct \
    | awk -F'\t' 'BEGIN{OFS="\t";
        print "chrom","start","end","ID","FILTER","SOURCE","REGION","n_TRID_components","min_motif_len",
              "locus_len","n_alt","alt_len_diffs","genic_context","genes"}
      {
        if ($13 != ".")      { c="coding";   g=$13 }
        else if ($14 != ".") { c="UTR";      g=$14 }
        else if ($15 != ".") { c="intronic"; g=$15 }
        else                 { c="intergenic"; g="." }
        print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,c,g
      }' \
    | gzip -c > "$out_tsv"

gzip -cd "$out_tsv" \
    | awk -F'\t' 'NR>1 {n[$13]++; if ($5=="PASS") p[$13]++; tot++}
        END {print "genic_context\tn_sites\tn_PASS"; for (c in n) print c"\t"n[c]"\t"p[c]+0; print "TOTAL\t"tot}'
