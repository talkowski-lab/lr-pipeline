#!/usr/bin/env bash
# Regions low-coverage in srGS but well covered in lrGS, and their overlap with
# TR catalogs (STR = shortest motif <= MAX_STR_MOTIF bp, VNTR = longer).
#
# A 100bp bin qualifies if srGS low5_frac >= SR_MIN and lrGS low-coverage
# sample fraction < LR_MAX. Qualifying bins are merged into regions; each
# region is then labelled by the TR loci it overlaps.
#
# Usage:
#   srGS_low_lrGS_ok_TR_overlap.sh <srGS_bins.bed.gz> <lrGS_coverage_counts.tsv> \
#       <out_dir> <SR_MIN> <LR_MAX> <MAX_STR_MOTIF> <exclude_chroms(comma)> \
#       <segdup.bed.gz> <simprep.bed.gz> <centromere.bed.gz> LABEL=catalog.bed.gz [LABEL=catalog.bed.gz ...]
#
# Outputs (in out_dir):
#   regions.bed.gz         chrom start end size str_bp vntr_bp tr_bp segdup_bp simprep_bp cen_bp class majority
#                          class: STR | VNTR | STR+VNTR | none (any catalog overlap); *_bp are
#                          union bp of each track inside the region; majority: first track
#                          covering >= 50% of the region, in priority order centromere >
#                          catalog TR (STR or VNTR, whichever has more bp) > TRF simple repeat >
#                          SegDup > other
#   tr_loci.bed.gz         catalog loci used: chrom start end type catalog min_motif
#   summary.tsv            counts / bp by class, with any-overlap and >=50%-covered definitions
set -euo pipefail

sr_bins=$1
lr_counts=$2
out_dir=$3
sr_min=$4
lr_max=$5
max_str=$6
exclude=$7
segdup=$8
simprep=$9
centromere=${10}
shift 10

mkdir -p "$out_dir"
tmp=$(mktemp -d)
trap 'rm -rf "$tmp"' EXIT

# Qualifying bins: join the two bin grids line by line (identical bin order).
paste <(gzcat "$sr_bins") <(tail -n +2 "$lr_counts") \
  | awk -F'\t' -v OFS='\t' -v s="$sr_min" -v l="$lr_max" -v ex="$exclude" '
    BEGIN { n = split(ex, e, ","); for (i = 1; i <= n; i++) skip[e[i]] = 1 }
    $1 != $8 || $2 != $9 { print "bin grid mismatch at " $1 ":" $2 > "/dev/stderr"; exit 1 }
    ($1 in skip) || $12 == 0 { next }
    $6 >= s && $11 / $12 < l { print $1, $2, $3 }' \
  | bedtools merge > "$tmp/regions.bed"

# TR loci from each catalog, typed by shortest motif length.
for spec in "$@"; do
  label=${spec%%=*}
  gzcat "${spec#*=}" | awk -F'\t' -v OFS='\t' -v lab="$label" -v m="$max_str" '{
      motifs = $4; sub(/.*MOTIFS=/, "", motifs); sub(/;.*/, "", motifs)
      k = split(motifs, a, ","); mn = length(a[1])
      for (i = 2; i <= k; i++) if (length(a[i]) < mn) mn = length(a[i])
      print $1, $2, $3, (mn <= m ? "STR" : "VNTR"), lab, mn }'
done | sort -k1,1 -k2,2n > "$tmp/tr.bed"
gzip -c "$tmp/tr.bed" > "$out_dir/tr_loci.bed.gz"

awk '$4 == "STR"' "$tmp/tr.bed" | cut -f1-3 | bedtools merge > "$tmp/str.bed"
awk '$4 == "VNTR"' "$tmp/tr.bed" | cut -f1-3 | bedtools merge > "$tmp/vntr.bed"
cut -f1-3 "$tmp/tr.bed" | bedtools merge > "$tmp/tr_all.bed"
gzcat "$segdup" | cut -f1-3 | sort -k1,1 -k2,2n | bedtools merge > "$tmp/segdup.bed"
gzcat "$simprep" | cut -f1-3 | sort -k1,1 -k2,2n | bedtools merge > "$tmp/simprep.bed"
gzcat "$centromere" | cut -f1-3 | sort -k1,1 -k2,2n | bedtools merge > "$tmp/cen.bed"

cov_bp() { bedtools coverage -a "$tmp/regions.bed" -b "$1" -sorted -g "$tmp/genome" | cut -f5; }
cut -f1 "$tmp/regions.bed" "$tmp/tr.bed" "$tmp/segdup.bed" "$tmp/simprep.bed" "$tmp/cen.bed" \
  | awk '!seen[$1]++ { print $1 "\t" 1e10 }' > "$tmp/genome"
sort -k1,1 "$tmp/genome" -o "$tmp/genome"
sort -k1,1 -k2,2n "$tmp/regions.bed" -o "$tmp/regions.bed"
for t in str vntr tr_all segdup simprep cen; do
  sort -k1,1 -k2,2n "$tmp/$t.bed" -o "$tmp/$t.bed"
  cov_bp "$tmp/$t.bed" > "$tmp/$t.cov"
done

paste "$tmp/regions.bed" "$tmp/str.cov" "$tmp/vntr.cov" "$tmp/tr_all.cov" "$tmp/segdup.cov" "$tmp/simprep.cov" "$tmp/cen.cov" \
  | awk -F'\t' -v OFS='\t' '{
      size = $3 - $2; h = size / 2
      c = ($4 > 0 && $5 > 0) ? "STR+VNTR" : ($4 > 0 ? "STR" : ($5 > 0 ? "VNTR" : "none"))
      if ($9 >= h) m = "centromere"
      else if ($6 >= h) m = ($4 >= $5) ? "STR" : "VNTR"
      else if ($8 >= h) m = "TRF_simple_repeat(non-catalog)"
      else if ($7 >= h) m = "SegDup"
      else m = "other"
      print $1, $2, $3, size, $4, $5, $6, $7, $8, $9, c, m }' \
  | sort -k1,1V -k2,2n > "$tmp/annot.bed"
gzip -c "$tmp/annot.bed" > "$out_dir/regions.bed.gz"

# Summary: any-overlap class, and majority-bp class (see header).
awk -F'\t' -v OFS='\t' '
  { n++; bp += $4; cn[$11]++; cb[$11] += $4; hn[$12]++; hb[$12] += $4
    trbp += $7; sdbp += $8; srbp += $9; cenbp += $10 }
  END {
    print "definition", "class", "regions", "pct_regions", "bp", "pct_bp"
    print "all", "all", n, 100, bp, 100
    split("STR STR+VNTR VNTR none", k1, " ")
    for (i = 1; i <= 4; i++) print "any_overlap", k1[i], cn[k1[i]] + 0, 100 * cn[k1[i]] / n, cb[k1[i]] + 0, 100 * cb[k1[i]] / bp
    split("centromere STR VNTR TRF_simple_repeat(non-catalog) SegDup other", k2, " ")
    for (i = 1; i <= 6; i++) print "majority_bp", k2[i], hn[k2[i]] + 0, 100 * hn[k2[i]] / n, hb[k2[i]] + 0, 100 * hb[k2[i]] / bp
    print "bp_total", "TR_union_bp", "", "", trbp, 100 * trbp / bp
    print "bp_total", "SimpleRepeat_bp", "", "", srbp, 100 * srbp / bp
    print "bp_total", "SegDup_bp", "", "", sdbp, 100 * sdbp / bp
    print "bp_total", "centromere_bp", "", "", cenbp, 100 * cenbp / bp }' "$tmp/annot.bed" > "$out_dir/summary.tsv"

column -t -s $'\t' "$out_dir/summary.tsv"
