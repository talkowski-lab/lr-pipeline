#!/bin/bash
# Generate the hprc_hgsvc long-read low-coverage blacklist at three
# sample-fraction thresholds: >=90%, >=70%, >=50% of eligible samples.
#
# Input:
#   hgsvc_hprc.cohort.coverage_counts.tsv
#     chrom, start, end, low_coverage_sample_count, eligible_sample_count
#     (100bp bins genome-wide; low_coverage flag per sample is already
#     computed using each sample's cutoff = 20% of that sample's overall
#     median coverage on autosomes/chrX, 10% of median on chrX/chrY for
#     male samples; eligible_sample_count excludes females from chrY)
#
# Output (for p in 90, 70, 50):
#   hgsvc_hprc.blacklist.p${p}.raw.bed(.gz)        -- merged bins where >=p%
#                                                      of eligible samples
#                                                      show low coverage
#   hgsvc_hprc.blacklist.p${p}.annotated.bed(.gz)  -- + size, + overlap bp
#                                                      with telomere/
#                                                      centromere/N-masked/
#                                                      SegDup/SimpRep,
#                                                      + primary_annotation
set -euo pipefail

DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
COUNTS="$DIR/hgsvc_hprc.cohort.coverage_counts.tsv"

# 1. threshold + merge for each of the 3 fractions
for p in 90 70 50; do
  frac=$(python3 -c "print($p/100)")
  awk -F'\t' -v OFS='\t' -v frac="$frac" 'NR>1 && ($4/$5) >= frac {print $1,$2,$3}' "$COUNTS" \
    | sort -k1,1 -k2,2n \
    | bedtools merge -i - \
    > "$DIR/hgsvc_hprc.blacklist.p${p}.raw.bed"
  echo "p${p}: $(wc -l < "$DIR/hgsvc_hprc.blacklist.p${p}.raw.bed") regions"
done

# 2. annotate all three with the 5 tracks + primary_annotation
python3 "$DIR/scripts/annotate_blacklist.py" 90 70 50

# 3. regenerate the downsampled coverage curve used for plotting
python3 "$DIR/scripts/downsample_coverage.py"

# 4. per-chromosome 3-panel plots (90/70/50, top to bottom)
python3 "$DIR/scripts/plot_blacklist_per_chr.py"

# 5. bgzip the bed outputs (bgzipped copies are what's tracked in git)
for p in 90 70 50; do
  bgzip -k -f "$DIR/hgsvc_hprc.blacklist.p${p}.raw.bed"
  bgzip -k -f "$DIR/hgsvc_hprc.blacklist.p${p}.annotated.bed"
done
