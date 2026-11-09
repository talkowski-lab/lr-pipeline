#!/usr/bin/env bash
# Extract the depth (RD), BAF and genotype (GT) tables for every breakpoint SD of a GD cluster, plus the depth across the
# whole span, with samples renamed (prepare_region_data.sh per window). Example: 22q11.21 breakpoints A-D.
#
# Windows (from the GD NAHR breakpoint-SD table, breakpoint_sd_windows.py):
#   each SD +/- --flank-frac x its length (default half)      -> <prefix>.<cluster>_<BP>.renamed.{RD.txt,baf.tsv,GT.tsv}.gz
#   whole span: first SD start - flank of the first SD to last SD end + flank of the last SD
#                                                              -> <prefix>.<cluster>_<BP1>-<BPn>.region.renamed.RD.txt.gz
#   and <prefix>.<cluster>.windows.bed listing all windows.
#
# Usage:
#   extract_gd_cluster_data.sh --cluster 22q11 --breakpoints A,B,C,D --bp-sd GD_NAHR.breakpoint_SDs.bed.gz \
#     --rd aou_lrGS_GD_carriers.RD.txt.gz --vcf "gs://.../{chr}.g.vcf.bgz" --id-map id_map.tsv \
#     [--prefix aou_lrGS] [--outdir .] [--flank-frac 0.5] [--samples samples.txt] [--scripts-dir <helpers>]
#   --rd / --vcf / --id-map / --samples: as in prepare_region_data.sh (samples default to column 1 of the ID map)
# Needs the helper scripts next to this one (or --scripts-dir): breakpoint_sd_windows.py, prepare_region_data.sh,
# vcf_to_baf_table.py, extract_gt.sh, subset_rename_samples.py.
set -euo pipefail

usage() { sed -n '2,20p' "$0" | sed 's/^# \{0,1\}//'; exit 1; }
cluster=""; bps=""; bpsd=""; rd=""; vcf=""; idmap=""; samples=""; prefix="gd"; outdir="."; frac=0.5
sdir="$(cd "$(dirname "$0")" && pwd)"
while [[ $# -gt 0 ]]; do
  case $1 in
    --cluster) cluster=$2; shift 2 ;; --breakpoints) bps=$2; shift 2 ;; --bp-sd) bpsd=$2; shift 2 ;;
    --rd) rd=$2; shift 2 ;; --vcf) vcf=$2; shift 2 ;; --id-map) idmap=$2; shift 2 ;; --samples) samples=$2; shift 2 ;;
    --prefix) prefix=$2; shift 2 ;; --outdir) outdir=$2; shift 2 ;; --flank-frac) frac=$2; shift 2 ;;
    --scripts-dir) sdir=$2; shift 2 ;; -h|--help) usage ;; *) echo "unknown option $1" >&2; usage ;;
  esac
done
for v in cluster bps bpsd rd vcf idmap; do [[ -z "${!v}" ]] && { echo "missing --$v" >&2; usage; }; done
mkdir -p "$outdir"
bed="$outdir/$prefix.$cluster.windows.bed"

# SD windows: chrom, win_start (0-based), win_end, name (<cluster>_<BP>), sd_start, sd_end, sd_length, flank
python3 "$sdir/breakpoint_sd_windows.py" --bed "$bpsd" --cluster "$cluster" --breakpoints "$bps" --flank-frac "$frac" \
  --out "$bed"

extra=(); [[ -n "$samples" ]] && extra=(--samples "$samples")
while IFS=$'\t' read -r chrom ws we name _; do
  bp=${name#"${cluster}_"}
  bash "$sdir/prepare_region_data.sh" --chr "$chrom" --start $((ws + 1)) --end "$we" --name "$prefix.${cluster}_$bp" \
    --rd "$rd" --vcf "$vcf" --id-map "$idmap" --outdir "$outdir" ${extra[@]+"${extra[@]}"}
done < "$bed"

# Whole span: from the first window start to the last window end (flank of the first and last SD)
read -r chrom span_start span_end < <(awk -F'\t' 'NR == 1 {c = $1; s = $2} {e = $3} END {print c, s, e}' "$bed")
first=$(head -1 "$bed" | cut -f4); last=$(tail -1 "$bed" | cut -f4)
span="${cluster}_${first#"${cluster}_"}-${last#"${cluster}_"}"
bash "$sdir/prepare_region_data.sh" --chr "$chrom" --start $((span_start + 1)) --end "$span_end" \
  --name "$prefix.$span.region" --rd "$rd" --id-map "$idmap" --outdir "$outdir" --rd-only
printf "%s\t%s\t%s\t%s\n" "$chrom" "$span_start" "$span_end" "$span.region" >> "$bed"
echo "done: windows in $bed; tables $outdir/$prefix.${cluster}_*.renamed.*.gz" >&2
