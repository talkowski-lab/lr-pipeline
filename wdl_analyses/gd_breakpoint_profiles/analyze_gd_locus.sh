#!/usr/bin/env bash
# Breakpoint profile of one genomic-disorder (GD) locus from the prepared tables of its two breakpoint-SD windows
# (prepare_region_data.sh outputs): BAF loss scan, LD and one figure with the left SD on the left and the right SD on the
# right (SD tracks, copy number, folded BAF, het-only BAF, LD; DEL carriers red to orange, DUP carriers blue to green,
# arrows at the BAF breakpoint estimates of DEL carriers).
#
# For the left window the deletion lies to the right (find_baf_terminal_loh.py --side right); for the right window it
# lies to the left (--side left). LD is computed without the locus carriers (--highlight-pattern).
#
# Usage:
#   analyze_gd_locus.sh --name 16p11.2-BP4-5 --highlight-pattern "GD_16p11\.2\-BP4\-5_(DEL|DUP)_chr16_29638676_30188531" \
#     --left-prefix 16p11.2-BP4-5.left --left-title "16p11:BP4" --left-window 29166925,29795926 \
#     --right-prefix 16p11.2-BP4-5.right --right-title "16p11:BP5" --right-window 30111067,30420922 \
#     --medians GD_samples.RD_median --segdup hg38.SegDup.sorted.merged.bed.gz --bp-sd GD_NAHR.breakpoint_SDs.bed.gz \
#     [--gd-start 29638676 --gd-end 30188531] [--outdir .] [--scripts-dir <helper scripts; default: this folder>]
#   --left-prefix/--right-prefix: path prefix of <prefix>.renamed.{RD.txt,baf.tsv,GT.tsv}.gz
#   --left-window/--right-window: start,end (bp) of each plotted window; --gd-start/--gd-end: dashed lines at the GD edges
#   --overview-rd: optional depth matrix across the whole locus (left window start to right window end), drawn as a
#   full-width copy-number panel with SD tracks on top of the figure
#   --cn-bin: copy-number bin size (bp) for all depth panels (default: automatic per window)
# Outputs (in --outdir, prefix <name>): .breakpoint_profile.pdf, .left_loh.ranked.tsv, .right_loh.ranked.tsv,
#   .left.ld.{sites,binned}.tsv, .right.ld.{sites,binned}.tsv, .{left,right}_loh.clean_sites.baf.tsv.gz
# Needs python3 (numpy, pandas), R (ggplot2, data.table, optparse, patchwork) and the helper scripts
# find_baf_terminal_loh.py, compute_ld.py, plot_two_regions.R.
set -euo pipefail

usage() { sed -n '2,24p' "$0" | sed 's/^# \{0,1\}//'; exit 1; }
name=""; pat=""; lp=""; lt="left SD"; lw=""; rp=""; rt="right SD"; rw=""; med=""; sd=""; bpsd=""; gs=""; ge=""
ovrd=""; cnbin=0; outdir="."; sdir="$(cd "$(dirname "$0")" && pwd)"
while [[ $# -gt 0 ]]; do
  case $1 in
    --name) name=$2; shift 2 ;; --highlight-pattern) pat=$2; shift 2 ;;
    --left-prefix) lp=$2; shift 2 ;; --left-title) lt=$2; shift 2 ;; --left-window) lw=$2; shift 2 ;;
    --right-prefix) rp=$2; shift 2 ;; --right-title) rt=$2; shift 2 ;; --right-window) rw=$2; shift 2 ;;
    --medians) med=$2; shift 2 ;; --segdup) sd=$2; shift 2 ;; --bp-sd) bpsd=$2; shift 2 ;;
    --gd-start) gs=$2; shift 2 ;; --gd-end) ge=$2; shift 2 ;; --overview-rd) ovrd=$2; shift 2 ;; --cn-bin) cnbin=$2; shift 2 ;;
    --outdir) outdir=$2; shift 2 ;; --scripts-dir) sdir=$2; shift 2 ;; -h|--help) usage ;;
    *) echo "unknown option $1" >&2; usage ;;
  esac
done
for v in name pat lp lw rp rw med sd bpsd; do [[ -z "${!v}" ]] && { echo "missing a required option ($v)" >&2; usage; }; done
mkdir -p "$outdir"
o="$outdir/$name"

# BAF loss scan on clean biallelic sites; the deletion runs right of the left SD and left of the right SD
python3 "$sdir/find_baf_terminal_loh.py" --baf "$lp.renamed.baf.tsv.gz" --out-prefix "$o.left_loh" --side right --top 15
python3 "$sdir/find_baf_terminal_loh.py" --baf "$rp.renamed.baf.tsv.gz" --out-prefix "$o.right_loh" --side left --top 15

# LD without the locus carriers; window size scales with the window (about 100 windows per region, at least 5 kb)
ld_bin() { awk -F, '{b = int(($2 - $1) / 100 / 1000 + 0.5) * 1000; print (b < 5000 ? 5000 : b)}' <<< "$1"; }
python3 "$sdir/compute_ld.py" --gt "$lp.renamed.GT.tsv.gz" --out-prefix "$o.left.ld" --exclude-pattern "$pat" \
  --bin-size "$(ld_bin "$lw")"
python3 "$sdir/compute_ld.py" --gt "$rp.renamed.GT.tsv.gz" --out-prefix "$o.right.ld" --exclude-pattern "$pat" \
  --bin-size "$(ld_bin "$rw")"
rm -f "$o.left.ld.pairs.tsv.gz" "$o.right.ld.pairs.tsv.gz"

# Figure: dashed lines at the GD edges
mark1=(); mark2=()
[[ -n "$gs" ]] && mark1=(--mark1 "$(awk -v x="$gs" 'BEGIN{print x / 1e6}')")
[[ -n "$ge" ]] && mark2=(--mark2 "$(awk -v x="$ge" 'BEGIN{print x / 1e6}')")
overview=()
[[ -n "$ovrd" ]] && overview=(--overview-rd "$ovrd" --overview-window "${lw%%,*},${rw##*,}")
Rscript "$sdir/plot_two_regions.R" \
  --baf1 "$o.left_loh.clean_sites.baf.tsv.gz" --rd1 "$lp.renamed.RD.txt.gz" --ld1 "$o.left.ld" --bp1 "$o.left_loh.ranked.tsv" \
  --title1 "left: $lt" --win1 "$lw" ${mark1[@]+"${mark1[@]}"} \
  --baf2 "$o.right_loh.clean_sites.baf.tsv.gz" --rd2 "$rp.renamed.RD.txt.gz" --ld2 "$o.right.ld" --bp2 "$o.right_loh.ranked.tsv" \
  --title2 "right: $rt" --win2 "$rw" ${mark2[@]+"${mark2[@]}"} ${overview[@]+"${overview[@]}"} \
  --bin-size "$cnbin" --overview-bin "$cnbin" \
  --medians "$med" --highlight-pattern "$pat" --del-pattern "_DEL_" --dup-pattern "_DUP_" --height 22 \
  --segdup "SD=$sd,SD breakpoint=$bpsd" --out "$o.breakpoint_profile.pdf" --title "$name"
echo "done: $o.breakpoint_profile.pdf" >&2
