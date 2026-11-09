#!/usr/bin/env bash
# Extract genotypes (GT) of selected samples from a VCF/BCF into a table: CHROM POS ID REF ALT + one GT column per sample.
#
# Usage: extract_gt.sh -v <vcf> -s <samples> -o <out.tsv[.gz]> [-r <region>] [-n]
#   -v  input VCF/BCF (indexed if -r is used)
#   -s  samples: comma-separated IDs, or a file with one ID per line (first column used)
#   -o  output table; gzipped when it ends in .gz; '-' for stdout
#   -r  region, e.g. chr16:29300000-29700000 (optional)
#   -n  keep only sites where at least one selected sample is non-reference (GT has an ALT allele)
# Needs bcftools on PATH. Samples missing from the VCF stop the run.
set -euo pipefail

usage() { sed -n '2,11p' "$0" | sed 's/^# \{0,1\}//'; exit 1; }
vcf=""; samples=""; out=""; region=""; nonref=0
while getopts "v:s:o:r:nh" o; do
  case $o in
    v) vcf=$OPTARG ;; s) samples=$OPTARG ;; o) out=$OPTARG ;; r) region=$OPTARG ;; n) nonref=1 ;; *) usage ;;
  esac
done
[[ -z "$vcf" || -z "$samples" || -z "$out" ]] && usage
command -v bcftools >/dev/null || { echo "bcftools not found on PATH" >&2; exit 1; }

tmp=$(mktemp -d); trap 'rm -rf "$tmp"' EXIT
# unique IDs, first occurrence order (bcftools rejects a sample listed twice)
if [[ -f "$samples" ]]; then awk 'NF && !seen[$1]++ {print $1}' "$samples" > "$tmp/s.txt"; else tr ',' '\n' <<< "$samples" | awk 'NF && !seen[$0]++' > "$tmp/s.txt"; fi
missing=$(grep -vxFf <(bcftools query -l "$vcf") "$tmp/s.txt" || true)
[[ -n "$missing" ]] && { echo "not in the VCF: $(paste -sd, - <<< "$missing")" >&2; exit 1; }

args=(-S "$tmp/s.txt" -f '%CHROM\t%POS\t%ID\t%REF\t%ALT[\t%GT]\n')
[[ -n "$region" ]] && args+=(-r "$region")
[[ $nonref -eq 1 ]] && args+=(-i 'GT="alt"')
{ printf 'CHROM\tPOS\tID\tREF\tALT\t%s\n' "$(paste -sd'\t' "$tmp/s.txt")"
  bcftools query "${args[@]}" "$vcf"
} | if [[ "$out" == "-" ]]; then cat; elif [[ "$out" == *.gz ]]; then gzip -c > "$out"; else cat > "$out"; fi
[[ "$out" != "-" ]] && echo "wrote $out ($(wc -l < "$tmp/s.txt" | tr -d ' ') samples)" >&2
exit 0
