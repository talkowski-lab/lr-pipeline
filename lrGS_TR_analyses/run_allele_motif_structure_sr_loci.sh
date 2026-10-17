#!/usr/bin/env bash
# Run plot_allele_motif_structure.py for every gnomAD SR STR locus with a matched long-read TRV, then merge the PDFs.
# Main panels: the matched TRV in the long-read table (per-sample genotypes for the ancestry panel); top panel: gnomAD SR and,
# when it has a matched TRV, the second cohort.
# Usage: run_allele_motif_structure_sr_loci.sh <main loci.tsv> <main TRV table> <genotypes.tsv> <ancestry.tsv[.gz]>
#        <SR allele_sizes.tsv> <second loci.tsv> <second label> <second TRV table> <out_dir> <out_tag>
#   loci.tsv = split_sr_matched_strs.py .NAME.loci.tsv (columns id, match_flag, TRV); loci with match_flag no_TRV are skipped.
# Writes <out_dir>/<id>.<out_tag>.{pdf,alleles.tsv,motifs.tsv} and <out_dir>/all_SR_loci.<out_tag>.pdf (one page per locus).
set -euo pipefail

if [ $# -ne 10 ]; then
    sed -n 2,9p "$0" >&2
    exit 1
fi
main_loci=$1; main_trv=$2; genotypes=$3; ancestry=$4; sr_sizes=$5
second_loci=$6; second_label=$7; second_trv=$8; out_dir=$9; tag=${10}
here=$(cd "$(dirname "$0")" && pwd)
mkdir -p "$out_dir"

awk -F'\t' 'NR==1{for(i=1;i<=NF;i++) c[$i]=i; next} $c["match_flag"]!="no_TRV" && $c["TRV"]!="" {print $c["id"]"\t"$c["TRV"]}' \
    "$main_loci" | sort > "$out_dir/.main_ids"
awk -F'\t' 'NR==1{for(i=1;i<=NF;i++) c[$i]=i; next} $c["match_flag"]!="no_TRV" && $c["TRV"]!="" {print $c["id"]"\t"$c["TRV"]}' \
    "$second_loci" | sort > "$out_dir/.second_ids"

pdfs=()
while IFS=$'\t' read -r id trv; do
    top=()
    second=$(awk -F'\t' -v i="$id" '$1==i{print $2}' "$out_dir/.second_ids")
    if [ -n "$second" ]; then
        top=(--top-trv "$second_label=$second_trv=$second")
    fi
    python "$here/plot_allele_motif_structure.py" --trv "$main_trv" --trv-id "$trv" --label "$id" \
        --genotypes "$genotypes" --ancestry "$ancestry" --sr-allele-sizes "$sr_sizes" --sr-id "$id" "${top[@]}" \
        --out-prefix "$out_dir/$id.$tag" > "$out_dir/$id.$tag.log" 2>&1 || { echo "FAILED: $id (see $out_dir/$id.$tag.log)"; continue; }
    pdfs+=("$out_dir/$id.$tag.pdf")
done < "$out_dir/.main_ids"

python - "$out_dir/all_SR_loci.$tag.pdf" "${pdfs[@]}" <<'PY'
import sys
from pypdf import PdfWriter
w = PdfWriter()
for f in sys.argv[2:]:
    w.append(f)
w.write(sys.argv[1])
print(f"{len(sys.argv) - 2} loci -> {sys.argv[1]}")
PY
rm -f "$out_dir/.main_ids" "$out_dir/.second_ids"
