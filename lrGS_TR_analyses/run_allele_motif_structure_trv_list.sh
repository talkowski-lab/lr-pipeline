#!/usr/bin/env bash
# Run plot_allele_motif_structure.py for a list of TRVs (no SR data), then merge the pages in list order.
# Usage: run_allele_motif_structure_trv_list.sh <list.tsv: TRV_ID label> <main TRV table> <genotypes.tsv> <ancestry.tsv[.gz]>
#        <second label> <second TRV table> <out_dir> <out_tag>
#   The second cohort's record with the same TRV ID (if present) is shown in the top histogram.
# Writes <out_dir>/<label>.<out_tag>.{pdf,alleles.tsv,motifs.tsv,log} and <out_dir>/all.<out_tag>.pdf.
set -euo pipefail

if [ $# -ne 8 ]; then
    sed -n 2,7p "$0" >&2
    exit 1
fi
list=$1; main_trv=$2; genotypes=$3; ancestry=$4; second_label=$5; second_trv=$6; out_dir=$7; tag=$8
here=$(cd "$(dirname "$0")" && pwd)
mkdir -p "$out_dir"

pdfs=()
while IFS=$'\t' read -r id label; do
    if python "$here/plot_allele_motif_structure.py" --trv "$main_trv" --trv-id "$id" --label "$label" \
        --genotypes "$genotypes" --ancestry "$ancestry" --top-trv "$second_label=$second_trv=$id" \
        --out-prefix "$out_dir/$label.$tag" > "$out_dir/$label.$tag.log" 2>&1; then
        pdfs+=("$out_dir/$label.$tag.pdf")
    else
        echo "FAILED: $id $label (see $out_dir/$label.$tag.log)"
    fi
done < "$list"

python - "$out_dir/all.$tag.pdf" "${pdfs[@]}" <<'PY'
import sys
from pypdf import PdfWriter
w = PdfWriter()
for f in sys.argv[2:]:
    w.append(f)
w.write(sys.argv[1])
print(f"{len(sys.argv) - 2} TRVs -> {sys.argv[1]}")
PY
