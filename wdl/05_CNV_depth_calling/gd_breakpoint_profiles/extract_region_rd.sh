#!/usr/bin/env bash
# Depth (RD) matrix across one region, header kept, samples kept and renamed with an ID map.
#
# Usage: extract_region_rd.sh <rd_matrix.txt.gz> <chrom:start-end> <id_map.tsv> <out.renamed.RD.txt.gz>
#   rd_matrix: bgzipped + tabix-indexed '#Chr Start End <samples>' (local); id_map: real ID, new ID (header skipped)
# Example (16p11.2 BP1-BP5 with flanks):
#   extract_region_rd.sh aou_lrGS_GD_carriers.RD.txt.gz chr16:28255730-30420922 id_map.tsv 16p11/16p11_BP1-5.region.renamed.RD.txt.gz
# Needs tabix, bgzip and subset_rename_samples.py next to this script.
set -euo pipefail
if [[ $# -ne 4 || "$1" == "-h" || "$1" == "--help" ]]; then sed -n '2,9p' "$0" | sed 's/^# \{0,1\}//'; exit 1; fi
rd=$1; region=$2; idmap=$3; out=$4
sdir="$(cd "$(dirname "$0")" && pwd)"
mkdir -p "$(dirname "$out")"
raw="${out%.renamed.RD.txt.gz}.RD.txt.gz"; [[ "$raw" == "$out" ]] && raw="$out.raw.gz"
{ gzip -cd "$rd" | head -1 || true; tabix "$rd" "$region"; } | bgzip -c > "$raw"
echo "RD: $(( $(gzip -cd "$raw" | wc -l) - 1 )) bins in $region" >&2
python3 "$sdir/subset_rename_samples.py" --input "$raw" --id-map "$idmap" --out "$out" --fixed-cols 3 --allow-missing
rm -f "$raw"
echo "done: $out" >&2
