#!/usr/bin/env bash
# BAF and genotype (GT) tables for named windows of a VCF (local or gs://, read with pysam, no remote tabix needed),
# samples taken from column 1 of an ID map (header skipped, unique, limited to the VCF) and renamed with it.
#
# Usage: extract_breakpoint_baf_gt.sh <vcf> <id_map.tsv> <out_prefix> NAME=chrom:start-end [NAME=chrom:start-end ...]
#   outputs per window: <out_prefix>_<NAME>.vcf.gz (+ .tbi), <out_prefix>_<NAME>.renamed.baf.tsv.gz,
#                       <out_prefix>_<NAME>.renamed.GT.tsv.gz
# Example (16p11.2 breakpoint SDs, each +/- half its length):
#   extract_breakpoint_baf_gt.sh "gs://.../chr16.g.vcf.bgz" id_map.tsv 16p11/16p11 \
#     BP1=chr16:28255730-28543402 BP2=chr16:28482918-28920780 BP3=chr16:29000136-29140307 \
#     BP4=chr16:29166926-29795926 BP5=chr16:30111068-30420922
# For gs:// the token comes from $GCS_OAUTH_TOKEN or `gcloud auth print-access-token`, the requester-pays project from
# $GCS_REQUESTER_PAYS_PROJECT or $GOOGLE_PROJECT. Needs python3 with pysam and, next to this script:
# fetch_vcf_regions.py, vcf_to_baf_table.py, extract_gt.py, subset_rename_samples.py.
set -euo pipefail
if [[ $# -lt 4 || "$1" == "-h" || "$1" == "--help" ]]; then sed -n '2,15p' "$0" | sed 's/^# \{0,1\}//'; exit 1; fi
vcf=$1; idmap=$2; prefix=$3; shift 3
sdir="$(cd "$(dirname "$0")" && pwd)"
mkdir -p "$(dirname "$prefix")"

# 1. Fetch every window from the VCF (one pass over the remote file per window)
regions=(); for w in "$@"; do regions+=(--region "$w"); done
python3 "$sdir/fetch_vcf_regions.py" --vcf "$vcf" --out-prefix "$prefix" "${regions[@]}"

# 2. Samples: column 1 of the ID map, header skipped, unique
samples="$prefix.samples.txt"
tail -n +2 "$idmap" | awk 'NF && !s[$1]++ {print $1}' > "$samples"

# 3. BAF and GT per window, then rename
for w in "$@"; do
  name=${w%%=*}; region=${w#*=}; p="${prefix}_$name"
  python3 "$sdir/vcf_to_baf_table.py" "$p.vcf.gz" "$samples" "$p.baf.tsv.gz" --region "$region"
  python3 "$sdir/extract_gt.py" --vcf "$p.vcf.gz" --samples "$samples" --out "$p.GT.tsv.gz" --region "$region"
  for t in baf GT; do
    python3 "$sdir/subset_rename_samples.py" --input "$p.$t.tsv.gz" --id-map "$idmap" --out "$p.renamed.$t.tsv.gz" \
      --fixed-cols 5 --allow-missing
    rm -f "$p.$t.tsv.gz"
  done
done
echo "done: ${prefix}_<NAME>.renamed.{baf,GT}.tsv.gz for $*" >&2
