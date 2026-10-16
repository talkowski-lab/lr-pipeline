#!/usr/bin/env bash
# Local driver for analyze_TRV_vcf.py (same steps as the AnalyzeTRVariants WDL).
# For each VCF in <vcf_list> (one per contig; gs:// needs GCS_OAUTH_TOKEN): extract TRVs with bcftools, dump the
# query table, run `analyze_TRV_vcf.py analyze`, then merge the per-contig per-sample summaries.
# Usage: run_analyze_TRV_vcf.sh <vcf_list> <gencode.gtf.gz> <out_dir>
set -euo pipefail

if [ $# -ne 3 ]; then
    echo "Usage: $0 <vcf_list> <gencode.gtf.gz> <out_dir>" >&2
    exit 1
fi
vcf_list=$1
gtf=$2
out_dir=$3
script_dir=$(cd "$(dirname "$0")" && pwd)
query_format='%CHROM\t%POS\t%ID\t%REF\t%ALT\t%FILTER\t%INFO/TRID\t%INFO/MOTIFS\t%INFO/AC\t%INFO/AF\t%INFO/AN\t%INFO/MC_allele[\t%GT]\n'

mkdir -p "$out_dir"
summaries=()
while read -r vcf; do
    [ -z "$vcf" ] && continue
    prefix=$out_dir/$(basename "$vcf" | sed -E 's/\.vcf(\.gz|\.bgz)?$//')
    bcftools view -i 'INFO/allele_type="trv"' -Oz -o "$prefix.TRV.vcf.gz" "$vcf"
    bcftools index -t -f "$prefix.TRV.vcf.gz"
    bcftools query -l "$prefix.TRV.vcf.gz" > "$prefix.samples.txt"
    bcftools query -f "$query_format" "$prefix.TRV.vcf.gz" | gzip > "$prefix.TRV.query.tsv.gz"
    python3 "$script_dir/analyze_TRV_vcf.py" analyze \
        --query-tsv "$prefix.TRV.query.tsv.gz" \
        --samples "$prefix.samples.txt" \
        --gtf "$gtf" \
        --prefix "$prefix"
    summaries+=("$prefix.TRV.per_sample_summary.tsv")
done < "$vcf_list"

python3 "$script_dir/analyze_TRV_vcf.py" merge-summaries --summaries "${summaries[@]}" --out "$out_dir/TRV.per_sample_summary.tsv"
