#!/bin/bash
# Extract INFO fields for insertions >= MIN_LEN bp from an annotated gnomAD-LR VCF (local or gs://).
usage() {
    echo "Usage: $0 <vcf> <min_len> <out.tsv.gz>"
    echo "  vcf      annotated VCF (bgzipped + tbi); gs:// needs GCS_OAUTH_TOKEN=\$(gcloud auth print-access-token)"
    echo "  min_len  minimum |SVLEN| (use the same value as the remap run's min_length)"
    echo "  out      output TSV (gzipped, with header)"
}
if [ "${1:-}" = "-h" ] || [ "${1:-}" = "--help" ] || [ $# -ne 3 ]; then usage; exit 0; fi
set -euo pipefail

VCF=$1
MIN_LEN=$2
OUT=$3

{
    echo -e "ID\tFILTER\tSVLEN\tAF\tAC\tSOURCE\tREGION\tTR_ENVELOPED\tSUB_FAMILY\tORIGIN\tallele_subtype\tTRID\tgnomAD_V4_match_type"
    bcftools query \
        -i "INFO/allele_type=\"ins\" && abs(INFO/SVLEN)>=${MIN_LEN}" \
        -f '%ID\t%FILTER\t%INFO/SVLEN\t%INFO/AF\t%INFO/AC\t%INFO/SOURCE\t%INFO/REGION\t%INFO/TR_ENVELOPED\t%INFO/SUB_FAMILY\t%INFO/ORIGIN\t%INFO/allele_subtype\t%INFO/TRID\t%INFO/gnomAD_V4_match_type\n' \
        "${VCF}"
} | gzip > "${OUT}"
