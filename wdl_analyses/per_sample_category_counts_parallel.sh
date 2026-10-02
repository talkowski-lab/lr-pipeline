#!/bin/bash
# Splits one chromosome-shard VCF into N disjoint genomic-position chunks
# (using bcftools view -r, which is tabix-index-based -- a true seek, not a
# linear scan -- so this is genuine data-parallelism: N workers each handle
# 1/N of the data, rather than N workers each rescanning 100% of it).
# Each chunk is processed independently and concurrently by the existing
# per_sample_category_counts.py, then the chunk results are
# summed/unioned back together with concat_sample_category_counts.py --
# correct because disjoint position ranges can never double-count a variant
# or a sample's genotype at it.
#
# (An earlier version of this ran the 9 count categories + 3 gene-list
# categories as 9-12 concurrent `bcftools view -i <category filter>` passes
# instead. That measured SLOWER in practice (~167s vs ~66s single-pass on a
# real chr22 test) because each category's `-i` filter still requires a full
# linear scan/decompression of the entire file -- there's no index to skip
# non-matching INFO-field values, so parallelizing across categories just
# multiplies total I/O instead of dividing it. Position-range splitting
# avoids this because -r IS index-seekable.)
#
# Usage: per_sample_category_counts_parallel.sh <in.vcf.gz> <out_prefix> <n_chunks> <per_sample_category_counts.py> <concat_sample_category_counts.py> [category]
# category (default: all) is passed through to per_sample_category_counts.py.
set -euo pipefail

VCF=$1
OUT_PREFIX=$2
N_CHUNKS=$3
PER_SAMPLE_SCRIPT=$4
CONCAT_SCRIPT=$5
CATEGORY=${6:-all}

WORKDIR=$(mktemp -d)
trap 'rm -rf "$WORKDIR"' EXIT

# tabix -l lists only the contig(s) actually present in the index -- an
# index lookup, not a data scan -- unlike `bcftools view -H | head -1`,
# which is both wasteful (reads/decompresses from the start of the file)
# and unsafe under `set -o pipefail` (head's early exit sends SIGPIPE
# upstream, which pipefail then reports as the pipeline's exit status).
CONTIG=$(tabix -l "$VCF" | awk 'NR==1')
LENGTH=$(bcftools view -h "$VCF" | grep -m1 "^##contig=<ID=${CONTIG},length=" | sed -E 's/.*length=([0-9]+).*/\1/')

if [ -z "$CONTIG" ] || [ -z "$LENGTH" ]; then
    echo "ERROR: could not determine contig/length from $VCF" >&2
    exit 1
fi

CHUNK_SIZE=$(( (LENGTH + N_CHUNKS - 1) / N_CHUNKS ))

# --- slice into N disjoint position-range chunks, in parallel (index-seekable) ---
for (( i=0; i<N_CHUNKS; i++ )); do
    start=$(( i * CHUNK_SIZE + 1 ))
    end=$(( (i + 1) * CHUNK_SIZE ))
    if [ "$end" -gt "$LENGTH" ]; then end=$LENGTH; fi
    (
        bcftools view -r "${CONTIG}:${start}-${end}" --regions-overlap pos "$VCF" -Ob -o "$WORKDIR/chunk_${i}.bcf"
    ) &
done
wait

# --- process each chunk concurrently with the existing per-sample script ---
for (( i=0; i<N_CHUNKS; i++ )); do
    (
        python3 "$PER_SAMPLE_SCRIPT" "$WORKDIR/chunk_${i}.bcf" "$WORKDIR/chunk_${i}" "$CATEGORY"
    ) &
done
wait

# --- sum counts / union gene sets across chunks (disjoint, so this is exact) ---
python3 "$CONCAT_SCRIPT" "${OUT_PREFIX}.category_counts.tsv" "$WORKDIR"/chunk_*.category_counts.tsv

echo "Wrote ${OUT_PREFIX}.category_counts.tsv"
