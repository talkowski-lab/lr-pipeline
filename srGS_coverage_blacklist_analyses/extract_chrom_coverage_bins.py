#!/usr/bin/env python3
"""Localize only the Hail Table partitions covering a single chromosome from a
gnomAD-style coverage Hail Table (gs://.../*.ht), bin into fixed-size windows,
and export per-bin coverage summaries as a BED file.

Loci absent from the table are treated as zero coverage in every sample, so
every bin across the full chromosome length is emitted (matching mosdepth's
per-bin output used for lrGS). Output columns (no header):
  chrom, start, end,
  mean_cov       - sum of per-locus cohort mean depth / bin length
  n_loci         - number of loci in the bin present in the table
  low<X>_frac    - per-locus fraction of samples NOT over_<X> (1 - over_<X>),
                   averaged over the bin (absent loci count as 1.0); one
                   column per --low-thresholds value, in the order given

Uses the table's own rows/metadata.json.gz (_jRangeBounds) to identify which
partitions overlap the requested chromosome, so only those partitions (plus
the small globals/references/index components) are downloaded -- not the
whole table.
"""
import argparse
import gzip
import json
import os
from concurrent.futures import ThreadPoolExecutor

import hail as hl
from google.cloud import storage


def parse_gs_path(gs_path):
    rest = gs_path[len("gs://"):]
    bucket, _, prefix = rest.partition("/")
    return bucket, prefix.rstrip("/")


def download_blob(bucket, blob_name, local_path):
    os.makedirs(os.path.dirname(local_path), exist_ok=True)
    bucket.blob(blob_name).download_to_filename(local_path)


def download_prefix(client, bucket, bucket_name, table_prefix, sub_dir, local_root):
    for blob in client.list_blobs(bucket_name, prefix=f"{table_prefix}/{sub_dir}/"):
        rel = blob.name[len(table_prefix) + 1:]
        download_blob(bucket, blob.name, os.path.join(local_root, rel))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--hail-table-path", required=True)
    ap.add_argument("--chrom", required=True)
    ap.add_argument("--bin-size", type=int, default=100)
    ap.add_argument("--out-bed", required=True)
    ap.add_argument("--work-dir", default="ht_local")
    ap.add_argument("--reference-genome", default="GRCh38")
    ap.add_argument("--low-thresholds", default="5,10",
                    help="comma-separated over_<X> fields to report as low<X>_frac")
    ap.add_argument("--driver-memory", default="12g")
    ap.add_argument("--download-threads", type=int, default=16)
    args = ap.parse_args()

    client = storage.Client.create_anonymous_client()
    bucket_name, table_prefix = parse_gs_path(args.hail_table_path)
    bucket = client.bucket(bucket_name)
    local_root = os.path.join(args.work_dir, os.path.basename(table_prefix))

    download_blob(bucket, f"{table_prefix}/metadata.json.gz", os.path.join(local_root, "metadata.json.gz"))
    download_prefix(client, bucket, bucket_name, table_prefix, "globals", local_root)
    download_prefix(client, bucket, bucket_name, table_prefix, "references", local_root)

    rows_meta_path = os.path.join(local_root, "rows", "metadata.json.gz")
    download_blob(bucket, f"{table_prefix}/rows/metadata.json.gz", rows_meta_path)
    with gzip.open(rows_meta_path) as f:
        rows_meta = json.load(f)

    bounds = rows_meta["_jRangeBounds"]
    partfiles = rows_meta["_partFiles"]
    idxs = [
        i for i, b in enumerate(bounds)
        if b["start"]["locus"]["contig"] == args.chrom or b["end"]["locus"]["contig"] == args.chrom
    ]
    if not idxs:
        raise SystemExit(f"No partitions found overlapping {args.chrom}")
    print(f"{args.chrom}: {len(idxs)} partitions (index {min(idxs)}-{max(idxs)})")

    def fetch_partition(i):
        pf = partfiles[i]
        download_blob(bucket, f"{table_prefix}/rows/parts/{pf}", os.path.join(local_root, "rows", "parts", pf))
        for blob in client.list_blobs(bucket_name, prefix=f"{table_prefix}/index/{pf}.idx/"):
            rel = blob.name[len(table_prefix) + 1:]
            download_blob(bucket, blob.name, os.path.join(local_root, rel))

    with ThreadPoolExecutor(max_workers=args.download_threads) as pool:
        list(pool.map(fetch_partition, idxs))

    thresholds = [int(x) for x in args.low_thresholds.split(",")]
    hl.init(spark_conf={"spark.driver.memory": args.driver_memory})
    chrom_len = hl.get_reference(args.reference_genome).lengths[args.chrom]
    n_bins = (chrom_len + args.bin_size - 1) // args.bin_size

    ht = hl.read_table(local_root)
    ht = ht.filter(ht.locus.contig == args.chrom)
    ht = ht.annotate(bin_start=(ht.locus.position - 1) // args.bin_size * args.bin_size)
    agg_exprs = {"sum_mean": hl.agg.sum(ht.mean), "n_loci": hl.agg.count()}
    for x in thresholds:
        agg_exprs[f"sum_low{x}"] = hl.agg.sum(1.0 - ht[f"over_{x}"])
    grouped = ht.group_by(bin_start=ht.bin_start).aggregate(**agg_exprs)

    bins = hl.utils.range_table(n_bins)
    bins = bins.annotate(bin_start=hl.int32(bins.idx * args.bin_size)).key_by("bin_start")
    bins = bins.annotate(end=hl.min(bins.bin_start + args.bin_size, chrom_len))
    bins = bins.annotate(**grouped[bins.bin_start])
    bin_len = bins.end - bins.bin_start
    n_loci = hl.or_else(bins.n_loci, 0)
    out_exprs = {
        "chrom": args.chrom,
        "start": bins.bin_start,
        "end": bins.end,
        "mean_cov": hl.or_else(bins.sum_mean, 0.0) / bin_len,
        "n_loci": n_loci,
    }
    for x in thresholds:
        out_exprs[f"low{x}_frac"] = (hl.or_else(bins[f"sum_low{x}"], 0.0) + hl.float64(bin_len - n_loci)) / bin_len
    bins = bins.key_by().select(**out_exprs)
    bins.export(args.out_bed, header=False)


if __name__ == "__main__":
    main()
