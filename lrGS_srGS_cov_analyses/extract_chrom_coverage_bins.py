#!/usr/bin/env python3
"""Localize only the Hail Table partitions covering a single chromosome from a
gnomAD-style coverage Hail Table (gs://.../*.ht), bin into fixed-size windows,
and export the mean of the 'mean' column per bin as a BED file.

Uses the table's own rows/metadata.json.gz (_jRangeBounds) to identify which
partitions overlap the requested chromosome, so only those partitions (plus
the small globals/references/index components) are downloaded -- not the
whole table.
"""
import argparse
import gzip
import json
import os

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

    for i in idxs:
        pf = partfiles[i]
        download_blob(bucket, f"{table_prefix}/rows/parts/{pf}", os.path.join(local_root, "rows", "parts", pf))
        for blob in client.list_blobs(bucket_name, prefix=f"{table_prefix}/index/{pf}.idx/"):
            rel = blob.name[len(table_prefix) + 1:]
            download_blob(bucket, blob.name, os.path.join(local_root, rel))

    hl.init()
    ht = hl.read_table(local_root)
    ht = ht.filter(ht.locus.contig == args.chrom)
    ht = ht.annotate(bin_start=(ht.locus.position - 1) // args.bin_size * args.bin_size)
    grouped = ht.group_by(contig=ht.locus.contig, bin_start=ht.bin_start).aggregate(
        mean_cov=hl.agg.mean(ht.mean)
    )
    grouped = grouped.annotate(
        locus=hl.locus(grouped.contig, grouped.bin_start + 1, reference_genome=args.reference_genome)
    )
    grouped = grouped.key_by("locus")
    grouped = grouped.select(
        chrom=grouped.contig,
        start=grouped.bin_start,
        end=grouped.bin_start + args.bin_size,
        mean_cov=grouped.mean_cov,
    )
    grouped = grouped.key_by().drop("locus")
    grouped.export(args.out_bed, header=False)


if __name__ == "__main__":
    main()
