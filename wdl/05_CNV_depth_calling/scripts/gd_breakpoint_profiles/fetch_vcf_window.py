#!/usr/bin/env python3
"""Fetch one region of a VCF (local or gs://, requester pays) with pysam into a bgzipped, tabix-indexed VCF, robust to
long remote reads: an OAuth access token lasts about an hour, so a large region of a many-sample VCF can outlive it and
the read then fails ("truncated file" / failed to read BGZF block). For gs:// paths the file is reopened with a fresh
token every --refresh-minutes, and after a read error it is reopened with a fresh token and the fetch resumes right
after the last record written (at most --retries read errors in all, waiting 15 s x error count), so no record is
duplicated or lost. Records overlapping the region are written, as with a plain fetch.

Token: VM metadata server, else `gcloud auth print-access-token`, else $GCS_OAUTH_TOKEN as given. Requester-pays
project: --requester-pays-project, else $GCS_REQUESTER_PAYS_PROJECT, else $GOOGLE_PROJECT.

Usage:
  fetch_vcf_window.py --vcf gs://.../chr1.g.vcf.bgz --region chr1:145541830-147568500 --out window.vcf.gz \\
      [--refresh-minutes 40] [--retries 6] [--requester-pays-project PROJECT]
"""
import argparse
import json
import os
import subprocess
import sys
import time
import urllib.request

import pysam

METADATA_TOKEN_URL = "http://metadata.google.internal/computeMetadata/v1/instance/service-accounts/default/token"


def refresh_token():
    """Put a fresh access token in $GCS_OAUTH_TOKEN (htslib reads it when a gs:// file is opened)."""
    try:
        req = urllib.request.Request(METADATA_TOKEN_URL, headers={"Metadata-Flavor": "Google"})
        with urllib.request.urlopen(req, timeout=10) as fh:
            os.environ["GCS_OAUTH_TOKEN"] = json.load(fh)["access_token"]
        return
    except Exception:
        pass
    try:
        os.environ["GCS_OAUTH_TOKEN"] = subprocess.check_output(["gcloud", "auth", "print-access-token"], text=True,
                                                                stderr=subprocess.DEVNULL).strip()
    except Exception:
        pass  # keep whatever $GCS_OAUTH_TOKEN holds


class Reopen(Exception):
    pass


def fetch_window(path, chrom, start, end, out_path, refresh_seconds=2400, retries=6):
    remote = path.startswith("gs://")
    if remote:
        refresh_token()
    vcf = pysam.VariantFile(path)
    if chrom not in vcf.header.contigs:
        sys.exit(f"{chrom} is not a contig of {path}")
    n_samples = len(vcf.header.samples)
    out = pysam.VariantFile(out_path, "wz", header=vcf.header)
    n = 0
    last_pos, at_last = None, 0  # position of the last record written and how many were written at it
    failures = 0
    while True:
        opened = time.time()
        fetch_from = start - 1 if last_pos is None else last_pos - 1
        skip_at_last = at_last
        try:
            for rec in vcf.fetch(chrom, fetch_from, end):
                if last_pos is not None and rec.pos < last_pos:
                    continue  # written before the reopen (records overlapping the resume point)
                if last_pos is not None and rec.pos == last_pos and skip_at_last > 0:
                    skip_at_last -= 1
                    continue
                out.write(rec)
                n += 1
                if rec.pos == last_pos:
                    at_last += 1
                else:
                    last_pos, at_last = rec.pos, 1
                if remote and time.time() - opened > refresh_seconds:
                    raise Reopen()
            break
        except Reopen:
            print(f"  {n} records; reopening with a fresh token at {chrom}:{last_pos}", file=sys.stderr)
        except OSError as e:
            failures += 1
            print(f"  read error after {n} records at {chrom}:{last_pos} ({e}); retry {failures}/{retries}",
                  file=sys.stderr)
            if failures > retries:
                raise
            time.sleep(15 * failures)
        if remote:
            refresh_token()
        vcf.close()
        vcf = pysam.VariantFile(path)
    out.close()
    vcf.close()
    pysam.tabix_index(out_path, preset="vcf", force=True)
    print(f"VCF: {n} records, {n_samples} samples", file=sys.stderr)
    return n


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--vcf", required=True, help="VCF/BCF path (local or gs://) with its index next to it")
    p.add_argument("--region", required=True, help="chrom:start-end (1-based, inclusive)")
    p.add_argument("--out", required=True, help="output .vcf.gz (indexed)")
    p.add_argument("--refresh-minutes", type=float, default=40)
    p.add_argument("--retries", type=int, default=6)
    p.add_argument("--requester-pays-project")
    if len(sys.argv) == 1:
        p.print_help()
        sys.exit(1)
    a = p.parse_args()
    rp = a.requester_pays_project or os.environ.get("GCS_REQUESTER_PAYS_PROJECT") or os.environ.get("GOOGLE_PROJECT")
    if rp:
        os.environ["GCS_REQUESTER_PAYS_PROJECT"] = rp
    chrom, rng = a.region.rsplit(":", 1)
    start, end = (int(x.replace(",", "")) for x in rng.split("-"))
    fetch_window(a.vcf, chrom, start, end, a.out, a.refresh_minutes * 60, a.retries)


if __name__ == "__main__":
    main()
