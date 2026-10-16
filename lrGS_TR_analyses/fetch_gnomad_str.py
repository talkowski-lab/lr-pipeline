#!/usr/bin/env python3
"""Download the gnomAD short-read STR loci and their allele-size distributions from the gnomAD browser GraphQL API.

For every locus in short_tandem_repeats(dataset): id, gene, diseases, main reference region (as returned by the API;
start is 0-based), reference repeat unit, and the allele-size distribution summed over all cohorts (ancestry, sex,
repeat-unit variant, quality): allele count per repeat-unit count.
Outputs (prefix --out-prefix):
  .loci.tsv          id gene chrom start stop reference_repeat_unit unit_len ref_units diseases n_alleles
  .allele_sizes.tsv  id repunit repunit_count allele_count (summed over ancestry / sex / quality cohorts)
  .raw.json          raw API responses

Usage:
  fetch_gnomad_str.py --dataset gnomad_r4 --out-prefix gnomad_srSTR/gnomad_r4_STR
"""
import argparse
import collections
import json
import time
import urllib.error
import urllib.request

import pandas as pd

API = "https://gnomad.broadinstitute.org/api"
LIST_Q = ("{ short_tandem_repeats(dataset: %s) { id gene { symbol } associated_diseases { name } "
          "main_reference_region { chrom start stop } reference_repeat_unit } }")
DETAIL_Q = ('{ short_tandem_repeat(id: "%s", dataset: %s) { id allele_size_distribution { ancestry_group sex repunit '
            "quality_description distribution { repunit_count frequency } } } }")


def gql(query, retries=8):
    """POST a GraphQL query; on HTTP 429 / 5xx wait (Retry-After or exponential backoff) and retry."""
    req = urllib.request.Request(API, data=json.dumps({"query": query}).encode(), headers={"Content-Type": "application/json"})
    for attempt in range(retries):
        try:
            with urllib.request.urlopen(req, timeout=60) as r:
                d = json.load(r)
            break
        except urllib.error.HTTPError as e:
            if e.code not in (429, 500, 502, 503, 504) or attempt == retries - 1:
                raise
            wait = float(e.headers.get("Retry-After") or 0) or min(5 * 2 ** attempt, 120)
            print(f"HTTP {e.code}; retrying in {wait:.0f} s")
            time.sleep(wait)
    if "errors" in d:
        raise SystemExit(f"gnomAD API error: {d['errors']}")
    return d["data"]


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--dataset", default="gnomad_r4")
    p.add_argument("--out-prefix", required=True)
    p.add_argument("--sleep", type=float, default=2.0, help="Seconds between API calls (default 2)")
    args = p.parse_args()

    loci = gql(LIST_Q % args.dataset)["short_tandem_repeats"]
    raw, rows, sizes = {"list": loci, "details": {}}, [], []
    for x in loci:
        det = gql(DETAIL_Q % (x["id"], args.dataset))["short_tandem_repeat"]
        raw["details"][x["id"]] = det
        counts = collections.Counter()
        for cohort in det["allele_size_distribution"]:
            for item in cohort["distribution"]:
                counts[(cohort["repunit"], item["repunit_count"])] += item["frequency"]
        for (ru, n), c in sorted(counts.items()):
            sizes.append({"id": x["id"], "repunit": ru, "repunit_count": n, "allele_count": c})
        reg, unit = x["main_reference_region"], x["reference_repeat_unit"]
        rows.append({"id": x["id"], "gene": x["gene"]["symbol"], "chrom": "chr" + reg["chrom"], "start": reg["start"],
                     "stop": reg["stop"], "reference_repeat_unit": unit, "unit_len": len(unit),
                     "ref_units": (reg["stop"] - reg["start"]) / len(unit),
                     "diseases": "; ".join(d["name"] for d in x["associated_diseases"]), "n_alleles": sum(counts.values())})
        time.sleep(args.sleep)
    pd.DataFrame(rows).to_csv(f"{args.out_prefix}.loci.tsv", sep="\t", index=False)
    pd.DataFrame(sizes).to_csv(f"{args.out_prefix}.allele_sizes.tsv", sep="\t", index=False)
    with open(f"{args.out_prefix}.raw.json", "w") as f:
        json.dump(raw, f)
    print(f"{len(rows)} loci; {len(sizes)} (locus, repunit, size) rows; alleles per locus median "
          f"{pd.DataFrame(rows)['n_alleles'].median():,.0f}")


if __name__ == "__main__":
    main()
