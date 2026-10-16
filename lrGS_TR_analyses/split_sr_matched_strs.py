#!/usr/bin/env python3
"""Write the gnomAD short-read STR loci and their matched long-read TRVs into parallel per-dataset files.

Inputs: fetch_gnomad_str.py outputs (--loci, --allele-sizes) and compare_gnomad_str_to_lr.py outputs (--comparison .tsv,
--long-read-alleles .long_read_alleles.tsv), plus the TRV tables used in the comparison (--trv LABEL=NAME=path; LABEL = dataset
label in the comparison, NAME = file name tag; tables tabix-indexed).
Outputs in --out-dir, all named <prefix>.<tag>.* so the datasets sort next to each other:
  <prefix>.sr.loci.tsv              gnomAD loci (id gene chrom start stop repeat unit ref_units diseases n_alleles)
  <prefix>.sr.allele_dist.tsv       id gene repunit_count allele_count (summed over cohorts and repunits)
  <prefix>.NAME.loci.tsv            per gnomAD locus: matched TRV ID, match_flag, span ratio, alleles, median / p99 / max, KS
  <prefix>.NAME.allele_dist.tsv     id gene repunit_count allele_count (long-read units, same scale as SR)
  <prefix>.NAME.TRV.bed             full TRV records of the matched TRVs + gnomAD_id, match_flag (sorted by position)

Usage:
  split_sr_matched_strs.py --loci gnomad_r4_STR.loci.tsv --allele-sizes gnomad_r4_STR.allele_sizes.tsv \\
      --comparison cmp.tsv --long-read-alleles cmp.long_read_alleles.tsv \\
      --trv hprc_hgsvc=hgsvc_hprc=gnomAD_LR.hprc_hgsvc.TRV.bed.gz AoU_I=aou_i=gnomAD_LR.AoU_I.TRV.sites.bed.gz \\
      --prefix STR_SRmatch --out-dir SR_matched
"""
import argparse
import os
import subprocess

import pandas as pd


def trv_records(path, ids):
    header = subprocess.run(f"gzip -cd {path} | sed -n 1p", shell=True, capture_output=True, text=True).stdout
    cols = header.lstrip("#").rstrip("\n").split("\t")
    rows = []
    for trv_id, gid, flag, region in ids:
        out = subprocess.run(["tabix", path, region], capture_output=True, text=True, check=True).stdout
        for line in out.splitlines():
            f = line.split("\t")
            if f[cols.index("ID")] == trv_id:
                rows.append(f + [gid, flag])
                break
        else:
            raise SystemExit(f"{trv_id} ({gid}) not found in {path} at {region}")
    t = pd.DataFrame(rows, columns=cols + ["gnomAD_id", "match_flag"])
    t["_s"] = t[cols[1]].astype(int)
    t = t.sort_values([cols[0], "_s"]).drop(columns="_s")
    return t.rename(columns={cols[0]: "#" + cols[0]})


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--loci", required=True)
    p.add_argument("--allele-sizes", required=True)
    p.add_argument("--comparison", required=True)
    p.add_argument("--long-read-alleles", required=True)
    p.add_argument("--trv", nargs="+", required=True, help="LABEL=NAME=path")
    p.add_argument("--prefix", required=True)
    p.add_argument("--out-dir", required=True)
    args = p.parse_args()

    os.makedirs(args.out_dir, exist_ok=True)
    out = os.path.join(args.out_dir, args.prefix)
    loci = pd.read_csv(args.loci, sep="\t")
    gene = loci.set_index("id")["gene"]
    loci.to_csv(f"{out}.sr.loci.tsv", sep="\t", index=False)
    sr = pd.read_csv(args.allele_sizes, sep="\t").groupby(["id", "repunit_count"], as_index=False)["allele_count"].sum()
    sr.insert(1, "gene", sr["id"].map(gene))
    sr.to_csv(f"{out}.sr.allele_dist.tsv", sep="\t", index=False)
    print(f"sr: {len(loci)} loci, {int(sr['allele_count'].sum()):,} alleles")

    cmp_ = pd.read_csv(args.comparison, sep="\t")
    lra = pd.read_csv(args.long_read_alleles, sep="\t")
    for spec in args.trv:
        label, name, path = spec.split("=", 2)
        lcols = [c for c in cmp_.columns if c.startswith(f"{label}_")]
        t = cmp_[["id", "gene", "region", "repeat_unit", "ref_units"] + lcols].rename(
            columns={c: c[len(label) + 1:] for c in lcols})
        t.to_csv(f"{out}.{name}.loci.tsv", sep="\t", index=False)
        d = lra[lra["dataset"] == label].drop(columns="dataset")
        d.insert(1, "gene", d["id"].map(gene))
        d.to_csv(f"{out}.{name}.allele_dist.tsv", sep="\t", index=False)
        m = t[t["TRV"].notna() & (t["match_flag"] != "no_TRV")]
        bed = trv_records(path, list(zip(m["TRV"], m["id"], m["match_flag"], m["region"])))
        bed.to_csv(f"{out}.{name}.TRV.bed", sep="\t", index=False)
        print(f"{name}: {len(m)} matched loci ({(m['match_flag'] == 'ok').sum()} ok), {bed['ID'].nunique()} distinct TRVs, "
              f"{int(d['allele_count'].sum()):,} alleles")


if __name__ == "__main__":
    main()
