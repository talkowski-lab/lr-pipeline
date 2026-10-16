#!/usr/bin/env python3
"""Pull per-sample genotypes for a set of TRVs from per-contig VCFs (local or gs://; for gs:// export GCS_OAUTH_TOKEN first).

Inputs:
  --trv-bed   TRV records with columns #chrom start end ID (e.g. split_sr_matched_strs.py .TRV.bed); optional gnomAD_id column
  --vcf-map   two columns: contig, VCF path (indexed)
  --loci      optional fetch_gnomad_str.py .loci.tsv (id, reference_repeat_unit, ref_units) to add repeat units on the gnomAD
              short-read scale: round(ref_units) + (allele length - REF length) / unit length
Outputs (prefix --out-prefix):
  .vcf.gz (+ .tbi)   the TRV records with all samples
  .genotypes.tsv     one row per sample and TRV: gnomAD_id, ID, sample, GT, per haplotype allele index / allele length (bp,
                     padding base removed) / SR-scale units, and the FORMAT fields AL, ALLR, SD, MC, AP, PS
Usage:
  export GCS_OAUTH_TOKEN=$(gcloud auth print-access-token)
  extract_trv_genotypes.py --trv-bed STR_SRmatch.hgsvc_hprc.TRV.bed --vcf-map contigs.tsv --loci STR_SRmatch.sr.loci.tsv \\
      --out-prefix STR_SRmatch.hgsvc_hprc.genotypes
"""
import argparse
import os
import subprocess
import tempfile

import pandas as pd

FMT = ["GT", "AL", "ALLR", "SD", "MC", "AP", "PS"]


def run(cmd, **kw):
    return subprocess.run(cmd, check=True, capture_output=True, text=True, **kw).stdout


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--trv-bed", required=True)
    p.add_argument("--vcf-map", required=True)
    p.add_argument("--loci")
    p.add_argument("--out-prefix", required=True)
    args = p.parse_args()

    bed = pd.read_csv(args.trv_bed, sep="\t", dtype=str)
    bed["pos"] = bed["ID"].str.split("-").str[1].astype(int)
    vcfs = dict(line.split()[:2] for line in open(args.vcf_map) if line.strip())
    tmp = tempfile.mkdtemp()
    parts = []
    for chrom, g in bed.groupby("#chrom", sort=False):
        ids = os.path.join(tmp, f"{chrom}.ids")
        reg = os.path.join(tmp, f"{chrom}.regions")
        with open(ids, "w") as f:
            f.write("\n".join(g["ID"]) + "\n")
        with open(reg, "w") as f:
            for r in g.itertuples():
                f.write(f"{chrom}\t{min(r.pos, int(r.start) + 1)}\t{int(r.end)}\n")
        out = os.path.join(tmp, f"{chrom}.vcf.gz")
        run(["bcftools", "view", "-R", reg, "-i", f"ID=@{ids}", "-Oz", "-o", out, vcfs[chrom]])
        parts.append(out)
    vcf = f"{args.out_prefix}.vcf.gz"
    cat = os.path.join(tmp, "concat.vcf.gz")
    run(["bcftools", "concat", "-Oz", "-o", cat] + parts)
    run(["bcftools", "sort", "-T", tmp, "-Oz", "-o", vcf, cat])
    run(["tabix", "-f", "-p", "vcf", vcf])
    n = int(run(["bcftools", "view", "-H", vcf]).count("\n"))
    missing = set(bed["ID"]) - set(run(["bcftools", "query", "-f", "%ID\n", vcf]).split())
    print(f"{n} records written to {vcf}; missing: {sorted(missing) if missing else 'none'}")

    q = run(["bcftools", "query", "-f", "%ID\t%REF\t%ALT[\t%SAMPLE=" + ":".join(f"%{x}" for x in FMT) + "]\n", vcf])
    units = {}
    if args.loci:
        loci = pd.read_csv(args.loci, sep="\t").set_index("id")
        units = {gid: (round(loci.at[gid, "ref_units"]), len(loci.at[gid, "reference_repeat_unit"])) for gid in loci.index}
    gid_of = dict(zip(bed["ID"], bed["gnomAD_id"])) if "gnomAD_id" in bed else {}
    rows = []
    for line in q.splitlines():
        f = line.split("\t")
        trv, alleles = f[0], [f[1]] + f[2].split(",")
        lens = [len(a) - 1 for a in alleles]
        gid = gid_of.get(trv, ".")
        for s in f[3:]:
            sample, vals = s.split("=", 1)
            d = dict(zip(FMT, vals.split(":")))
            row = {"gnomAD_id": gid, "ID": trv, "sample": sample, "GT": d["GT"]}
            gt = d["GT"].replace("|", "/").split("/")
            for h in (0, 1):
                a = gt[h] if h < len(gt) else "."
                row[f"hap{h + 1}_allele"] = a
                row[f"hap{h + 1}_len_bp"] = lens[int(a)] if a != "." else ""
                if gid in units and a != ".":
                    base, k = units[gid]
                    row[f"hap{h + 1}_SR_units"] = round(base + (lens[int(a)] - lens[0]) / k, 2)
                elif gid in units:
                    row[f"hap{h + 1}_SR_units"] = ""
            row.update({x: d[x] for x in FMT[1:]})
            rows.append(row)
    t = pd.DataFrame(rows)
    t.to_csv(f"{args.out_prefix}.genotypes.tsv", sep="\t", index=False)
    print(f"{len(t):,} sample x TRV rows; {t['sample'].nunique()} samples; {t['ID'].nunique()} TRVs; "
          f"missing GT {(t['GT'].str.contains(r'\.')).sum():,}")


if __name__ == "__main__":
    main()
