#!/usr/bin/env python3
"""Allele distributions at one STR locus: gnomAD short-read (ExpansionHunter) repeat-unit counts vs long-read TRV alleles.

Long-read input: one or more TRV tables (--trv LABEL=path; bgzipped + tabix-indexed, header with columns ID, REF, ALT, AC, AN;
e.g. the merged TRV table or the VEP-parsed TRV bed). The record is chosen by --trv-id. ALT alleles are weighted by AC and REF
by AN - sum(AC). Per allele (padding base removed):
  length_units      allele length / motif length (what a length-based comparison uses)
  longest_pure_run  longest uninterrupted run of the motif (any rotation, any phase), in units
  motif_units       number of motif copies in the allele (any rotation, non-overlapping, greedy left to right)
Short-read input: fetch_gnomad_str.py .allele_sizes.tsv (--sr-sizes) and locus id (--sr-id); allele counts summed over repunits.
Outputs (prefix --out-prefix):
  .alleles.tsv  per dataset and distinct allele: AC, length_bp, length_units, longest_pure_run, motif_units, sequence
  .summary.tsv  per dataset and measure: n alleles, median, p1, p99, min, max
  .pdf          AF distribution: SR repunit count (grey) vs each long-read measure (one panel per measure)

Usage:
  str_locus_allele_distribution.py --trv hprc_hgsvc=gnomAD_LR.hprc_hgsvc.TRV.bed.gz AoU_I=gnomAD_LR.AoU_I.vep_parsed.annotated.TRV.bed.gz \
      --trv-id chr16-17470907-TRV-15 --motif GCC --sr-sizes gnomad_r4_STR.allele_sizes.tsv --sr-id XYLT1 --out-prefix XYLT1/XYLT1
"""
import argparse
import re
import subprocess

import matplotlib
import numpy as np
import pandas as pd

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

plt.rcParams["font.family"] = "Arial"
plt.rcParams["pdf.fonttype"] = 42

INK = "#1d1f24"
INK2 = "#5b5e66"
GRID = "#e6e4df"
SR_COLOR = "#8a8984"
MEASURES = {"length_units": "Allele length / motif length", "longest_pure_run": "Longest pure motif run (units)",
            "motif_units": "Motif copies in allele (interruptions allowed)"}


def rotations(motif):
    return {motif[i:] + motif[:i] for i in range(len(motif))}


def longest_pure_run(seq, motif):
    k, best = len(motif), 0
    for r in rotations(motif):
        for m in re.finditer(f"(?:{r})+", seq):
            best = max(best, len(m.group()) // k)
    return best


def motif_units(seq, motif):
    pat = "|".join(sorted(rotations(motif)))
    return len(re.findall(f"(?:{pat})", seq))


def read_record(path, trv_id):
    header = subprocess.run(f"gzip -cd {path} | sed -n 1p", shell=True, capture_output=True, text=True).stdout
    cols = header.lstrip("#").rstrip("\n").split("\t")
    chrom, pos = trv_id.split("-")[:2]
    region = f"{chrom}:{max(1, int(pos) - 5)}-{int(pos) + 5}"
    out = subprocess.run(["tabix", path, region], capture_output=True, text=True, check=True).stdout
    for line in out.splitlines():
        rec = dict(zip(cols, line.split("\t")))
        if rec["ID"] == trv_id:
            return rec
    raise SystemExit(f"{trv_id} not found in {path}")


def long_read_alleles(rec, motif):
    alts, acs = rec["ALT"].split(","), [int(x) if x not in (".", "") else 0 for x in rec["AC"].split(",")]
    seqs = [rec["REF"]] + alts
    counts = [int(rec["AN"]) - sum(acs)] + acs
    rows = []
    for i, (s, c) in enumerate(zip(seqs, counts)):
        if c <= 0:
            continue
        a = s[1:]
        rows.append({"allele": "REF" if i == 0 else f"ALT{i}", "AC": c, "length_bp": len(a),
                     "length_units": round(len(a) / len(motif), 2), "longest_pure_run": longest_pure_run(a, motif),
                     "motif_units": motif_units(a, motif), "sequence": a})
    return pd.DataFrame(rows)


def wquant(v, w, q):
    o = np.argsort(v)
    v, c = np.asarray(v)[o], np.cumsum(np.asarray(w)[o])
    return v[np.searchsorted(c, q * c[-1])]


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--trv", nargs="+", required=True, help="LABEL=path (one or more long-read TRV tables)")
    p.add_argument("--colors", nargs="+", default=["#2a78d6", "#eb6834", "#1baf7a"])
    p.add_argument("--trv-id", required=True)
    p.add_argument("--motif", required=True)
    p.add_argument("--sr-sizes", required=True)
    p.add_argument("--sr-id", required=True)
    p.add_argument("--out-prefix", required=True)
    args = p.parse_args()

    sr = pd.read_csv(args.sr_sizes, sep="\t")
    sr = sr[sr["id"] == args.sr_id].groupby("repunit_count")["allele_count"].sum()
    datasets = {}
    for spec in args.trv:
        label, path = spec.split("=", 1)
        datasets[label] = long_read_alleles(read_record(path, args.trv_id), args.motif)
    allele_rows = pd.concat([d.assign(dataset=k) for k, d in datasets.items()])
    allele_rows = allele_rows.sort_values(["dataset", "AC"], ascending=[True, False])
    allele_rows[["dataset"] + [c for c in allele_rows.columns if c != "dataset"]].to_csv(
        f"{args.out_prefix}.alleles.tsv", sep="\t", index=False)

    summ = [{"dataset": "gnomAD_SR", "measure": "repunit_count", "n": int(sr.sum()), "median": wquant(sr.index, sr, 0.5),
             "p1": wquant(sr.index, sr, 0.01), "p99": wquant(sr.index, sr, 0.99), "min": sr.index.min(),
             "max": sr.index.max()}]
    for k, d in datasets.items():
        for m in MEASURES:
            summ.append({"dataset": k, "measure": m, "n": int(d["AC"].sum()), "median": wquant(d[m], d["AC"], 0.5),
                         "p1": wquant(d[m], d["AC"], 0.01), "p99": wquant(d[m], d["AC"], 0.99), "min": d[m].min(),
                         "max": d[m].max()})
    summ = pd.DataFrame(summ)
    summ.to_csv(f"{args.out_prefix}.summary.tsv", sep="\t", index=False)
    print(summ.to_string(index=False))

    fig, axes = plt.subplots(len(MEASURES), 1, figsize=(11, 3.2 * len(MEASURES)), sharex=True)
    for ax, (m, mlabel) in zip(axes, MEASURES.items()):
        ax.bar(sr.index, sr / sr.sum(), width=0.9, color=SR_COLOR, alpha=0.8,
               label=f"gnomAD SR repunit count (n = {int(sr.sum()):,})")
        for (k, d), col in zip(datasets.items(), args.colors):
            g = d.groupby(np.round(d[m]).astype(int))["AC"].sum()
            ax.plot(g.index, g / g.sum(), color=col, marker="o", markersize=3.5, linewidth=0, label=f"{k} LR (n = {int(g.sum()):,})")
            ax.vlines(g.index, 0, g / g.sum(), color=col, linewidth=1.2)
        ax.set_yscale("log")
        ax.set_ylabel("Allele frequency (log)", fontsize=9, color=INK)
        ax.set_title(f"Long-read measure: {mlabel}", fontsize=9.5, color=INK, loc="left")
        ax.legend(frameon=False, fontsize=8, labelcolor=INK, loc="upper right")
        for sp in ["top", "right"]:
            ax.spines[sp].set_visible(False)
        ax.tick_params(colors=INK2, labelsize=8)
        ax.grid(color=GRID, linewidth=0.6)
        ax.set_axisbelow(True)
    axes[-1].set_xlabel(f"Repeat units ({args.motif})", fontsize=9, color=INK)
    fig.suptitle(f"{args.sr_id} ({args.trv_id}): gnomAD short-read vs long-read allele distributions", fontsize=10.5, color=INK)
    fig.tight_layout()
    fig.savefig(f"{args.out_prefix}.pdf")
    plt.close(fig)


if __name__ == "__main__":
    main()
