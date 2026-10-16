#!/usr/bin/env python3
"""Compare gnomAD short-read STR allele-size distributions with long-read TRV alleles at the same loci.

Inputs:
  --loci, --allele-sizes   fetch_gnomad_str.py outputs (locus coordinates, reference repeat unit, allele counts per
                           repeat-unit count summed over cohorts)
  --trv                    merged TRV table (merge_trv_tables.sh output; tabix-indexed; columns #chrom start end ID FILTER
                           SOURCE motifs AC AN alt_size_diff)
Matching: TRV records overlapping the gnomAD region; the one with the largest overlap is used (ties: shorter span).
Long-read repeat units per allele = round(gnomAD reference units) + (ALT size - REF size) / unit length; REF alleles weighted
by AN - sum(AC), ALT alleles by AC; units rounded to the nearest integer. match_flag marks approximate comparisons:
  wider_TRV    TRV span > 1.5 x the gnomAD region (compound locus: size changes in other components are counted as units)
  partial_TRV  TRV covers < 80% of the gnomAD region
  ok           otherwise
Per locus: allele counts, median / p99 / max units in each dataset, KS distance (max |CDF_SR - CDF_LR|), and the fraction
of long-read alleles above the short-read maximum.
Outputs (prefix --out-prefix):
  .pdf   one panel per locus (allele frequency vs repeat units, log y; dotted line = reference units), then summary pages
  .tsv   per-locus statistics and the matched TRV
  .long_read_alleles.tsv   per locus and repeat-unit count: long-read allele counts

Usage:
  compare_gnomad_str_to_lr.py --loci gnomad_r4_STR.loci.tsv --allele-sizes gnomad_r4_STR.allele_sizes.tsv \
      --trv gnomAD_LR.hprc_hgsvc.TRV.bed.gz --label hprc_hgsvc --out-prefix out/gnomad_STR_vs_LR
"""
import argparse
import math
import subprocess

import matplotlib
import numpy as np
import pandas as pd

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.backends.backend_pdf import PdfPages  # noqa: E402

plt.rcParams["font.family"] = "Arial"
plt.rcParams["pdf.fonttype"] = 42

INK = "#1d1f24"
INK2 = "#5b5e66"
GRID = "#e6e4df"
SR = "#5b5e66"
LR = "#2a78d6"


def trv_header(path):
    out = subprocess.run(f"gzip -cd '{path}' | sed -n '1p'", shell=True, check=True, capture_output=True, text=True).stdout
    return out.rstrip("\n").split("\t")


def overlapping_trvs(path, header, chrom, start, stop):
    out = subprocess.run(["tabix", path, f"{chrom}:{start + 1}-{stop}"], check=True, capture_output=True, text=True).stdout
    rows = [line.split("\t") for line in out.strip().split("\n") if line]
    return pd.DataFrame(rows, columns=header) if rows else pd.DataFrame(columns=header)


def summarize(units, counts):
    units, counts = np.asarray(units, float), np.asarray(counts, float)
    o = np.argsort(units)
    units, counts = units[o], counts[o]
    cdf = np.cumsum(counts) / counts.sum()
    med = units[np.searchsorted(cdf, 0.5)]
    p99 = units[min(np.searchsorted(cdf, 0.99), len(units) - 1)]
    return med, p99, units[counts > 0].max()


def ks(a_units, a_counts, b_units, b_counts):
    grid = np.union1d(a_units, b_units)

    def cdf(u, c):
        s = pd.Series(c, index=u).groupby(level=0).sum().reindex(grid, fill_value=0)
        return np.cumsum(s.to_numpy()) / s.sum()
    return float(np.max(np.abs(cdf(a_units, a_counts) - cdf(b_units, b_counts))))


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--loci", required=True)
    p.add_argument("--allele-sizes", required=True)
    p.add_argument("--trv", required=True)
    p.add_argument("--label", default="long-read")
    p.add_argument("--out-prefix", required=True)
    args = p.parse_args()

    loci = pd.read_csv(args.loci, sep="\t")
    sizes = pd.read_csv(args.allele_sizes, sep="\t")
    header = trv_header(args.trv)
    stats, lr_rows, panels = [], [], []
    for r in loci.itertuples():
        sr = sizes[sizes["id"] == r.id].groupby("repunit_count")["allele_count"].sum()
        row = {"id": r.id, "gene": r.gene, "region": f"{r.chrom}:{r.start + 1}-{r.stop}", "repeat_unit": r.reference_repeat_unit,
               "ref_units": r.ref_units, "SR_alleles": int(sr.sum())}
        if len(sr):
            row["SR_median"], row["SR_p99"], row["SR_max"] = summarize(sr.index, sr.values)
        cand = overlapping_trvs(args.trv, header, r.chrom, r.start, r.stop)
        if cand.empty:
            row["LR_TRV"] = "none"
            row["match_flag"] = "no_TRV"
            stats.append(row)
            panels.append((r, sr, None, row))
            continue
        cand["start"] = cand["start"].astype(int)
        cand["end"] = cand["end"].astype(int)
        cand["overlap"] = np.minimum(cand["end"], r.stop) - np.maximum(cand["start"], r.start)
        cand["span"] = cand["end"] - cand["start"]
        best = cand.sort_values(["overlap", "span"], ascending=[False, True]).iloc[0]
        ac = np.array([int(x) for x in best["AC"].split(",")])
        diff = np.array([int(x) for x in best["alt_size_diff"].split(",")])
        an = int(best["AN"])
        anchor = round(r.ref_units)
        units = np.round(np.concatenate([[anchor], anchor + diff / r.unit_len])).astype(int)
        g_span = r.stop - r.start
        flag = "wider_TRV" if best["span"] > 1.5 * g_span else ("partial_TRV" if best["overlap"] < 0.8 * g_span else "ok")
        counts = np.concatenate([[an - ac.sum()], ac])
        lr = pd.Series(counts, index=units).groupby(level=0).sum()
        lr = lr[lr > 0]
        for u, c in lr.items():
            lr_rows.append({"id": r.id, "repunit_count": int(u), "allele_count": int(c)})
        row.update({"match_flag": flag, "TRV_span_over_gnomAD_span": round(best["span"] / g_span, 2),
                    "LR_TRV": best["ID"], "LR_FILTER": best["FILTER"], "LR_SOURCE": best["SOURCE"], "LR_motifs": best["motifs"],
                    "LR_TRV_span": f"{best['start']}-{best['end']}", "n_TRV_overlapping": len(cand), "LR_alleles": int(lr.sum()),
                    "LR_distinct_ALT": int((ac > 0).sum())})
        row["LR_median"], row["LR_p99"], row["LR_max"] = summarize(lr.index, lr.values)
        if len(sr):
            row["KS_distance"] = ks(sr.index.to_numpy(), sr.values, lr.index.to_numpy(), lr.values)
            row["LR_frac_above_SR_max"] = float(lr[lr.index > row["SR_max"]].sum() / lr.sum())
        stats.append(row)
        panels.append((r, sr, lr, row))

    t = pd.DataFrame(stats)
    t.to_csv(f"{args.out_prefix}.tsv", sep="\t", index=False, float_format="%.4g")
    pd.DataFrame(lr_rows).to_csv(f"{args.out_prefix}.long_read_alleles.tsv", sep="\t", index=False)

    with PdfPages(f"{args.out_prefix}.pdf") as pdf:
        per_page, ncol = 16, 4
        for page in range(math.ceil(len(panels) / per_page)):
            chunk = panels[page * per_page:(page + 1) * per_page]
            nrow = math.ceil(len(chunk) / ncol)
            fig, axes = plt.subplots(nrow, ncol, figsize=(16, 3.0 * nrow), squeeze=False)
            for k, (r, sr, lr, row) in enumerate(chunk):
                ax = axes[k // ncol][k % ncol]
                for dist, color, lab in [(sr, SR, f"gnomAD SR (n={int(sr.sum()):,})"),
                                         (lr, LR, f"{args.label} LR (n={int(lr.sum()):,})" if lr is not None else "")]:
                    if dist is None or not len(dist):
                        continue
                    grid = np.arange(int(dist.index.min()), int(dist.index.max()) + 1)
                    y = (dist.reindex(grid, fill_value=0) / dist.sum()).to_numpy(dtype=float)
                    y[y == 0] = np.nan  # gaps break the line; isolated alleles show as markers
                    ax.plot(grid, y, color=color, linewidth=1.3, marker="o", markersize=2.2, label=lab)
                ax.axvline(r.ref_units, color=INK2, linestyle=":", linewidth=0.8)
                ax.set_yscale("log")
                ks_txt = f"; KS {row['KS_distance']:.2f}" if "KS_distance" in row and not pd.isna(row.get("KS_distance")) else ""
                flag = row.get("match_flag", "")
                note = " (no TRV in long-read data)" if lr is None else ("" if flag == "ok" else f" [{flag}]")
                ax.set_title(f"{r.id} ({r.reference_repeat_unit}){note}{ks_txt}", fontsize=8.5, color=INK)
                ax.set_xlabel("Repeat units", fontsize=7.5, color=INK)
                if k % ncol == 0:
                    ax.set_ylabel("Allele frequency (log)", fontsize=7.5, color=INK)
                for sp in ["top", "right"]:
                    ax.spines[sp].set_visible(False)
                ax.tick_params(colors=INK2, labelsize=7)
                ax.grid(color=GRID, linewidth=0.5)
                ax.legend(frameon=False, fontsize=6.5, labelcolor=INK, loc="upper right")
            for k in range(len(chunk), nrow * ncol):
                axes[k // ncol][k % ncol].axis("off")
            fig.suptitle(f"gnomAD r4 short-read STR loci vs {args.label} long-read TRVs: allele size distribution "
                         f"(dotted = GRCh38 reference units); page {page + 1}", fontsize=10, color=INK)
            fig.tight_layout()
            pdf.savefig(fig)
            plt.close(fig)
        m = t[t["match_flag"] == "ok"].dropna(subset=["LR_median", "SR_median"])
        fig, axes = plt.subplots(1, 3, figsize=(15, 4.6))
        for ax, (a, b, title) in zip(axes, [("SR_median", "LR_median", "Median repeat units"),
                                            ("SR_p99", "LR_p99", "99th percentile repeat units"),
                                            ("SR_max", "LR_max", "Maximum repeat units")]):
            ax.scatter(m[a], m[b], s=18, color=LR)
            for rr in m.itertuples():
                if abs(getattr(rr, b) - getattr(rr, a)) > 0.25 * max(getattr(rr, a), 1):
                    ax.annotate(rr.id, (getattr(rr, a), getattr(rr, b)), fontsize=6.5, color=INK2, xytext=(3, 2),
                                textcoords="offset points")
            lim = [min(m[a].min(), m[b].min()) * 0.8, max(m[a].max(), m[b].max()) * 1.2]
            ax.plot(lim, lim, color=INK2, linestyle="--", linewidth=0.8)
            ax.set_xscale("log")
            ax.set_yscale("log")
            ax.set_xlabel("gnomAD short-read", fontsize=9, color=INK)
            ax.set_ylabel(f"{args.label} long-read", fontsize=9, color=INK)
            ax.set_title(title, fontsize=10, color=INK)
            for sp in ["top", "right"]:
                ax.spines[sp].set_visible(False)
            ax.grid(color=GRID, linewidth=0.5)
        fig.suptitle(f"Summary across {len(m)} loci with a clean match (match_flag ok; labels: loci differing by > 25%)",
                     fontsize=10, color=INK)
        fig.tight_layout()
        pdf.savefig(fig)
        plt.close(fig)
    cols = ["id", "repeat_unit", "match_flag", "SR_alleles", "LR_alleles", "SR_median", "LR_median", "SR_p99", "LR_p99",
            "SR_max", "LR_max", "KS_distance"]
    print(t[[c for c in cols if c in t]].round(2).to_string(index=False))
    print(f"\nloci: {len(t)}; match_flag: {t['match_flag'].value_counts().to_dict()}")


if __name__ == "__main__":
    main()
