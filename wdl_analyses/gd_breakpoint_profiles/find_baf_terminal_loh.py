#!/usr/bin/env python3
"""Rank samples by evidence for loss of heterozygosity (LOH, e.g. a heterozygous deletion) over the right (or left)
end of the region in a wide BAF table (CHROM POS ID REF ALT + one column per sample).

Clean sites: biallelic SNVs that behave like ordinary polymorphisms, not paralog (SD) differences:
  het (--het-lo < BAF < --het-hi) in >= --min-het samples, hom-alt (BAF >= 0.95) in >= --min-homalt samples,
  median BAF of het-looking values within 0.5 +/- --max-het-dev, and >= --min-clean-frac of het-looking values
  inside the het band.
Per sample (on clean sites, sorted by position): the het-free segment runs from just after the sample's last
(--side right) or before its first (--side left) clean het to the region end. Expected hets in that segment =
sum of the cohort het frequency of its sites, scaled by the sample's own het rate over the rest of the region
(observed / expected), giving lambda; P(no het in the segment) is about exp(-lambda), so high lambda = strong LOH.

Outputs:
  <out-prefix>.ranked.tsv     one row per sample, sorted by lambda (sample, obs_outside, rate, loh_start/loh_end,
                              clean_sites_in_loh, exp_in_loh, lambda, homalt_in_loh, bp_het_bound, bp_est, bp_ci95_bound)
Breakpoint estimate (the inner edge of the LOH segment, i.e. where the deletion would begin):
  bp_het_bound   position of the sample's last (right side) / first (left side) clean het; the breakpoint cannot be on the
                 far side of it
  bp_est         moving from bp_het_bound into the LOH segment, the clean site where the expected het count
                 (sample rate x cohort het frequency, summed) reaches ln 2: the median breakpoint if any position
                 between het-bearing sites is equally likely a priori
  bp_ci95_bound  where that expected count reaches ln 20 (about 3): with 95% probability the breakpoint lies between
                 bp_het_bound and this site (it would otherwise have left at least one het visible)
  <out-prefix>.clean_sites.baf.tsv.gz  the input table restricted to clean sites (for plotting)

Usage:
  find_baf_terminal_loh.py --baf aou_II/chr16_bp4.GD_carriers.baf.tsv.gz --out-prefix aou_II/chr16_bp4.right_loh \\
      [--side right|left] [--top 20] [--format markdown|tsv]
Works on any BAF table of this layout; sample names are taken from the header (columns 6 onward).
Printed to stdout: the top --top samples as a table (rank, sample, het-free start in Mb, expected hets missing,
lambda, hom-alt sites in the het-free stretch, last het, breakpoint estimate, 95% interval), as markdown (default) or TSV.
"""
import argparse
import sys
import warnings

import numpy as np
import pandas as pd


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--baf", required=True)
    p.add_argument("--out-prefix", required=True)
    p.add_argument("--side", choices=["right", "left"], default="right")
    p.add_argument("--het-lo", type=float, default=0.2)
    p.add_argument("--het-hi", type=float, default=0.8)
    p.add_argument("--min-het", type=int, default=10)
    p.add_argument("--min-homalt", type=int, default=3)
    p.add_argument("--max-het-dev", type=float, default=0.1)
    p.add_argument("--min-clean-frac", type=float, default=0.8)
    p.add_argument("--top", type=int, default=20, help="rows printed to stdout")
    p.add_argument("--format", choices=["markdown", "tsv"], default="markdown", help="printed table format")
    if len(sys.argv) == 1:
        p.print_help()
        sys.exit(1)
    a = p.parse_args()

    d = pd.read_csv(a.baf, sep="\t")
    samples = list(d.columns[5:])
    x = d[samples].to_numpy(dtype=float)
    snv = (d["REF"].str.len().eq(1) & d["ALT"].str.len().eq(1)).to_numpy()
    anyhet = (x > 0.05) & (x < 0.95)
    het = (x > a.het_lo) & (x < a.het_hi)
    homalt = x >= 0.95
    n_any = anyhet.sum(1)
    with np.errstate(invalid="ignore"), warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        med = np.nanmedian(np.where(anyhet, x, np.nan), axis=1)
    clean = (snv & (het.sum(1) >= a.min_het) & (homalt.sum(1) >= a.min_homalt)
             & (np.abs(med - 0.5) < a.max_het_dev) & (het.sum(1) / np.maximum(n_any, 1) >= a.min_clean_frac))
    print(f"{clean.sum()} clean sites of {len(d)}", file=sys.stderr)
    dc = d[clean].sort_values("POS")
    dc.to_csv(a.out_prefix + ".clean_sites.baf.tsv.gz", sep="\t", index=False, compression="gzip")

    pos = dc["POS"].to_numpy()
    xc = dc[samples].to_numpy(dtype=float)
    hc = (xc > a.het_lo) & (xc < a.het_hi)
    ac = xc >= 0.95
    if a.side == "left":
        pos, hc, ac = pos[::-1], hc[::-1], ac[::-1]
    f = hc.mean(1)
    rows = []
    for j, s in enumerate(samples):
        nz = np.where(hc[:, j])[0]
        k = nz.max() + 1 if len(nz) else 0
        exp_out, obs_out = f[:k].sum(), hc[:k, j].sum()
        rate = obs_out / exp_out if exp_out > 0 else np.nan
        exp_in = f[k:].sum()
        seg = pos[k:]
        cum = rate * np.cumsum(f[k:]) if np.isfinite(rate) else np.full(len(seg), np.nan)

        def reach(level):
            i = np.where(cum >= level)[0]
            return int(seg[i[0]]) if len(i) else (int(seg[-1]) if len(seg) else np.nan)
        bp = dict(bp_het_bound=int(pos[k - 1]) if k > 0 else np.nan, bp_est=reach(np.log(2)), bp_ci95_bound=reach(np.log(20)))
        rows.append(dict(sample=s, obs_outside=int(obs_out), rate=round(rate, 2),
                         loh_start=int(seg.min()) if len(seg) else np.nan, loh_end=int(seg.max()) if len(seg) else np.nan,
                         clean_sites_in_loh=len(seg), exp_in_loh=round(exp_in, 1),
                         **{"lambda": round(rate * exp_in, 1)}, homalt_in_loh=int(ac[k:, j].sum()), **bp))
    r = pd.DataFrame(rows).sort_values("lambda", ascending=False)
    r.to_csv(a.out_prefix + ".ranked.tsv", sep="\t", index=False)

    # Printed table: rank, sample, start of the het-free end segment, expected hets missing there, score, hom-alt sites
    t = r.head(a.top).reset_index(drop=True)
    side = "from" if a.side == "right" else "up to"
    header = ["Rank", "Sample", f"Het-free {side} (Mb)", "Expected hets missing", "Score (lambda)", "Hom-alt in het-free stretch",
              "Last het (Mb)" if a.side == "right" else "First het (Mb)", "Breakpoint est. (Mb)", "95% interval (Mb)"]
    loh_edge = t["loh_start"] if a.side == "right" else t["loh_end"]

    def mb(v):
        return "NA" if pd.isna(v) else f"{v / 1e6:.3f}"

    def ci(i):
        lo, hi = sorted([t.at[i, "bp_het_bound"], t.at[i, "bp_ci95_bound"]])
        return f"{mb(lo)}-{mb(hi)}"
    body = [[str(i + 1), t.at[i, "sample"], mb(loh_edge[i]), f"{t.at[i, 'exp_in_loh']:.1f}",
             f"{t.at[i, 'lambda']:.1f}", str(t.at[i, "homalt_in_loh"]), mb(t.at[i, "bp_het_bound"]), mb(t.at[i, "bp_est"]),
             ci(i)]
            for i in range(len(t))]
    if a.format == "markdown":
        print("| " + " | ".join(header) + " |")
        print("|" + "|".join("---" for _ in header) + "|")
        for row in body:
            print("| " + " | ".join(row) + " |")
    else:
        print("\t".join(header))
        for row in body:
            print("\t".join(row))


if __name__ == "__main__":
    main()
