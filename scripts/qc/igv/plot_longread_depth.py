#!/usr/bin/env python3
# Derived from scripts/plot_longread_depth.py of yuliamostovoy/variant-interpretation
# https://github.com/yuliamostovoy/variant-interpretation/blob/main/scripts/plot_longread_depth.py
"""
Long-read read-depth plot for CNV (DEL/DUP) loci.

For each variant in the per-family BED it draws a normalized-depth profile across the
region +/- flank, one line per family sample (carriers highlighted, others gray), shades
the called interval, and marks the expected het DEL (0.5x) / DUP (1.5x) levels. One PNG per
variant is written to --outdir, named {family}_{ID}.png so it pairs with the IGV plot of
the same variant.

Depth comes from mosdepth `--by <windows>` output: one {sample}.regions.bed.gz per sample,
columns chrom,start,end,depth (mean depth per window). Normalization is local and per
sample: divide by the median depth over that sample's windows lying in the flanks (outside
the called interval), which needs no genome-wide coverage stats.
"""

import argparse
import glob
import gzip
import os
import sys
from statistics import median

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt


def read_ped(ped_file):
    """Return dict sample -> {'family','affected'} using the first 6 standard columns."""
    info = {}
    with open(ped_file) as fh:
        for line in fh:
            if line.startswith(("#", "FamilyID", "family_id")):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 6:
                continue
            fam, iid, _mo, _fa, _sex, aff = f[:6]
            info[iid] = {"family": fam, "affected": aff}
    return info


def load_depth(path):
    """Load a mosdepth regions.bed.gz -> list of (chrom, start, end, depth)."""
    rows = []
    with gzip.open(path, "rt") as fh:
        for line in fh:
            c, s, e, d = line.rstrip("\n").split("\t")[:4]
            rows.append((c, int(s), int(e), float(d)))
    return rows


def sample_from_filename(path):
    base = os.path.basename(path)
    return base.split(".regions.bed.gz")[0]


# colorblind-safe (Okabe-Ito) palette for annotation tracks
ANNOTATION_PALETTE = ["#E69F00", "#56B4E9", "#009E73", "#D55E00", "#CC79A7", "#0072B2", "#F0E442", "#999999"]


def load_bed_intervals(path):
    """Load a (optionally gzipped) BED -> {chrom: [(start, end), ...]}."""
    by_chrom = {}
    opener = gzip.open if path.endswith(".gz") else open
    with opener(path, "rt") as fh:
        for line in fh:
            if not line.strip() or line.startswith(("#", "track", "browser")):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 3:
                continue
            try:
                by_chrom.setdefault(f[0], []).append((int(f[1]), int(f[2])))
            except ValueError:
                continue
    return by_chrom


def norm_factor(windows, chrom, cnv_start, cnv_end):
    """Median depth over flank windows (outside the called interval) on this chrom."""
    flank = [d for (c, s, e, d) in windows if c == chrom and (e <= cnv_start or s >= cnv_end)]
    pool = flank if len(flank) >= 3 else [d for (c, s, e, d) in windows if c == chrom]
    pool = [d for d in pool if d > 0]
    return median(pool) if pool else 1.0


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--bed", required=True, help="per-family variant BED: chrom,start,end,ID,svtype,samples")
    ap.add_argument("--ped", required=True)
    ap.add_argument("--family", required=True)
    ap.add_argument("--flank", type=int, default=5000, help="minimum bp of flank shown on each side of the interval (floor)")
    ap.add_argument(
        "--flank-frac",
        type=float,
        default=0.1,
        help="flank as a fraction of event length; per-variant flank = max(--flank, frac*SVLEN)",
    )
    ap.add_argument("--depth-dir", default=".", help="directory containing {sample}.regions.bed.gz files")
    ap.add_argument("--outdir", default="rd_plots")
    ap.add_argument(
        "--annotation-beds", nargs="*", default=[], help="optional BED files of regions to highlight (e.g. N-gaps, segdups)"
    )
    ap.add_argument("--annotation-names", nargs="*", default=[], help="labels for --annotation-beds, in the same order")
    args = ap.parse_args()

    os.makedirs(args.outdir, exist_ok=True)
    ped = read_ped(args.ped)
    affected = {s for s, info in ped.items() if str(info.get("affected")) == "2"}

    annotations = []
    for i, bed_path in enumerate(args.annotation_beds):
        name = args.annotation_names[i] if i < len(args.annotation_names) else os.path.basename(bed_path)
        annotations.append(
            {"name": name, "color": ANNOTATION_PALETTE[i % len(ANNOTATION_PALETTE)], "by_chrom": load_bed_intervals(bed_path)}
        )

    depth = {}
    for path in glob.glob(os.path.join(args.depth_dir, "*.regions.bed.gz")):
        depth[sample_from_filename(path)] = load_depth(path)
    if not depth:
        sys.exit(f"ERROR: no *.regions.bed.gz found in {args.depth_dir}")

    n_plots = 0
    with open(args.bed) as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            chrom, start, end, vid, svtype = f[0], int(f[1]), int(f[2]), f[3], f[4].upper()
            carriers = set(f[5].split(",")) if len(f) > 5 and f[5] else set()
            if svtype not in ("DEL", "DUP"):
                continue

            flank = max(args.flank, int(args.flank_frac * (end - start)))
            lo, hi = max(0, start - flank), end + flank
            # 16in x 120dpi = 1920px wide, matching the IGV snapshot width so the combined
            # figure needs no rescaling
            fig, ax = plt.subplots(figsize=(16, 4))

            for sample, windows in sorted(depth.items()):
                factor = norm_factor(windows, chrom, start, end)
                xs, ys = [], []
                for c, s, e, d in windows:
                    if c == chrom and e > lo and s < hi:
                        xs.append((s + e) / 2.0)
                        ys.append(d / factor if factor else 0.0)
                if not xs:
                    continue
                pts = sorted(zip(xs, ys))
                xs, ys = [p[0] for p in pts], [p[1] for p in pts]
                is_carrier = sample in carriers
                tags = (["carrier"] if is_carrier else []) + (["affected"] if sample in affected else [])
                label = sample + (" (" + ", ".join(tags) + ")" if tags else "")
                ax.plot(
                    xs,
                    ys,
                    color="crimson" if is_carrier else "0.7",
                    lw=2.0 if is_carrier else 1.0,
                    zorder=3 if is_carrier else 1,
                    label=label,
                )

            ax.axvspan(start, end, color="0.6", alpha=0.35, zorder=0)

            # shade any annotation intervals overlapping the visible window
            for track in annotations:
                first = True
                for s, e in track["by_chrom"].get(chrom, []):
                    if e > lo and s < hi:
                        ax.axvspan(
                            max(s, lo),
                            min(e, hi),
                            color=track["color"],
                            alpha=0.18,
                            zorder=0,
                            label=track["name"] if first else "_nolegend_",
                        )
                        first = False

            ax.axhline(1.0, color="black", ls="--", lw=0.8, zorder=2)
            ax.axhline(0.5 if svtype == "DEL" else 1.5, color="green", ls=":", lw=0.8, zorder=2)
            ax.set_ylim(0, 3)
            ax.set_xlabel(f"{chrom} position")
            ax.set_ylabel("normalized depth")
            # SV info is rendered as a header on the combined figure (concat_igv_depth.py)
            ax.legend(fontsize=6, ncol=2, loc="upper right")
            fig.tight_layout()
            fig.savefig(os.path.join(args.outdir, f"{args.family}_{vid}.png"), dpi=120)
            plt.close(fig)
            n_plots += 1

    sys.stderr.write(f"Wrote {n_plots} depth plots to {args.outdir}\n")


if __name__ == "__main__":
    main()
