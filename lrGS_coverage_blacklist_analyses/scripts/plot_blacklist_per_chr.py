import os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np

SCR = "/private/tmp/claude-502/-Users-xzhao-Downloads-github/ce274515-e939-4f8a-9d3d-67c167dd0e2b/scratchpad"
DIR = "/Users/xzhao/Downloads/github/lr-annotation/lrGS_coverage_blacklist_analyses"
OUTDIR = f"{DIR}/blacklist_region.per_chr"
os.makedirs(OUTDIR, exist_ok=True)

CHROMS = [f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY"]
THRESHOLDS = [90, 70, 50]

# load downsampled coverage curve
cov = {c: {"pos": [], "low": [], "elig": []} for c in CHROMS}
with open(f"{SCR}/coverage_downsampled.tsv") as f:
    for line in f:
        chrom, start, low, elig = line.rstrip("\n").split("\t")
        if chrom not in cov:
            continue
        cov[chrom]["pos"].append(int(start))
        cov[chrom]["low"].append(int(low))
        cov[chrom]["elig"].append(int(elig))

# load blacklist regions (full resolution) for each threshold
blacklist = {p: {c: [] for c in CHROMS} for p in THRESHOLDS}
for p in THRESHOLDS:
    with open(f"{DIR}/hgsvc_hprc.blacklist.p{p}.raw.bed") as f:
        for line in f:
            chrom, start, end = line.rstrip("\n").split("\t")
            if chrom in blacklist[p]:
                blacklist[p][chrom].append((int(start), int(end)))

blue_patch = mpatches.Patch(color="steelblue", alpha=0.8, label="Low coverage sample count")
red_patch = mpatches.Patch(color="red", alpha=0.35, label="Blacklist region")

for chrom in CHROMS:
    pos = np.array(cov[chrom]["pos"]) / 1e6
    low = np.array(cov[chrom]["low"])
    elig = np.array(cov[chrom]["elig"])
    if len(pos) == 0:
        continue
    eligible_mode = int(np.median(elig))

    fig, axes = plt.subplots(3, 1, figsize=(20, 15), sharex=True)

    for ax, p in zip(axes, THRESHOLDS):
        threshold = (p / 100) * eligible_mode
        ax.fill_between(pos, low, step="post", color="steelblue", alpha=0.8, linewidth=0.5)
        ax.axhline(threshold, color="black", linestyle=":", linewidth=1)

        for start, end in blacklist[p][chrom]:
            ax.axvspan(start / 1e6, end / 1e6, color="red", alpha=0.35, linewidth=0)

        n_regions = len(blacklist[p][chrom])
        total_mb = sum(e - s for s, e in blacklist[p][chrom]) / 1e6
        ax.set_ylabel("Low Coverage\nSample Count")
        ax.set_title(f"blacklist >= {p}% of eligible samples (n={eligible_mode}) -- "
                     f"{n_regions} regions, {total_mb:.2f} Mb", fontsize=11)
        ax.legend(handles=[blue_patch, red_patch], loc="upper right", fontsize=9)

    axes[-1].set_xlabel("Position (Mb)")
    axes[0].set_xlim(pos.min(), pos.max())
    fig.suptitle(chrom, fontsize=14)
    fig.tight_layout()
    fig.savefig(f"{OUTDIR}/{chrom}.blacklist.png", dpi=120)
    plt.close(fig)
    print(f"{chrom}: saved 3-panel plot")

print("all done")
