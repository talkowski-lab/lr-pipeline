#!/usr/bin/env python3
"""STR density (count per Mb) across coding, intronic, intergenic non-desert and intergenic desert, overall and by motif.

Inputs:
  --trv         merged TRV table (merge_trv_tables.sh output; columns #chrom FILTER motifs AC genic_context intergenic_category)
  --gtf         GENCODE GTF used for genic_context (CDS / exon / transcript)
  --intergenic  intergenic regions BED with header chr start end size is_desert ... size_stratum
STRs: FILTER PASS, summed AC > 0, shortest motif in [--min-motif, --max-motif] (default 1-6 bp), autosomes (chr1-22; the
intergenic region file is autosomal). Contexts and territory (denominator):
  coding                 merged CDS bases
  intronic               merged transcript spans minus merged exons
  intergenic_non_desert  Tiny / Small / Medium regions minus MANE transcript bases; STRs with that intergenic_category
  intergenic_desert      Desert regions minus MANE transcript bases; STRs with intergenic_category Desert
Intergenic STRs outside the listed regions are not counted.
Outputs: <out-prefix>.pdf (left: STRs per Mb by context; right: STRs per Mb by motif size per context, log y) and .tsv.

Usage:
  str_density_by_context.py --trv gnomAD_LR.hprc_hgsvc.TRV.bed.gz --gtf gencode.gtf.gz \
      --intergenic intergenic_regions_annotated.bed.gz --label hprc_hgsvc --out-prefix out/STR_density
"""
import argparse
import collections
import gzip

import matplotlib
import pandas as pd

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

plt.rcParams["font.family"] = "Arial"
plt.rcParams["pdf.fonttype"] = 42

INK = "#1d1f24"
INK2 = "#5b5e66"
GRID = "#e6e4df"
AUTOSOMES = {f"chr{i}" for i in range(1, 23)}
CONTEXTS = ["coding", "intronic", "intergenic_non_desert", "intergenic_desert"]
LABELS = {"coding": "Coding", "intronic": "Intronic", "intergenic_non_desert": "Intergenic\nnon-desert",
          "intergenic_desert": "Intergenic\ndesert"}
COLORS = {"coding": "#2a78d6", "intronic": "#1baf7a", "intergenic_non_desert": "#eda100", "intergenic_desert": "#8a6d3b"}
NON_DESERT = {"Tiny (<10kb)", "Small (10–100kb)", "Medium (100–460kb)"}
DESERT = {"Desert (≥460kb)"}


def merge(iv):
    out = []
    for s, e in sorted(iv):
        if out and s <= out[-1][1]:
            out[-1][1] = max(out[-1][1], e)
        else:
            out.append([s, e])
    return out


def subtract(a, b):
    """Merged intervals a minus merged intervals b."""
    out, j = [], 0
    for s, e in a:
        cur = s
        while j < len(b) and b[j][1] <= cur:
            j += 1
        k = j
        while k < len(b) and b[k][0] < e:
            if b[k][0] > cur:
                out.append([cur, b[k][0]])
            cur = max(cur, b[k][1])
            k += 1
        if cur < e:
            out.append([cur, e])
    return out


def total(iv):
    return sum(e - s for s, e in iv)


def territories(gtf, regions_path):
    feats = collections.defaultdict(lambda: collections.defaultdict(list))
    with gzip.open(gtf, "rt") as f:
        for line in f:
            if line.startswith("#"):
                continue
            t = line.split("\t", 5)
            if t[0] in AUTOSOMES and t[2] in ("CDS", "exon", "transcript"):
                feats[t[2]][t[0]].append((int(t[3]) - 1, int(t[4])))
    r = pd.read_csv(regions_path, sep="\t", dtype={"chr": str})
    r = r[r["chr"].isin(AUTOSOMES)]
    bp = collections.Counter()
    for chrom in AUTOSOMES:
        cds = merge(feats["CDS"][chrom])
        exons = merge(feats["exon"][chrom])
        tx = merge(feats["transcript"][chrom])
        bp["coding"] += total(cds)
        bp["intronic"] += total(subtract(tx, exons))
        sub = r[r["chr"] == chrom]
        for ctx, cats in [("intergenic_non_desert", NON_DESERT), ("intergenic_desert", DESERT)]:
            sel = sub[sub["size_stratum"].isin(cats)]
            reg = merge(list(zip(sel["start"], sel["end"])))
            bp[ctx] += total(subtract(reg, tx))
    return bp


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--trv", required=True)
    p.add_argument("--gtf", required=True)
    p.add_argument("--intergenic", required=True)
    p.add_argument("--min-motif", type=int, default=1)
    p.add_argument("--max-motif", type=int, default=6)
    p.add_argument("--label", default="")
    p.add_argument("--out-prefix", required=True)
    args = p.parse_args()

    bp = territories(args.gtf, args.intergenic)
    counts = collections.defaultdict(collections.Counter)
    cols = ["#chrom", "FILTER", "motifs", "AC", "genic_context", "intergenic_category"]
    for ch in pd.read_csv(args.trv, sep="\t", usecols=cols, dtype=str, chunksize=500000):
        ch = ch[(ch["FILTER"] == "PASS") & ch["#chrom"].isin(AUTOSOMES)]
        motif = ch["motifs"].map(lambda m: min(len(x) for x in m.split(",")))
        nonref = ch["AC"].map(lambda s: sum(int(x) for x in s.split(",") if x not in (".", "")) > 0)
        ch = ch.assign(motif=motif)[nonref & motif.between(args.min_motif, args.max_motif)]
        ctx = ch["genic_context"].copy()
        inter = ctx == "intergenic"
        ctx[inter] = "."
        ctx[inter & ch["intergenic_category"].isin(NON_DESERT)] = "intergenic_non_desert"
        ctx[inter & ch["intergenic_category"].isin(DESERT)] = "intergenic_desert"
        for (c, m), n in ch.assign(ctx=ctx).groupby(["ctx", "motif"]).size().items():
            if c in CONTEXTS:
                counts[c][m] += n

    motifs = list(range(args.min_motif, args.max_motif + 1))
    rows = []
    for c in CONTEXTS:
        mb = bp[c] / 1e6
        n_all = sum(counts[c].values())
        rows.append({"context": c, "motif": "all", "n_STR": n_all, "territory_Mb": round(mb, 2), "STR_per_Mb": n_all / mb})
        for m in motifs:
            rows.append({"context": c, "motif": m, "n_STR": counts[c][m], "territory_Mb": round(mb, 2),
                         "STR_per_Mb": counts[c][m] / mb})
    t = pd.DataFrame(rows)
    t.to_csv(f"{args.out_prefix}.tsv", sep="\t", index=False, float_format="%.4g")

    fig, axes = plt.subplots(1, 2, figsize=(13, 4.8), gridspec_kw={"width_ratios": [1, 1.4]})
    ax = axes[0]
    a = t[t["motif"] == "all"]
    bars = ax.bar([LABELS[c] for c in a["context"]], a["STR_per_Mb"], color=[COLORS[c] for c in a["context"]],
                  edgecolor="white", linewidth=1)
    for b, r in zip(bars, a.itertuples()):
        ax.text(b.get_x() + b.get_width() / 2, b.get_height(),
                f"{r.STR_per_Mb:,.0f} / Mb\n{r.n_STR:,} STRs\n{r.territory_Mb:,.0f} Mb", ha="center", va="bottom",
                fontsize=7.5, color=INK)
    ax.set_ylim(0, a["STR_per_Mb"].max() * 1.35)
    ax.set_ylabel("Non-ref PASS STRs per Mb", fontsize=9, color=INK)
    ax.set_title("All STR motifs (1–6 bp)", fontsize=10, color=INK)
    ax = axes[1]
    for c in CONTEXTS:
        s = t[(t["context"] == c) & (t["motif"] != "all")]
        ax.plot(s["motif"].astype(int), s["STR_per_Mb"], color=COLORS[c], linewidth=2, marker="o", markersize=5,
                label=LABELS[c].replace("\n", " "))
    ax.set_yscale("log")
    ax.set_xticks(motifs)
    ax.set_xlabel("Repeat motif size (bp)", fontsize=9, color=INK)
    ax.set_ylabel("Non-ref PASS STRs per Mb (log)", fontsize=9, color=INK)
    ax.set_title("By motif size", fontsize=10, color=INK)
    ax.legend(frameon=False, fontsize=8.5, labelcolor=INK)
    for a_ in axes:
        for sp in ["top", "right"]:
            a_.spines[sp].set_visible(False)
        a_.tick_params(colors=INK2, labelsize=8)
        a_.grid(axis="y", color=GRID, linewidth=0.6)
        a_.set_axisbelow(True)
    fig.suptitle(f"{args.label + ': ' if args.label else ''}STR density by genomic context (autosomes; territory = CDS, "
                 f"introns, and intergenic regions minus MANE transcripts)", fontsize=10, color=INK)
    fig.tight_layout()
    fig.savefig(f"{args.out_prefix}.pdf")
    plt.close(fig)
    wide = t.pivot_table(index="context", columns="motif", values="STR_per_Mb", aggfunc="first").reindex(CONTEXTS)
    print(wide.round(1).to_string())
    print(a[["context", "n_STR", "territory_Mb", "STR_per_Mb"]].round(1).to_string(index=False))


if __name__ == "__main__":
    main()
