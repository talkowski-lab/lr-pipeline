#!/usr/bin/env python3
"""Figures and summary numbers for the TRV / STR slide deck.

Inputs
  --overview LABEL=variant_overview.tsv ...     from count_variant_overview.sh
  --loci LABEL=TRV_loci.tsv.gz ...              from build_trv_locus_table.py (first label = main cohort for sections 4-5)
  --coding-alleles LABEL=alleles.tsv.gz ...     from classify_coding_TRV_alleles.py (±1 bp section)
  --constraint gnomad.v4.1.1.constraint_metrics.tsv.bgz   (LOEUF = lof.oe_ci.upper, MANE genes)
  --out-dir DIR
Filters for sections 4-5: FILTER PASS, single-component TRID, copies >= 2. STR = motif 1-6 bp, VNTR = motif >= 7 bp.
Section 5 per-motif panels use --str-source loci (default TRExplorer). Constraint model: Poisson GLM of n_alt (distinct ALT
alleles with AC > 0) on log(size), log(size)^2, motif class x source, motif class x log(size), GC and purity,
fitted on all coding TRVs
(locus inside one CDS block) of the main cohort; obs/exp is aggregated per gene (multi-gene loci attributed to each gene)
and per LOEUF decile (10 equal-size bins over MANE genes with LOEUF).
Outputs: PNG (deck) + PDF per figure, and summary.json with the numbers quoted on the slides.
"""
import argparse
import json
import os

import matplotlib
import numpy as np
import pandas as pd
import statsmodels.formula.api as smf
import statsmodels.api as sm

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

plt.rcParams["font.family"] = "Arial"
plt.rcParams["pdf.fonttype"] = 42

INK = "#1d1f24"
INK2 = "#5b5e66"
GRID = "#e6e4df"
CONTEXTS = ["coding", "UTR", "intronic", "intergenic"]
CTX_COLORS = {"coding": "#2a78d6", "UTR": "#eb6834", "intronic": "#1baf7a", "intergenic": "#eda100"}
COHORT_COLORS = ["#2a78d6", "#eb6834"]
MOTIF_COLORS = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300", "#4a3aa7", "#e34948"]
METRICS = [("nonref_AF", "Non-ref allele frequency"), ("n_alt", "Distinct ALT alleles per locus"),
           ("max_abs_delta", "Max |ALT − REF| (bp)")]
SUMMARY = {}


def style(ax, grid="y"):
    for s in ["top", "right"]:
        ax.spines[s].set_visible(False)
    for s in ["left", "bottom"]:
        ax.spines[s].set_color(INK2)
    ax.tick_params(colors=INK2, labelsize=8)
    if grid:
        ax.grid(axis=grid, color=GRID, linewidth=0.6)
    ax.set_axisbelow(True)


def save(fig, out_dir, name):
    fig.savefig(os.path.join(out_dir, f"{name}.png"), dpi=200, bbox_inches="tight")
    fig.savefig(os.path.join(out_dir, f"{name}.pdf"), bbox_inches="tight")
    plt.close(fig)


def kv(items):
    return dict(x.split("=", 1) for x in items)


def fmt_n(v):
    return f"{v / 1e6:.2f}M" if v >= 1e6 else (f"{v / 1e3:.0f}k" if v >= 1e3 else f"{v:.0f}")


# ---------- section 1: variant overview ----------
def fig_overview(paths, out_dir):
    tabs = {k: pd.read_csv(v, sep="\t") for k, v in paths.items()}
    order = ["snv", "del_indel", "ins_indel", "del_SV", "ins_SV", "trv"]
    others = sorted({c for t in tabs.values() for c in t["variant_class"]} - set(order))
    order += others
    names = {"snv": "SNV", "del_indel": "Deletion <50 bp", "ins_indel": "Insertion <50 bp", "del_SV": "Deletion SV ≥50 bp",
             "ins_SV": "Insertion SV ≥50 bp", "trv": "Tandem repeat (TRV)"}
    fig, ax = plt.subplots(figsize=(11, 4.2))
    x = np.arange(len(order))
    w = 0.38
    SUMMARY["overview"] = {}
    for i, (label, t) in enumerate(tabs.items()):
        t = t.set_index("variant_class").reindex(order).fillna(0)
        v = t["n_PASS"].to_numpy()
        bars = ax.bar(x + (i - 0.5) * w, v, width=w, color=COHORT_COLORS[i], label=label, edgecolor="white", linewidth=1)
        for b, val in zip(bars, v):
            if val > 0:
                ax.text(b.get_x() + b.get_width() / 2, val * 1.15, fmt_n(val), ha="center", va="bottom", fontsize=7, color=INK)
        SUMMARY["overview"][label] = {c: {"all": int(t.loc[c, "n_all"]), "PASS": int(t.loc[c, "n_PASS"])} for c in order}
    ax.set_yscale("log")
    ax.set_xticks(x)
    ax.set_xticklabels([names.get(c, c) for c in order], rotation=20, ha="right", fontsize=8)
    ax.set_ylabel("PASS variant sites (log)", fontsize=9, color=INK)
    ax.legend(frameon=False, fontsize=9, labelcolor=INK)
    style(ax)
    save(fig, out_dir, "s1_variant_overview")


# ---------- section 2: TRV counts ----------
def motif_class(m):
    return np.select([m <= 6, m <= 10, m <= 20, m <= 50, m <= 100], [m.astype(str), "7–10", "11–20", "21–50", "51–100"], ">100")


def fig_trv_counts(loci, out_dir):
    SUMMARY["trv"] = {}
    fig, axes = plt.subplots(1, 3, figsize=(15, 4.2))
    labels = list(loci)
    # by genic context
    ctx_order = ["coding", "coding_partial", "UTR", "intronic", "intergenic"]
    x = np.arange(len(ctx_order))
    w = 0.38
    for i, lab in enumerate(labels):
        d = loci[lab]
        p = d[d["FILTER"] == "PASS"]
        v = p["genic_context"].value_counts().reindex(ctx_order).fillna(0).to_numpy()
        bars = axes[0].bar(x + (i - 0.5) * w, v, width=w, color=COHORT_COLORS[i], label=lab, edgecolor="white", linewidth=1)
        for b, val in zip(bars, v):
            axes[0].text(b.get_x() + b.get_width() / 2, val * 1.15, f"{fmt_n(val)}\n{100 * val / v.sum():.1f}%",
                         ha="center", va="bottom", fontsize=6.5, color=INK)
        SUMMARY["trv"][lab] = {"n_sites": int(len(d)), "n_PASS": int(len(p)),
                               "by_context": {c: int(n) for c, n in zip(ctx_order, v)},
                               "by_source": p["source"].value_counts().to_dict()}
    axes[0].set_yscale("log")
    axes[0].set_xticks(x)
    axes[0].set_xticklabels(["coding\n(in CDS)", "coding\n(crosses CDS edge)", "UTR", "intronic", "intergenic"], fontsize=8)
    axes[0].set_ylabel("PASS TRV sites (log)", fontsize=9, color=INK)
    axes[0].set_title("By genic context", fontsize=10, color=INK)
    axes[0].legend(frameon=False, fontsize=8, labelcolor=INK)
    style(axes[0])
    # by locus size
    edges = np.logspace(0, 4.5, 46)
    for i, lab in enumerate(labels):
        p = loci[lab][loci[lab]["FILTER"] == "PASS"]
        h = np.histogram(p["size"].clip(1, edges[-1] - 1), edges)[0]
        axes[1].step(edges[:-1], h, where="post", color=COHORT_COLORS[i], linewidth=2, label=lab)
        SUMMARY["trv"][lab]["size_median"] = float(p["size"].median())
        SUMMARY["trv"][lab]["size_p90"] = float(p["size"].quantile(0.9))
    axes[1].set_xscale("log")
    axes[1].set_yscale("log")
    axes[1].set_xlabel("Repeat span (bp; TRID start–end)", fontsize=9, color=INK)
    axes[1].set_ylabel("PASS TRV sites per bin", fontsize=9, color=INK)
    axes[1].set_title("By repeat size", fontsize=10, color=INK)
    style(axes[1], grid="both")
    # by motif size
    mc_order = ["1", "2", "3", "4", "5", "6", "7–10", "11–20", "21–50", "51–100", ">100"]
    x = np.arange(len(mc_order))
    for i, lab in enumerate(labels):
        p = loci[lab][loci[lab]["FILTER"] == "PASS"]
        v = pd.Series(motif_class(p["motif_len"].to_numpy())).value_counts().reindex(mc_order).fillna(0).to_numpy()
        axes[2].bar(x + (i - 0.5) * w, v, width=w, color=COHORT_COLORS[i], label=lab, edgecolor="white", linewidth=0.5)
        SUMMARY["trv"][lab]["by_motif"] = {m: int(n) for m, n in zip(mc_order, v)}
    axes[2].set_yscale("log")
    axes[2].set_xticks(x)
    axes[2].set_xticklabels(mc_order, rotation=40, ha="right", fontsize=8)
    axes[2].set_xlabel("Motif length (bp; shortest TRID motif)", fontsize=9, color=INK)
    axes[2].set_ylabel("PASS TRV sites (log)", fontsize=9, color=INK)
    axes[2].set_title("By motif size", fontsize=10, color=INK)
    style(axes[2])
    fig.tight_layout()
    save(fig, out_dir, "s2_trv_counts")


# ---------- section 3: ±1 bp ----------
def fig_pm1(paths, loci, out_dir):
    SUMMARY["pm1"] = {}
    fig, axes = plt.subplots(1, 3, figsize=(15, 4.0))
    acbins = [("AC = 1", 1, 1), ("AC 2–5", 2, 5), ("AC > 5", 6, 10 ** 9)]
    x = np.arange(len(acbins))
    w = 0.38
    for i, (lab, path) in enumerate(paths.items()):
        a = pd.read_csv(path, sep="\t", usecols=["ID", "AC", "cds_contained", "delta"])
        a = a[a["cds_contained"] & (a["AC"] > 0)]
        nz = a[a["delta"] != 0]
        pm1, m3 = [], []
        for _, lo, hi in acbins:
            s = nz[(nz["AC"] >= lo) & (nz["AC"] <= hi)]
            pm1.append(100 * (s["delta"].abs() == 1).mean())
            m3.append(100 * (s["delta"] % 3 == 0).mean())
        for ax, v in zip(axes[:2], [pm1, m3]):
            bars = ax.bar(x + (i - 0.5) * w, v, width=w, color=COHORT_COLORS[i], label=lab, edgecolor="white", linewidth=1)
            for b, val in zip(bars, v):
                ax.text(b.get_x() + b.get_width() / 2, val + 1, f"{val:.0f}%", ha="center", va="bottom", fontsize=7, color=INK)
        per_site = a.groupby("ID").size()
        SUMMARY["pm1"][lab] = {"pct_pm1_by_AC": dict(zip([b[0] for b in acbins], np.round(pm1, 1))),
                               "pct_mult3_by_AC": dict(zip([b[0] for b in acbins], np.round(m3, 1))),
                               "coding_sites": int(per_site.size), "alt_alleles_per_site_mean": round(float(per_site.mean()), 2),
                               "pm1_alleles": int((nz["delta"].abs() == 1).sum())}
    for ax, title in zip(axes[:2], ["±1 bp share of non-zero length changes", "Multiple-of-3 share of non-zero changes"]):
        ax.set_xticks(x)
        ax.set_xticklabels([b[0] for b in acbins], fontsize=8)
        ax.set_ylim(0, 105)
        ax.set_ylabel("% of carried ALT alleles", fontsize=9, color=INK)
        ax.set_title(title, fontsize=10, color=INK)
        style(ax)
    axes[1].axhline(100 / 3, color=INK2, linestyle="--", linewidth=0.8)
    axes[1].text(2.45, 35, "chance (33%)", fontsize=7, color=INK2, ha="right")
    axes[0].legend(frameon=False, fontsize=8, labelcolor=INK)
    # alt alleles per site, coding vs intergenic STR loci
    for i, (lab, d) in enumerate(loci.items()):
        p = d[(d["FILTER"] == "PASS") & (d["motif_len"] <= 6)]
        v = [p.loc[p["genic_context"] == c, "n_alt"].mean() for c in CONTEXTS]
        axes[2].bar(np.arange(4) + (i - 0.5) * w, v, width=w, color=COHORT_COLORS[i], label=lab, edgecolor="white", linewidth=1)
        SUMMARY["pm1"][lab]["str_n_alt_by_context"] = dict(zip(CONTEXTS, np.round(v, 2)))
    axes[2].set_xticks(np.arange(4))
    axes[2].set_xticklabels(CONTEXTS, fontsize=8)
    axes[2].set_ylabel("Mean distinct ALT alleles per STR locus", fontsize=9, color=INK)
    axes[2].set_title("ALT alleles per STR locus (motif 1–6 bp)", fontsize=10, color=INK)
    style(axes[2])
    fig.tight_layout()
    save(fig, out_dir, "s3_pm1_issue")


# ---------- section 4: mutation rate vs size / motif ----------
def binned_lines(ax, d, xcol, edges, group_col, groups, colors, metric, min_n=100, median=False):
    d = d.assign(bin=pd.cut(d[xcol], edges, right=False))
    for g, color in zip(groups, colors):
        s = d[d[group_col] == g]
        t = s.groupby("bin", observed=True).agg(n=(metric, "size"), x=(xcol, "median"), mean=(metric, "mean"),
                                                sd=(metric, "std"), med=(metric, "median"),
                                                q25=(metric, lambda v: v.quantile(.25)),
                                                q75=(metric, lambda v: v.quantile(.75)))
        t = t[t["n"] >= min_n]
        if t.empty:
            continue
        if median:
            y, lo, hi = t["med"], t["q25"], t["q75"]
        else:
            ci = 1.96 * t["sd"].fillna(0) / np.sqrt(t["n"])
            y, lo, hi = t["mean"], t["mean"] - ci, t["mean"] + ci
        ax.fill_between(t["x"], lo, hi, color=color, alpha=0.15, linewidth=0)
        ax.plot(t["x"], y, color=color, linewidth=1.8, marker="o", markersize=3, label=str(g))


def fig_mutation_rate(d, label, out_dir):
    p = d[(d["FILTER"] == "PASS") & (d["n_components"] == 1) & (d["copies"] >= 2)]
    strs = p[p["motif_len"] <= 6]
    vntr = p[p["motif_len"] >= 7].copy()
    vntr["mclass"] = pd.cut(vntr["motif_len"], [7, 11, 21, 51, 101, np.inf], right=False,
                            labels=["7–10", "11–20", "21–50", "51–100", ">100"])
    SUMMARY["mutation_rate"] = {"n_STR_loci": int(len(strs)), "n_VNTR_loci": int(len(vntr)),
                                "STR_by_source": strs["source"].value_counts().to_dict(),
                                "VNTR_by_source": vntr["source"].value_counts().to_dict()}
    fig, axes = plt.subplots(3, 4, figsize=(17, 11))
    size_edges = [6, 10, 12, 15, 20, 25, 30, 40, 50, 75, 100, 150, 250, 500, 1000, 2500, np.inf]
    copy_edges = [2, 3, 4, 5, 7, 10, 15, 25, np.inf]
    for r, (metric, ylabel) in enumerate(METRICS):
        med = metric == "max_abs_delta"
        # STR vs size, by motif (TRExplorer)
        s_te = strs[strs["source"] == "TRExplorer"]
        binned_lines(axes[r][0], s_te, "size", size_edges, "motif_len", range(1, 7), MOTIF_COLORS[:6], metric, median=med)
        # STR vs motif at fixed copy number (TRExplorer)
        s_te = s_te.assign(cbin=pd.cut(s_te["copies"], copy_edges, right=False))
        cm = plt.get_cmap("Blues")
        bins = [b for b in s_te["cbin"].cat.categories]
        for k, b in enumerate(bins):
            t = s_te[s_te["cbin"] == b].groupby("motif_len").agg(n=(metric, "size"), y=(metric, "median" if med else "mean"))
            t = t[t["n"] >= 100]
            if len(t) >= 2:
                if not np.isfinite(b.right):
                    lab = f"≥{b.left:g}"
                else:
                    lab = f"{b.left:g}" if b.right - b.left == 1 else f"{b.left:g}–{b.right - 1:g}"
                axes[r][1].plot(t.index, t["y"], color=cm(0.3 + 0.7 * k / (len(bins) - 1)), linewidth=1.8, marker="o",
                                markersize=3, label=f"{lab} copies")
        # VNTR vs size, by motif class (all sources)
        binned_lines(axes[r][2], vntr, "size", size_edges, "mclass", ["7–10", "11–20", "21–50", "51–100", ">100"],
                     MOTIF_COLORS[:5], metric, median=med)
        # VNTR vs motif class at fixed size
        vs = vntr.assign(sbin=pd.cut(vntr["size"], [0, 100, 250, 500, 1000, np.inf], right=False,
                                     labels=["<100", "100–249", "250–499", "500–999", "≥1000"]))
        for k, b in enumerate(vs["sbin"].cat.categories):
            t = vs[vs["sbin"] == b].groupby("mclass", observed=True).agg(n=(metric, "size"),
                                                                         y=(metric, "median" if med else "mean"))
            t = t[t["n"] >= 50]
            if len(t) >= 2:
                axes[r][3].plot([str(c) for c in t.index], t["y"], color=cm(0.3 + 0.7 * k / 4), linewidth=1.8, marker="o",
                                markersize=3, label=f"{b} bp")
        for c in range(4):
            ax = axes[r][c]
            if c in (0, 2):
                ax.set_xscale("log")
            if med:
                ax.set_yscale("symlog", linthresh=10)
                ax.set_ylim(bottom=0)
            ax.set_ylabel(ylabel + (" (median)" if med else " (mean)"), fontsize=8.5, color=INK)
            style(ax, grid="both")
    titles = ["STR (1–6 bp motif) vs repeat size\nTRExplorer; lines = motif length",
              "STR vs motif length at fixed copy number\nTRExplorer; lines = copy-number bin",
              "VNTR (≥7 bp motif) vs repeat size\nall catalogs; lines = motif class",
              "VNTR vs motif class at fixed repeat size\nall catalogs; lines = size bin"]
    xlabels = ["Repeat size (bp)", "Motif length (bp)", "Repeat size (bp)", "Motif length class (bp)"]
    for c in range(4):
        axes[0][c].set_title(titles[c], fontsize=9.5, color=INK)
        axes[2][c].set_xlabel(xlabels[c], fontsize=9, color=INK)
        axes[0][c].legend(frameon=False, fontsize=7, labelcolor=INK, loc="upper left")
    fig.suptitle(f"{label}: variability vs repeat size and motif (PASS, single-component loci, ≥2 copies; "
                 f"mean ± 95% CI, max |Δ| as median with IQR)", fontsize=11, color=INK)
    fig.tight_layout()
    save(fig, out_dir, "s4_mutation_rate")


# ---------- section 5: constraint model ----------
def load_loeuf(path):
    c = pd.read_csv(path, sep="\t", usecols=["gene", "gene_id", "mane_select", "lof.oe_ci.upper"],
                    dtype={"gene": str, "gene_id": str, "mane_select": str}, compression="gzip")
    c = c[(c["mane_select"].str.lower() == "true") & c["gene_id"].str.startswith("ENSG") & c["lof.oe_ci.upper"].notna()]
    c = c.drop_duplicates("gene")[["gene", "lof.oe_ci.upper"]].rename(columns={"lof.oe_ci.upper": "LOEUF"})
    c["decile"] = pd.qcut(c["LOEUF"], 10, labels=False) + 1
    return c


def fit_model(d):
    coding = d[(d["genic_context"] == "coding") & (d["FILTER"] == "PASS") & (d["n_components"] == 1) & (d["copies"] >= 2)].copy()
    coding = coding[coding["gc"].notna() & coding["purity"].notna()]
    coding["log_size"] = np.log(coding["size"])
    coding["mclass"] = pd.cut(coding["motif_len"], [0, 1, 2, 3, 4, 5, 6, 12, 24, np.inf],
                              labels=["m1", "m2", "m3", "m4", "m5", "m6", "m7_12", "m13_24", "m25plus"]).astype(str)
    model = smf.glm("n_alt ~ log_size + I(log_size**2) + C(mclass) * C(source) + log_size:C(mclass) + gc + purity",
                    data=coding, family=sm.families.Poisson()).fit()
    coding["expected"] = model.predict(coding)
    dev_expl = 1 - model.deviance / model.null_deviance
    SUMMARY["model"] = {"n_loci": int(len(coding)), "deviance_explained": round(float(dev_expl), 3),
                        "coef": {k: round(float(v), 4) for k, v in model.params.items()
                                 if k in ("log_size", "I(log_size ** 2)", "gc", "purity")},
                        "pvalues": {k: float(f"{v:.2e}") for k, v in model.pvalues.items() if k in ("gc", "purity", "log_size")},
                        "spearman_obs_exp": round(float(coding[["n_alt", "expected"]].corr(method="spearman").iloc[0, 1]), 3)}
    return coding, model


def gene_obs_exp(coding, loeuf, rng, n_boot=500):
    g = coding.assign(gene=coding["genes"].str.split(",")).explode("gene").merge(loeuf, on="gene", how="inner")
    per_gene = g.groupby(["gene", "decile", "LOEUF"], as_index=False).agg(obs=("n_alt", "sum"), exp=("expected", "sum"),
                                                                          n_loci=("ID", "size"))
    per_gene["oe"] = per_gene["obs"] / per_gene["exp"]
    rows = []
    for dec, s in per_gene.groupby("decile"):
        pooled = s["obs"].sum() / s["exp"].sum()
        idx = rng.integers(0, len(s), size=(n_boot, len(s)))
        boots = s["obs"].to_numpy()[idx].sum(1) / s["exp"].to_numpy()[idx].sum(1)
        rows.append({"decile": int(dec), "LOEUF_median": float(s["LOEUF"].median()), "n_genes": len(s),
                     "n_loci": int(s["n_loci"].sum()), "obs_exp": pooled, "lo": float(np.percentile(boots, 2.5)),
                     "hi": float(np.percentile(boots, 97.5))})
    return per_gene, pd.DataFrame(rows)


def draw_oe(ax_box, ax_line, per_gene, dec_tab, color, label, min_exp=3):
    pg = per_gene[per_gene["exp"] >= min_exp]
    data = [pg.loc[pg["decile"] == k, "oe"].to_numpy() for k in range(1, 11)]
    bp = ax_box.boxplot([np.log2(v[v > 0]) if (v > 0).any() else [np.nan] for v in data], positions=range(1, 11), widths=0.6,
                        showfliers=False, patch_artist=True)
    for b in bp["boxes"]:
        b.set(facecolor=color, alpha=0.35, edgecolor=color)
    for k in ("medians", "whiskers", "caps"):
        for el in bp[k]:
            el.set(color=color)
    ax_box.axhline(0, color=INK2, linestyle="--", linewidth=0.8)
    ax_box.set_xlabel("LOEUF decile (1 = most constrained)", fontsize=8.5, color=INK)
    ax_box.set_ylabel(f"log2 obs/exp per gene (genes with ≥{min_exp} expected)", fontsize=8.5, color=INK)
    ax_line.fill_between(dec_tab["decile"], dec_tab["lo"], dec_tab["hi"], color=color, alpha=0.18, linewidth=0)
    ax_line.plot(dec_tab["decile"], dec_tab["obs_exp"], color=color, linewidth=2, marker="o", markersize=4, label=label)
    ax_line.axhline(1, color=INK2, linestyle="--", linewidth=0.8)
    ax_line.set_xticks(range(1, 11))
    ax_line.set_xlabel("LOEUF decile (1 = most constrained)", fontsize=8.5, color=INK)
    ax_line.set_ylabel("Pooled obs/exp ALT alleles (95% CI)", fontsize=8.5, color=INK)
    style(ax_box)
    style(ax_line, grid="both")


def fig_model(coding, model, loeuf, rng, out_dir):
    per_gene, dec = gene_obs_exp(coding, loeuf, rng)
    SUMMARY["model"]["all_coding_obs_exp_by_decile"] = dec.round(3).to_dict("records")
    SUMMARY["model"]["n_genes"] = int(per_gene["gene"].nunique())
    fig, axes = plt.subplots(1, 3, figsize=(15, 4.2))
    # calibration
    q = pd.qcut(coding["expected"], 20, duplicates="drop")
    cal = coding.groupby(q, observed=True).agg(exp=("expected", "mean"), obs=("n_alt", "mean"))
    axes[0].plot(cal["exp"], cal["obs"], color="#2a78d6", marker="o", markersize=4, linewidth=1.8)
    lim = [0, max(cal["exp"].max(), cal["obs"].max()) * 1.05]
    axes[0].plot(lim, lim, color=INK2, linestyle="--", linewidth=0.8)
    axes[0].set_xlabel("Expected ALT alleles (model; 20 quantile bins)", fontsize=8.5, color=INK)
    axes[0].set_ylabel("Observed ALT alleles (mean)", fontsize=8.5, color=INK)
    axes[0].set_title(f"Calibration on {len(coding):,} coding TRVs\n"
                      f"deviance explained {SUMMARY['model']['deviance_explained']:.2f}", fontsize=10, color=INK)
    style(axes[0], grid="both")
    draw_oe(axes[1], axes[2], per_gene, dec, "#2a78d6", "all coding TRVs")
    axes[1].set_title(f"Per-gene obs/exp ({int((per_gene['exp'] >= 3).sum()):,} genes with ≥3 expected)", fontsize=10, color=INK)
    axes[2].set_title("Pooled obs/exp by LOEUF decile", fontsize=10, color=INK)
    fig.tight_layout()
    save(fig, out_dir, "s5_model")
    per_gene.to_csv(os.path.join(out_dir, "s5_model.per_gene_obs_exp.tsv"), sep="\t", index=False, float_format="%.4g")


def fig_motif(d, coding, loeuf, motif, source, rng, out_dir):
    p = d[(d["FILTER"] == "PASS") & (d["n_components"] == 1) & (d["copies"] >= 2) & (d["motif_len"] == motif)
          & (d["source"] == source) & d["genic_context"].isin(CONTEXTS)]
    s = {"n_loci": int(len(p)), "by_context": p["genic_context"].value_counts().reindex(CONTEXTS).fillna(0).astype(int).to_dict()}
    fig, axes = plt.subplots(2, 3, figsize=(16, 8.6))
    # (a) counts
    v = [s["by_context"][c] for c in CONTEXTS]
    bars = axes[0][0].bar(CONTEXTS, v, color=[CTX_COLORS[c] for c in CONTEXTS], edgecolor="white", linewidth=1)
    for b, val in zip(bars, v):
        axes[0][0].text(b.get_x() + b.get_width() / 2, val * 1.15, f"{val:,}\n{100 * val / max(sum(v), 1):.1f}%",
                        ha="center", va="bottom", fontsize=7, color=INK)
    axes[0][0].set_yscale("log")
    axes[0][0].set_ylim(top=max(v) * 8)
    axes[0][0].set_ylabel("Loci (log)", fontsize=8.5, color=INK)
    axes[0][0].set_title("a  Loci by genic context", fontsize=10, color=INK, loc="left")
    style(axes[0][0])
    # (b, c) mutation rate vs copies
    edges = [2, 3, 4, 5, 6, 8, 10, 15, 25, np.inf]
    for ax, (metric, ylabel), tag in zip([axes[0][1], axes[0][2]], METRICS[:2], ["b", "c"]):
        binned_lines(ax, p, "copies", edges, "genic_context", CONTEXTS, [CTX_COLORS[c] for c in CONTEXTS], metric, min_n=20)
        ax.set_xscale("log")
        ax.set_xticks([2, 3, 4, 6, 10, 20, 40])
        ax.xaxis.set_major_formatter(matplotlib.ticker.ScalarFormatter())
        ax.xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
        ax.set_xlabel("Copy number", fontsize=8.5, color=INK)
        ax.set_ylabel(ylabel + " (mean ± 95% CI)", fontsize=8.5, color=INK)
        ax.set_title(f"{tag}  {ylabel} vs copy number", fontsize=10, color=INK, loc="left")
        style(ax, grid="both")
    axes[0][1].legend(frameon=False, fontsize=7.5, labelcolor=INK, loc="upper left")
    s["nonref_AF_by_context"] = {c: round(float(p.loc[p["genic_context"] == c, "nonref_AF"].mean()), 4) for c in CONTEXTS}
    s["n_alt_by_context"] = {c: round(float(p.loc[p["genic_context"] == c, "n_alt"].mean()), 3) for c in CONTEXTS}
    # (d) max |delta| distribution
    ax = axes[1][0]
    for c in CONTEXTS:
        v = np.sort(p.loc[p["genic_context"] == c, "max_abs_delta"].to_numpy())
        if len(v):
            ax.step(v + 1, np.arange(1, len(v) + 1) / len(v), where="post", color=CTX_COLORS[c], linewidth=1.8, label=c)
    ax.set_xscale("log")
    ticks = [1, 2, 4, 7, 11, 31, 101, 1001]
    ax.set_xticks(ticks)
    ax.set_xticklabels([str(t - 1) for t in ticks])
    ax.set_xlabel("Max |ALT − REF| (bp)", fontsize=8.5, color=INK)
    ax.set_ylabel("Cumulative fraction of loci", fontsize=8.5, color=INK)
    ax.set_title("d  Max size difference ALT vs REF", fontsize=10, color=INK, loc="left")
    ax.legend(frameon=False, fontsize=7.5, labelcolor=INK, loc="lower right")
    style(ax, grid="both")
    s["max_abs_delta_median_by_context"] = {c: float(p.loc[p["genic_context"] == c, "max_abs_delta"].median()) for c in CONTEXTS}
    s["pct_max_abs_delta_ge_motif"] = {c: round(100 * float((p.loc[p["genic_context"] == c, "max_abs_delta"] >= motif).mean()), 1)
                                       for c in CONTEXTS}
    # (e, f) obs/exp by LOEUF decile
    cm = coding[(coding["motif_len"] == motif) & (coding["source"] == source)]
    per_gene, dec = gene_obs_exp(cm, loeuf, rng)
    draw_oe(axes[1][1], axes[1][2], per_gene, dec, CTX_COLORS["coding"], f"motif {motif} bp")
    axes[1][1].set_title(f"e  Per-gene obs/exp ({int((per_gene['exp'] >= 3).sum()):,} genes with ≥3 expected)",
                         fontsize=10, color=INK, loc="left")
    axes[1][2].set_title("f  Pooled obs/exp vs LOEUF decile", fontsize=10, color=INK, loc="left")
    for k, r in dec.iterrows():
        axes[1][2].text(r["decile"], r["hi"], f"{int(r['n_loci'])}", fontsize=6, color=INK2, ha="center", va="bottom")
    s["obs_exp_by_decile"] = dec.round(3).to_dict("records")
    top = per_gene.assign(excess=per_gene["obs"] - per_gene["exp"]).nlargest(5, "excess")
    s["top_excess_genes"] = top[["gene", "decile", "obs", "exp", "n_loci"]].round(1).to_dict("records")
    fig.suptitle(f"{motif}-bp motif STRs ({source}, PASS, ≥2 copies): {len(p):,} loci; coding = inside one CDS block; "
                 f"numbers above f = coding loci per decile", fontsize=11, color=INK)
    fig.tight_layout()
    save(fig, out_dir, f"s5_motif{motif}")
    SUMMARY.setdefault("motif", {})[str(motif)] = s


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--overview", nargs="+", required=True)
    p.add_argument("--loci", nargs="+", required=True)
    p.add_argument("--coding-alleles", nargs="+", required=True)
    p.add_argument("--constraint", required=True)
    p.add_argument("--str-source", default="TRExplorer")
    p.add_argument("--out-dir", required=True)
    p.add_argument("--seed", type=int, default=1)
    args = p.parse_args()
    os.makedirs(args.out_dir, exist_ok=True)
    rng = np.random.default_rng(args.seed)

    loci = {k: pd.read_csv(v, sep="\t", dtype={"genes": str}) for k, v in kv(args.loci).items()}
    main_label = next(iter(loci))
    fig_overview(kv(args.overview), args.out_dir)
    fig_trv_counts(loci, args.out_dir)
    fig_pm1(kv(args.coding_alleles), loci, args.out_dir)
    fig_mutation_rate(loci[main_label], main_label, args.out_dir)
    loeuf = load_loeuf(args.constraint)
    SUMMARY["loeuf_decile_cutoffs"] = [round(float(x), 3) for x in loeuf["LOEUF"].quantile(np.linspace(0.1, 0.9, 9))]
    coding, model = fit_model(loci[main_label])
    with open(os.path.join(args.out_dir, "s5_model.summary.txt"), "w") as f:
        f.write(model.summary().as_text())
    fig_model(coding, model, loeuf, rng, args.out_dir)
    for m in range(1, 7):
        fig_motif(loci[main_label], coding, loeuf, m, args.str_source, rng, args.out_dir)
    with open(os.path.join(args.out_dir, "summary.json"), "w") as f:
        json.dump(SUMMARY, f, indent=1, default=lambda o: o.item() if hasattr(o, "item") else str(o))
    print(json.dumps(SUMMARY, indent=1, default=lambda o: o.item() if hasattr(o, "item") else str(o))[:6000])


if __name__ == "__main__":
    main()
