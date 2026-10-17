#!/usr/bin/env python3
"""Learn from gnomAD short-read (SR) STR loci how far above the long-read (LR) median the SR population maximum reaches, then
flag rare LR alleles beyond that predicted SR reach at every STR, including the ones without SR data.

LR sites (merged TRV table): FILTER PASS, summed AC > 0, one motif (compound TRVs skipped), motif length in
[--min-motif, --max-motif]. Allele tract length (bp) = ref_len (REF) or ref_len + alt_size_diff (ALT); AC-weighted (REF by
AN - sum AC); median in bp.
Training: SR loci with a clean TRV match (--sr-loci from split_sr_matched_strs.py, match_flag ok). Target = SR reach in bp,
log(1 + (SR_max - SR_median) x motif length) (SR repeat-unit counts from --sr-allele-dist). Model: linear quantile regression at
--quantile (default 0.95) on --features (default motif_len). Candidate feature sets are compared by leave-one-locus-out (LOO)
coverage and pinball loss (.model_selection.tsv); LR spread features (IQR, p90 excess, median, n distinct alleles) do not
improve on motif length alone in hprc_hgsvc.
Application: per site T = LR median + predicted reach (bp). Extreme outlier allele: length > T and AC <= --max-ac (rare).
Score = (length - median) / predicted reach. Each outlier allele's sequence (VCF allele without the padding base) is checked
for repeat content: motif_fraction = fraction of the allele in tandem runs (>= 3 copies, >= 12 bp) of the site motif (best
rotation / strand); best_repeat_motif / fraction = the most covering 1-6 bp repeat in the allele. allele_class:
motif_expansion (motif_fraction >= --min-motif-fraction), other_repeat (a different short repeat covers >= that fraction, e.g.
RFC1 AAGGG at an AAAAG site), homopolymer_run (the covering repeat is a single base; often an assembly / sequencing
artifact), non_repeat_insertion (none; e.g. an SV, MEI or long-period VNTR recorded on the STR). Sites whose motif length
is outside the training motif lengths are flagged (extrapolated).
Outputs (prefix --out-prefix):
  .model_selection.tsv  LOO coverage / pinball loss per candidate feature set
  .model.tsv            coefficients of the chosen model
  .training.tsv         training loci: SR max / median, LOO predicted reach and threshold, LR alleles beyond it
  .sites.tsv.gz         per STR site: motif, median, threshold, n outlier alleles, outlier AC, max length, max score,
                        n outlier alleles per allele_class
  .outliers.tsv.gz      per outlier allele: site, motif, genic_context, gene, length (bp), units, AC, score,
                        motif_fraction, best_repeat_motif / fraction, allele_class
  .pdf                  calibration (LOO predicted vs observed SR reach), training loci LR max vs threshold, genome-wide rate

Usage:
  sr_envelope_outliers.py --trv gnomAD_LR.hprc_hgsvc.TRV.bed.gz --sr-loci STR_SRmatch.hgsvc_hprc.loci.tsv \\
      --sr-allele-dist STR_SRmatch.sr.allele_dist.tsv --out-prefix sr_envelope/hprc_hgsvc
"""
import argparse
import collections
import re

import matplotlib
import numpy as np
import pandas as pd
import statsmodels.api as sm

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

plt.rcParams["font.family"] = "Arial"
plt.rcParams["pdf.fonttype"] = 42

INK = "#1d1f24"
INK2 = "#5b5e66"
GRID = "#e6e4df"
BLUE = "#2a78d6"
RED = "#d03b3b"
GREY = "#8a8984"
CLASSES = ["motif_expansion", "other_repeat", "homopolymer_run", "non_repeat_insertion"]
CANDIDATES = [[], ["motif_len"], ["log_iqr"], ["log_p90_excess"], ["motif_len", "log_iqr"], ["motif_len", "log_p90_excess"],
              ["motif_len", "log_median", "log_iqr", "log_p90_excess", "log_n_distinct"]]


def wq(v, w, qs):
    o = np.argsort(v, kind="stable")
    v, c = v[o], np.cumsum(w[o])
    return [v[min(np.searchsorted(c, q * c[-1]), len(v) - 1)] for q in qs]


def motif_fraction(seq, motif, min_run_bp=12):
    """Fraction of the sequence in tandem runs of the motif (any rotation, either strand; a run = >= 3 copies and
    >= min_run_bp, so scattered copies do not count)."""
    if not seq:
        return np.nan
    comp = str.maketrans("ACGT", "TGCA")
    k = len(motif)
    n = max(3, -(-min_run_bp // k))
    cover = 0
    for strand in {motif, motif.translate(comp)[::-1]}:
        for i in range(k):
            r = strand[i:] + strand[:i]
            cover = max(cover, sum(len(m) for m in re.findall(f"(?:{r}){{{n},}}", seq.upper())))
    return cover / len(seq)


def best_repeat(seq, max_k=6):
    """Most covering short repeat in the sequence: (motif, fraction) over k = 1..max_k, using the most frequent k-mer per k."""
    best = ("", 0.0)
    seq = seq.upper()
    for k in range(1, max_k + 1):
        if len(seq) < 2 * k:
            break
        km = collections.Counter(seq[i:i + k] for i in range(len(seq) - k + 1)).most_common(1)[0][0]
        f = motif_fraction(seq, km)
        if f > best[1] + 0.05:
            best = (km, f)
    return best


def site_alleles(path, min_motif, max_motif, keep_ids=None, with_seq=False, chunksize=200000):
    """Yield (site row, motif length, allele lengths bp, allele counts, allele sequences or None) for qualifying sites."""
    cols = ["#chrom", "start", "end", "ID", "FILTER", "motifs", "ref_len", "alt_size_diff", "AC", "AN", "genic_context", "genes"]
    cols += ["REF", "ALT"] if with_seq else []
    for ch in pd.read_csv(path, sep="\t", usecols=cols, dtype=str, chunksize=chunksize):
        ch = ch[ch["FILTER"] == "PASS"]
        for r in ch.itertuples(index=False):
            if keep_ids is not None:
                if r.ID not in keep_ids:
                    continue
                k = min(len(m) for m in r.motifs.split(","))
            else:
                if "," in r.motifs:
                    continue
                k = len(r.motifs)
                if not min_motif <= k <= max_motif:
                    continue
            ac = np.array([int(x) if x not in (".", "") else 0 for x in r.AC.split(",")])
            if ac.sum() == 0:
                continue
            ref_len = int(r.ref_len)
            lens = np.maximum(np.concatenate([[ref_len], [ref_len + int(x) for x in r.alt_size_diff.split(",")]]), 0)
            counts = np.concatenate([[max(int(r.AN) - ac.sum(), 0)], ac])
            keep = counts > 0
            seqs = np.array([r.REF[1:]] + [a[1:] for a in r.ALT.split(",")], dtype=object)[keep] if with_seq else None
            yield r, k, lens[keep].astype(float), counts[keep], seqs


def features(lens, counts, k):
    u = lens / k
    q1, med, q3, p90 = wq(u, counts, [0.25, 0.5, 0.75, 0.9])
    return {"motif_len": k, "median_bp": med * k, "median_units": med, "IQR_units": q3 - q1, "n_alleles": int(counts.sum()),
            "n_distinct": len(u), "max_bp": lens.max(), "log_median": np.log1p(med), "log_iqr": np.log1p(q3 - q1),
            "log_p90_excess": np.log1p(max(p90 - med, 0)), "log_n_distinct": np.log(len(u))}


def design(df, cols):
    return sm.add_constant(df[cols], has_constant="add") if cols else pd.DataFrame({"const": np.ones(len(df))}, index=df.index)


def loo_predict(y, X, q):
    pred = []
    for i in range(len(y)):
        m = np.ones(len(y), bool)
        m[i] = False
        f = sm.QuantReg(y[m], X[m]).fit(q=q, max_iter=5000)
        pred.append(float(np.asarray(f.predict(X.iloc[[i]]))[0]))
    return np.array(pred)


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--trv", required=True)
    p.add_argument("--sr-loci", required=True)
    p.add_argument("--sr-allele-dist", required=True)
    p.add_argument("--features", nargs="*", default=["motif_len"])
    p.add_argument("--quantile", type=float, default=0.95)
    p.add_argument("--min-motif", type=int, default=2)
    p.add_argument("--max-motif", type=int, default=6)
    p.add_argument("--max-ac", type=int, default=5, help="Outlier alleles must have AC <= this (rare)")
    p.add_argument("--min-motif-fraction", type=float, default=0.5,
                   help="Repeat coverage needed for motif_expansion / other_repeat")
    p.add_argument("--label", default="")
    p.add_argument("--out-prefix", required=True)
    args = p.parse_args()

    loci = pd.read_csv(args.sr_loci, sep="\t")
    loci = loci[loci["match_flag"] == "ok"].set_index("TRV")
    sr = pd.read_csv(args.sr_allele_dist, sep="\t")
    sr = sr[sr["allele_count"] > 0]
    sr_stats = {gid: (wq(g["repunit_count"].to_numpy(float), g["allele_count"].to_numpy(float), [0.5])[0],
                      g["repunit_count"].max(), len(loci.loc[loci["id"] == gid, "repeat_unit"].iloc[0])
                      if (loci["id"] == gid).any() else np.nan) for gid, g in sr.groupby("id")}

    # training loci (any motif structure; the matched TRV is the gnomAD locus)
    rows, train_alleles = [], {}
    for r, k, lens, c, _ in site_alleles(args.trv, 0, 0, keep_ids=set(loci.index)):
        f = features(lens, c, k)
        f["ID"] = r.ID
        rows.append(f)
        train_alleles[r.ID] = (lens, c)
    tr = pd.DataFrame(rows)
    tr["gnomAD_id"] = tr["ID"].map(loci["id"])
    tr["SR_median"] = tr["gnomAD_id"].map(lambda g: sr_stats[g][0])
    tr["SR_max"] = tr["gnomAD_id"].map(lambda g: sr_stats[g][1])
    tr["SR_unit_len"] = tr["gnomAD_id"].map(lambda g: sr_stats[g][2])
    tr = tr[(tr["SR_unit_len"] == tr["motif_len"]) & tr["motif_len"].between(args.min_motif, args.max_motif)]
    tr = tr.reset_index(drop=True)
    tr["SR_reach_bp"] = (tr["SR_max"] - tr["SR_median"]).clip(lower=0) * tr["motif_len"]
    tr["y"] = np.log1p(tr["SR_reach_bp"])

    sel = []
    for cols in CANDIDATES:
        pred = loo_predict(tr["y"], design(tr, cols), args.quantile)
        r = tr["y"] - pred
        sel.append({"features": "+".join(cols) or "intercept_only", "LOO_coverage": (tr["y"] <= pred).mean(),
                    "LOO_pinball": np.mean(np.maximum(args.quantile * r, (args.quantile - 1) * r)),
                    "median_pred_reach_bp": float(np.median(np.expm1(pred)))})
    sel = pd.DataFrame(sel)
    sel.to_csv(f"{args.out_prefix}.model_selection.tsv", sep="\t", index=False, float_format="%.4g")
    print(sel.round(3).to_string(index=False))

    X = design(tr, args.features)
    fit = sm.QuantReg(tr["y"], X).fit(q=args.quantile, max_iter=5000)
    coef = pd.DataFrame({"term": fit.params.index, "coef": fit.params.values})
    coef.to_csv(f"{args.out_prefix}.model.tsv", sep="\t", index=False, float_format="%.4g")
    tr["pred_reach_bp_loo"] = np.expm1(loo_predict(tr["y"], X, args.quantile))
    tr["threshold_bp_loo"] = tr["median_bp"] + tr["pred_reach_bp_loo"]
    tr["SR_max_within_reach"] = tr["SR_reach_bp"] <= tr["pred_reach_bp_loo"]
    n_b, ac_b = [], []
    for i, t in zip(tr["ID"], tr["threshold_bp_loo"]):
        lens, c = train_alleles[i]
        m = (lens > t) & (c <= args.max_ac)
        n_b.append(int(m.sum()))
        ac_b.append(int(c[m].sum()))
    tr["n_rare_LR_alleles_beyond"], tr["AC_beyond"] = n_b, ac_b
    tr = tr.sort_values(["AC_beyond", "gnomAD_id"], ascending=[False, True])
    tr.drop(columns=[c for c in tr.columns if c.startswith("log_")]).to_csv(
        f"{args.out_prefix}.training.tsv", sep="\t", index=False, float_format="%.4g")
    # with motif length as the only feature the reach is fixed per motif length
    reach = {}
    if args.features == ["motif_len"]:
        for k in range(args.min_motif, args.max_motif + 1):
            reach[k] = float(np.expm1(np.asarray(fit.predict(design(pd.DataFrame({"motif_len": [k]}), ["motif_len"])))[0]))
    print(coef.to_string(index=False))
    if reach:
        print("predicted SR reach (bp) by motif length: " + ", ".join(f"{k}: {v:.0f}" for k, v in reach.items()))
    print(f"training loci {len(tr)}; LOO: SR max within predicted reach at {tr['SR_max_within_reach'].mean():.0%}")
    print(tr[["gnomAD_id", "motif_len", "median_bp", "SR_reach_bp", "pred_reach_bp_loo", "max_bp",
              "n_rare_LR_alleles_beyond", "AC_beyond"]].head(15).round(0).to_string(index=False))

    # genome-wide application
    train_k = set(tr["motif_len"])
    lo, hi = tr[args.features].min(), tr[args.features].max()
    site_rows, outs = [], []
    for r, k, lens, c, seqs in site_alleles(args.trv, args.min_motif, args.max_motif, with_seq=True):
        f = features(lens, c, k)
        if reach:
            rch = reach[k]
        else:
            xf = pd.DataFrame([f])[args.features].clip(lo, hi, axis=1) if args.features else pd.DataFrame(index=[0])
            rch = float(np.expm1(np.asarray(fit.predict(design(xf, args.features)))[0]))
        t = f["median_bp"] + rch
        m = (lens > t) & (c <= args.max_ac)
        for L, cc, sq in zip(lens[m], c[m], seqs[m]):
            mf = motif_fraction(sq, r.motifs)
            bm, bf = best_repeat(sq)
            if mf >= args.min_motif_fraction:
                cls = "motif_expansion"
            elif bf >= args.min_motif_fraction:
                cls = "homopolymer_run" if len(bm) == 1 else "other_repeat"
            else:
                cls = "non_repeat_insertion"
            outs.append({"ID": r.ID, "chrom": r[0], "start": int(r.start), "motif": r.motifs, "motif_len": k,
                         "genic_context": r.genic_context, "genes": r.genes, "median_bp": f["median_bp"],
                         "threshold_bp": round(t, 1), "length_bp": L, "units": round(L / k, 2), "AC": int(cc),
                         "score": round((L - f["median_bp"]) / rch, 2), "allele_seq_len": len(sq),
                         "motif_fraction": round(mf, 3), "best_repeat_motif": bm, "best_repeat_fraction": round(bf, 3),
                         "allele_class": cls})
        site_rows.append({"ID": r.ID, "chrom": r[0], "start": int(r.start), "end": int(r.end), "motif": r.motifs,
                          "motif_len": k, "genic_context": r.genic_context, "genes": r.genes, "n_alleles": f["n_alleles"],
                          "median_bp": f["median_bp"], "IQR_units": f["IQR_units"], "max_bp": f["max_bp"],
                          "threshold_bp": round(t, 1), "n_outlier_alleles": int(m.sum()), "outlier_AC": int(c[m].sum()),
                          "extrapolated": k not in train_k})
    sites = pd.DataFrame(site_rows)
    outs = pd.DataFrame(outs).sort_values("score", ascending=False)
    sites = sites.join(outs.groupby("ID")["score"].max().rename("max_score"), on="ID")
    for cls in CLASSES:
        n = outs[outs["allele_class"] == cls].groupby("ID").size().rename(f"n_{cls}")
        sites = sites.join(n, on="ID")
        sites[f"n_{cls}"] = sites[f"n_{cls}"].fillna(0).astype(int)
    sites.to_csv(f"{args.out_prefix}.sites.tsv.gz", sep="\t", index=False, float_format="%.4g")
    outs.to_csv(f"{args.out_prefix}.outliers.tsv.gz", sep="\t", index=False)
    hit = sites["n_outlier_alleles"] > 0
    print(f"genome-wide: {hit.sum():,} of {len(sites):,} single-motif STR sites (motif {args.min_motif}-{args.max_motif} bp) "
          f"with a rare allele (AC <= {args.max_ac}) beyond the SR reach ({hit.mean():.2%}); {len(outs):,} alleles")
    summ = sites.groupby("motif_len").agg(n_sites=("ID", "size"), n_hit=("n_outlier_alleles", lambda x: int((x > 0).sum())))
    summ["pct"] = (100 * summ["n_hit"] / summ["n_sites"]).round(3)
    print(summ.to_string())
    tab = {"any_outlier": sites[hit].groupby("genic_context").size()}
    for cls in CLASSES:
        tab[cls] = sites[sites[f"n_{cls}"] > 0].groupby("genic_context").size()
        print(f"{cls}: {int((sites[f'n_{cls}'] > 0).sum()):,} sites, {int((outs['allele_class'] == cls).sum()):,} alleles")
    print("sites by genic context (a site can have alleles of several classes):")
    print(pd.DataFrame(tab).fillna(0).astype(int).to_string())

    fig, axes = plt.subplots(1, 3, figsize=(17, 5))
    a = axes[0]
    a.scatter(tr["pred_reach_bp_loo"], tr["SR_reach_bp"], s=18, color=BLUE)
    for r in tr[~tr["SR_max_within_reach"]].itertuples():
        a.annotate(r.gnomAD_id, (r.pred_reach_bp_loo, r.SR_reach_bp), fontsize=6.5, color=INK2, xytext=(2, 2),
                   textcoords="offset points")
    mx = max(tr["pred_reach_bp_loo"].max(), tr["SR_reach_bp"].max()) * 1.3
    a.plot([5, mx], [5, mx], color=GREY, linestyle="--", linewidth=0.8)
    a.set_xscale("log")
    a.set_yscale("log")
    a.set_xlabel("Predicted SR reach above median (bp; leave-one-out)", fontsize=9, color=INK)
    a.set_ylabel("Observed SR max − SR median (bp)", fontsize=9, color=INK)
    a.set_title(f"Calibration ({len(tr)} SR loci, q = {args.quantile}): observed ≤ predicted at "
                f"{tr['SR_max_within_reach'].mean():.0%}", fontsize=9, color=INK)
    a = axes[1]
    col = np.where(tr["n_rare_LR_alleles_beyond"] > 0, RED, BLUE)
    a.scatter(tr["threshold_bp_loo"], tr["max_bp"], s=18, color=col)
    for r in tr[tr["n_rare_LR_alleles_beyond"] > 0].itertuples():
        a.annotate(f"{r.gnomAD_id} ({r.n_rare_LR_alleles_beyond})", (r.threshold_bp_loo, r.max_bp), fontsize=6.5,
                   color=INK2, xytext=(2, 2), textcoords="offset points")
    mx = max(tr["threshold_bp_loo"].max(), tr["max_bp"].max()) * 1.3
    a.plot([10, mx], [10, mx], color=GREY, linestyle="--", linewidth=0.8)
    a.set_xscale("log")
    a.set_yscale("log")
    a.set_xlabel("Threshold = LR median + predicted SR reach (bp; leave-one-out)", fontsize=9, color=INK)
    a.set_ylabel("Largest LR allele (bp)", fontsize=9, color=INK)
    a.set_title(f"Training loci: red = rare LR alleles (AC ≤ {args.max_ac}) beyond SR reach (n)", fontsize=9, color=INK)
    a = axes[2]
    cls_col = dict(zip(CLASSES, [BLUE, "#eb6834", "#eda100", GREY]))
    bottom = np.zeros(len(summ))
    for cls, colr in cls_col.items():
        n = sites[sites[f"n_{cls}"] > 0].groupby("motif_len").size().reindex(summ.index).fillna(0)
        a.bar(summ.index, n.to_numpy(), bottom=bottom, color=colr, label=f"{cls.replace('_', ' ')} ({int(n.sum()):,} sites)")
        bottom += n.to_numpy()
    a.set_yscale("log")
    a.set_xlabel("Motif length (bp); motif lengths 2 and 4 not in training", fontsize=9, color=INK)
    a.set_ylabel("Sites with a rare allele beyond SR reach, by allele class (log)", fontsize=9, color=INK)
    a.legend(frameon=False, fontsize=7, labelcolor=INK, loc="upper right", bbox_to_anchor=(1, 0.98))
    a.set_ylim(1, bottom.max() * 30)
    a.set_title(f"Genome-wide: {hit.sum():,} of {len(sites):,} single-motif STR sites", fontsize=9, color=INK)
    for a in axes:
        for sp in ["top", "right"]:
            a.spines[sp].set_visible(False)
        a.tick_params(colors=INK2, labelsize=8)
        a.grid(color=GRID, linewidth=0.6)
        a.set_axisbelow(True)
    fig.suptitle(f"{args.label + ': ' if args.label else ''}extreme LR STR alleles beyond the SR-calibrated reach", fontsize=11,
                 color=INK)
    fig.tight_layout()
    fig.savefig(f"{args.out_prefix}.pdf")
    plt.close(fig)


if __name__ == "__main__":
    main()
