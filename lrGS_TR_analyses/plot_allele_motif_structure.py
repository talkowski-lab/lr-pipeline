#!/usr/bin/env python3
"""Repeat-unit structure of every allele at one TRV: discover the repeat units, paint each allele by motif, list AC.

Input: a TRV table with columns ID, REF, ALT, AC, AN (merged TRV table or VEP-parsed TRV bed; bgzipped + tabix-indexed) and a
TRV ID. Alleles = REF + carried ALT alleles (AC > 0); the VCF padding base is removed; sequences on the reference + strand.
Motif discovery: primitive k-mers, k in --k (default: the site's catalog motif length +/- 1, e.g. 4-6 for a 5-bp motif), that
occur in tandem (seq[i:i+k] == seq[i+k:i+2k]) in any allele; rotations are merged; motifs covering >= --min-frac of all
tandem bases (summed over distinct alleles) are kept, at most --max-motifs; other tandem k-mers are pooled as "other repeat".
Motifs are named by their alphabetically smallest rotation (e.g. AAGGG for GGAAG / GAAGG / ...).
Allele decomposition (left to right): at each position the kept motif (any rotation, any length) with the longest run of
>= --min-copies copies is painted; otherwise the base is "other repeat" (inside a run of another k-mer) or "non-repeat".
Outputs (prefix --out-prefix):
  .pdf          one row per allele (REF first, then ALT alleles grouped by dominant motif and sorted by length), x = bp,
                coloured by motif; right: allele name, length and AC
  .alleles.tsv  per allele: name, AC, length, bp per motif, dominant motif, sequence (+ carrier haplotypes per ancestry)
Layout: allele name, length and AC as row labels on the left; with --genotypes and --ancestry, the left 20% panel shows the
genetic ancestry of the haplotypes carrying each allele (bar length = AC, split into the alleles from each ancestry; gnomAD
GEN_ANC_COLORS) and the right 80% the
repeat-unit structure. With --sr-allele-sizes and --sr-id, a grey histogram of the gnomAD short-read allele lengths (repeat-unit
count x unit length, bp; summed over cohorts and repeat-unit variants) sits above the structure panel on the same x axis;
--top-trv LABEL=path=TRV_ID adds the AC-weighted allele-length histogram of another long-read cohort (e.g. AoU, light
orange) below it.
  .motifs.tsv   per motif: rotations seen, tandem bp in alleles, AC-weighted bp, n alleles containing it

Usage:
  plot_allele_motif_structure.py --trv gnomAD_LR.hprc_hgsvc.TRV.bed.gz --trv-id chr4-39348424-TRV-55 --label RFC1 \\
      --genotypes STR_SRmatch.hgsvc_hprc.genotypes.tsv --ancestry ancestry_label.HGSVC_HPRC.tsv.gz \\
      --sr-allele-sizes gnomad_r4_STR.allele_sizes.tsv --sr-id RFC1 \\
      --top-trv AoU_I=gnomAD_LR.AoU_I.vep_parsed.annotated.TRV.bed.gz=chr4-39348424-TRV-55 \\
      --out-prefix RFC1/RFC1.hprc_hgsvc
"""
import argparse
import collections
import subprocess

import matplotlib
import pandas as pd

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.patches import Patch  # noqa: E402

plt.rcParams["font.family"] = "Arial"
plt.rcParams["pdf.fonttype"] = 42

INK = "#1d1f24"
INK2 = "#5b5e66"
PALETTE = ["#2a78d6", "#eb6834", "#1baf7a", "#d03b3b", "#8e5bd0", "#eda100", "#17a2b8", "#e377c2", "#8c564b", "#00a6a6",
           "#bcbd22", "#1f3b73"]
# gnomAD genetic ancestry colours (gnomad_methods GEN_ANC_COLORS); ASK (Ashkenazi, HG002) shown as ASJ
ANC_COLORS = {"AFR": "#941494", "AMR": "#ED1E24", "ASJ": "#FF7F50", "EAS": "#108C44", "EUR": "#6AA5CD", "SAS": "#FF9912",
              "OTH": "#ABB9B9"}
ANC_ALIAS = {"ASK": "ASJ"}
OTHER_REP = "#c9c6bf"
NON_REP = "#f1efea"


def read_record(path, trv_id):
    header = subprocess.run(f"gzip -cd {path} | sed -n 1p", shell=True, capture_output=True, text=True).stdout
    cols = header.lstrip("#").rstrip("\n").split("\t")
    chrom, pos = trv_id.split("-")[:2]
    out = subprocess.run(["tabix", path, f"{chrom}:{max(1, int(pos) - 10)}-{int(pos) + 10}"], capture_output=True, text=True,
                         check=True).stdout
    for line in out.splitlines():
        rec = dict(zip(cols, line.split("\t")))
        if rec["ID"] == trv_id:
            return rec
    raise SystemExit(f"{trv_id} not found in {path}")


def canon(kmer):
    return min(kmer[i:] + kmer[:i] for i in range(len(kmer)))


def primitive(km):
    k = len(km)
    return not any(k % d == 0 and km == km[:d] * (k // d) for d in range(1, k))


def discover(seqs, ks, min_frac, max_motifs):
    """Tandem k-mers (k in ks; primitive only) pooled by rotation; keep those covering >= min_frac of tandem bases."""
    tandem = collections.Counter()
    rot_seen = collections.defaultdict(collections.Counter)
    for s in seqs:
        for k in ks:
            for i in range(len(s) - 2 * k + 1):
                km = s[i:i + k]
                if km == s[i + k:i + 2 * k] and primitive(km):
                    c = canon(km)
                    tandem[c] += 1
                    rot_seen[c][km] += 1
    total = sum(tandem.values()) or 1
    keep = [c for c, n in tandem.most_common(max_motifs) if n / total >= min_frac]
    return keep, rot_seen, tandem


def decompose(s, motifs, min_copies):
    """List of (start, end, label) segments; label = canonical motif, 'other_repeat' or 'non_repeat'. At each position the
    kept motif (any rotation, any length) with the longest run of >= min_copies copies wins."""
    rots = {c: {c[i:] + c[:i] for i in range(len(c))} for c in motifs}
    ks = sorted({len(c) for c in motifs})
    segs, i = [], 0
    while i < len(s):
        best, best_len = None, 0
        for c, rs in rots.items():
            k = len(c)
            km = s[i:i + k]
            if km in rs:
                j = i
                while s[j:j + k] == km:
                    j += k
                if (j - i) // k >= min_copies and j - i > best_len:
                    best, best_len = c, j - i
        if best:
            segs.append((i, i + best_len, best))
            i += best_len
            continue
        lab, step = "non_repeat", 1
        for k in ks:
            km = s[i:i + k]
            if len(km) == k and primitive(km) and s[i + k:i + 2 * k] == km:
                lab, step = "other_repeat", 2 * k
                break
        segs.append((i, min(i + step, len(s)), lab))
        i += step
    merged = []
    for a, b, lab in segs:
        if merged and merged[-1][2] == lab and merged[-1][1] == a:
            merged[-1] = (merged[-1][0], b, lab)
        else:
            merged.append((a, b, lab))
    return merged


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--trv", required=True)
    p.add_argument("--trv-id", required=True)
    p.add_argument("--k", type=int, nargs="+",
                   help="Repeat unit lengths to search (default: catalog motif length of the site +/- 1)")
    p.add_argument("--min-frac", type=float, default=0.001)
    p.add_argument("--max-motifs", type=int, default=12)
    p.add_argument("--min-copies", type=int, default=2)
    p.add_argument("--genotypes", help="extract_trv_genotypes.py .genotypes.tsv (ID, sample, hap1_allele, hap2_allele)")
    p.add_argument("--ancestry", help="Sample ancestry table (TSV, optionally gzipped)")
    p.add_argument("--sample-col", default="sample")
    p.add_argument("--ancestry-col", default="ethnicity")
    p.add_argument("--sr-allele-sizes", help="fetch_gnomad_str.py .allele_sizes.tsv (id repunit repunit_count allele_count)")
    p.add_argument("--sr-id", help="gnomAD STR locus id for the SR histogram (e.g. RFC1)")
    p.add_argument("--top-trv", nargs="*", default=[],
                   help="LABEL=path=TRV_ID: other long-read cohorts shown as allele-length histograms on top (e.g. AoU)")
    p.add_argument("--top-colors", nargs="*", default=["#f5b88a", "#9ec5ef", "#9fd8c0"])
    p.add_argument("--label", default="")
    p.add_argument("--out-prefix", required=True)
    args = p.parse_args()

    rec = read_record(args.trv, args.trv_id)
    k0 = min(len(m) for m in rec.get("motifs", rec.get("TRID", "").split("-")[-1]).split(","))
    ks = args.k or [x for x in (k0 - 1, k0, k0 + 1) if x >= 2]
    alts = rec["ALT"].split(",")
    acs = [int(x) if x not in (".", "") else 0 for x in rec["AC"].split(",")]
    alleles = [("REF", rec["REF"][1:].upper(), int(rec["AN"]) - sum(acs))]
    alleles += [(f"ALT{i + 1}", a[1:].upper(), ac) for i, (a, ac) in enumerate(zip(alts, acs)) if ac > 0]
    motifs, rot_seen, tandem = discover([s for _, s, _ in alleles], ks, args.min_frac, args.max_motifs)
    names = {c: c for c in motifs}
    color = {c: PALETTE[i % len(PALETTE)] for i, c in enumerate(motifs)}
    color.update({"other_repeat": OTHER_REP, "non_repeat": NON_REP})

    rows = []
    for name, s, ac in alleles:
        segs = decompose(s, motifs, args.min_copies)
        comp = collections.Counter()
        for a, b, lab in segs:
            comp[lab] += b - a
        dom = max(motifs, key=lambda c: comp[c]) if any(comp[c] for c in motifs) else "non_repeat"
        rows.append({"allele": name, "AC": ac, "length_bp": len(s), "dominant_motif": names.get(dom, dom),
                     **{f"bp_{names[c]}": comp[c] for c in motifs}, "bp_other_repeat": comp["other_repeat"],
                     "bp_non_repeat": comp["non_repeat"], "sequence": s, "_segs": segs, "_dom": dom})
    d = pd.DataFrame(rows)
    order = {c: i for i, c in enumerate(motifs + ["non_repeat"])}
    ref, alt = d.iloc[:1], d.iloc[1:].copy()
    alt["_o"] = alt["_dom"].map(order)
    d = pd.concat([ref, alt.sort_values(["_o", "length_bp"]).drop(columns="_o")], ignore_index=True)
    anc_order = []
    if args.genotypes and args.ancestry:
        anc = pd.read_csv(args.ancestry, sep="\t", dtype=str)
        anc = dict(zip(anc[args.sample_col], anc[args.ancestry_col].str.upper().replace(ANC_ALIAS)))
        g = pd.read_csv(args.genotypes, sep="\t", dtype=str)
        g = g[g["ID"] == args.trv_id]
        carriers = collections.defaultdict(collections.Counter)
        for smp, h1, h2 in zip(g["sample"], g["hap1_allele"], g["hap2_allele"]):
            for a in (h1, h2):
                if isinstance(a, str) and a != ".":
                    carriers["REF" if a == "0" else f"ALT{a}"][anc.get(smp, "OTH")] += 1
        anc_order = [a for a in ANC_COLORS if any(carriers[x][a] for x in carriers)]
        anc_order += sorted({a for c in carriers.values() for a in c} - set(anc_order))
        for a in anc_order:
            d[f"hap_{a}"] = d["allele"].map(lambda x: carriers[x][a])
        mism = d[d[[f"hap_{a}" for a in anc_order]].sum(axis=1) != d["AC"]]
        if len(mism):
            print(f"warning: carrier haplotype count != AC for {len(mism)} alleles: {list(mism['allele'][:5])}")
    d.drop(columns=["_segs", "_dom"]).to_csv(f"{args.out_prefix}.alleles.tsv", sep="\t", index=False)
    mt = pd.DataFrame([{"motif": names[c], "rotations_seen": ",".join(r for r, _ in rot_seen[c].most_common()),
                        "tandem_kmers": tandem[c], "bp_in_alleles": int(d[f"bp_{names[c]}"].sum()),
                        "AC_weighted_bp": int((d[f"bp_{names[c]}"] * d["AC"]).sum()),
                        "n_alleles_with_motif": int((d[f"bp_{names[c]}"] > 0).sum()),
                        "AC_alleles_with_motif": int(d.loc[d[f"bp_{names[c]}"] > 0, "AC"].sum())} for c in motifs])
    mt.to_csv(f"{args.out_prefix}.motifs.tsv", sep="\t", index=False)
    print(mt.drop(columns="rotations_seen").to_string(index=False))
    print(d.groupby("dominant_motif").agg(n_alleles=("allele", "size"), AC=("AC", "sum"),
                                          median_len=("length_bp", "median"), max_len=("length_bp", "max")).to_string())

    n = len(d)
    h = max(4, 0.075 * n + 2.2)
    fs = 5 if n > 30 else 7
    maxlen = d["length_bp"].max()
    # histograms above the structure panel: gnomAD SR (grey) and other long-read cohorts (--top-trv; light colours)
    tops = []
    unit = k0
    if args.sr_allele_sizes and args.sr_id:
        sr = pd.read_csv(args.sr_allele_sizes, sep="\t")
        sr = sr[(sr["id"] == args.sr_id) & (sr["allele_count"] > 0)].copy()
        unit = int(sr["repunit"].str.len().min())
        # repeat-unit count x unit length = allele length in bp (summed over cohorts and repeat-unit variants)
        sr["bp"] = sr["repunit_count"] * sr["repunit"].str.len()
        sr = sr.groupby("bp")["allele_count"].sum()
        tops.append((f"gnomAD short-read allele length ({args.sr_id}; n = {int(sr.sum()):,}; repeat units x {unit} bp)",
                     "SR alleles", sr, "#8a8984"))
    for spec, col in zip(args.top_trv, args.top_colors):
        lab, path, tid = spec.split("=", 2)
        r2 = read_record(path, tid)
        a2 = r2["ALT"].split(",")
        c2 = [int(x) if x not in (".", "") else 0 for x in r2["AC"].split(",")]
        lens = [len(r2["REF"]) - 1] + [len(x) - 1 for x in a2]
        cnts = [int(r2["AN"]) - sum(c2)] + c2
        # binned to the repeat unit (bp, centre of each unit bin)
        t = pd.Series(cnts, index=[round(L / unit) * unit for L in lens])
        t = t[t > 0].groupby(level=0).sum()
        tops.append((f"{lab} long-read allele length ({tid}; n = {int(t.sum()):,}; bins of {unit} bp)", f"{lab} alleles", t,
                     col))
    top_h = 1.4 * len(tops)
    ncol = 2 if anc_order else 1
    fig = plt.figure(figsize=(16 if anc_order else 14, h + top_h))
    gs = fig.add_gridspec(len(tops) + 1, ncol, width_ratios=[20, 80] if anc_order else [1],
                          height_ratios=[1.4] * len(tops) + [h])
    row = len(tops)
    ax = fig.add_subplot(gs[row, ncol - 1])
    axa = fig.add_subplot(gs[row, 0], sharey=ax) if anc_order else None
    axl = axa if anc_order else ax
    for ti, (title, ylab, dist, col) in enumerate(tops):
        axh = fig.add_subplot(gs[ti, ncol - 1], sharex=ax)
        axh.bar(dist.index, dist.to_numpy(), width=unit, color=col, linewidth=0, align="center")
        axh.set_yscale("log")
        axh.set_ylabel(f"{ylab} (log)", fontsize=8, color=INK)
        axh.set_title(title, fontsize=9, color=INK, loc="left")
        axh.tick_params(axis="x", labelbottom=False)
        axh.tick_params(colors=INK2, labelsize=7)
        for sp in ["top", "right"]:
            axh.spines[sp].set_visible(False)
        maxlen = max(maxlen, dist.index.max())
    for y, segs in enumerate(d["_segs"]):
        for a, b, lab in segs:
            ax.barh(y, b - a, left=a, height=0.85, color=color[lab], linewidth=0)
    ax.set_ylim(n - 0.5, -0.5)
    ax.set_xlim(0, maxlen * 1.02)
    ax.set_xlabel("Position in allele (bp, reference + strand)", fontsize=9, color=INK)
    ax.set_title("Repeat-unit structure", fontsize=9.5, color=INK, loc="left")
    handles = [Patch(color=color[c], label=f"{names[c]} ({int(d[f'bp_{names[c]}'].gt(0).sum())} alleles)") for c in motifs]
    handles += [Patch(color=OTHER_REP, label="other repeat"), Patch(color=NON_REP, label="non-repeat")]
    ax.legend(handles=handles, frameon=False, fontsize=7.5, labelcolor=INK, loc="upper center",
              bbox_to_anchor=(0.5, -0.045 * 12 / h), ncol=min(4, len(handles)), title="Repeat unit", title_fontsize=8)
    axes = [ax]
    if anc_order:
        # carrier haplotypes of each allele by genetic ancestry: bar length = AC, each colour = alleles from that ancestry
        maxac = 0
        for y, r in enumerate(d.itertuples(index=False)):
            x0 = 0
            for a in anc_order:
                w = getattr(r, f"hap_{a}")
                if w:
                    axa.barh(y, w, left=x0, height=0.85, color=ANC_COLORS.get(a, "#ABB9B9"), linewidth=0)
                x0 += w
            maxac = max(maxac, x0)
        axa.set_xlim(0, maxac * 1.02)
        axa.set_xlabel("Allele count (carrier haplotypes)", fontsize=9, color=INK)
        axa.set_title("Carrier genetic ancestry", fontsize=9.5, color=INK, loc="left")
        ah = [Patch(color=ANC_COLORS.get(a, "#ABB9B9"), label=f"{a} ({int(d[f'hap_{a}'].sum())} haplotypes)")
              for a in anc_order]
        axa.legend(handles=ah, frameon=False, fontsize=7.5, labelcolor=INK, loc="upper center",
                   bbox_to_anchor=(0.5, -0.045 * 12 / h), ncol=2, title="Genetic ancestry (gnomAD colours)",
                   title_fontsize=8)
        axes.append(axa)
    # allele name, length and AC as row labels on the far left
    axl.set_yticks(range(n))
    axl.set_yticklabels([f"{al}  {ln:,} bp  AC {ac}" for al, ln, ac in zip(d["allele"], d["length_bp"], d["AC"])],
                        fontsize=fs, color=INK)
    if axl is not ax:
        ax.tick_params(labelleft=False)
    for a in axes:
        for sp in ["top", "right", "left"]:
            a.spines[sp].set_visible(False)
        a.tick_params(axis="x", colors=INK2, labelsize=8)
        a.tick_params(axis="y", left=False)
    # genic context of the STR (merged TRV table columns genic_context / genes; omitted when absent)
    ctx = rec.get("genic_context", "")
    where = ""
    if ctx and ctx != ".":
        genes = rec.get("genes", ".")
        if ctx == "intergenic":
            cat = rec.get("intergenic_category", ".")
            where = "intergenic" + (f", {cat} region" if cat not in ("", ".") else "")
        else:
            where = {"coding": "coding exon", "UTR": "UTR", "intronic": "intron"}.get(ctx, ctx)
            where += f" of {genes}" if genes not in ("", ".") else ""
        where = f"; STR in {where}"
    fig.suptitle(f"{args.label + ' ' if args.label else ''}{args.trv_id}{where}: {n} alleles (REF first; ALT grouped by "
                 f"dominant motif, sorted by length; AN {rec['AN']})", fontsize=10.5, color=INK)
    fig.tight_layout(rect=(0, 0, 1, 1 - 0.3 / (h + top_h)))
    # about 0.5 inch between the SR histogram and the allele panels
    fig.subplots_adjust(wspace=0.02, hspace=0.5 / ((h + top_h) / 2))
    fig.savefig(f"{args.out_prefix}.pdf")
    plt.close(fig)


if __name__ == "__main__":
    main()
