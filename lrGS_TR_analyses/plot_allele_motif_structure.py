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
  .alleles.tsv  per allele: name, AC, length, bp per motif, dominant motif, sequence
  .motifs.tsv   per motif: rotations seen, tandem bp in alleles, AC-weighted bp, n alleles containing it

Usage:
  plot_allele_motif_structure.py --trv gnomAD_LR.hprc_hgsvc.TRV.bed.gz --trv-id chr4-39348424-TRV-55 --label RFC1 \\
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
    h = max(4, 0.075 * n + 1.8)
    fig, ax = plt.subplots(figsize=(13, h))
    for y, segs in enumerate(d["_segs"]):
        for a, b, lab in segs:
            ax.barh(y, b - a, left=a, height=0.85, color=color[lab], linewidth=0)
    maxlen = d["length_bp"].max()
    for y, (al, ln, ac) in enumerate(zip(d["allele"], d["length_bp"], d["AC"])):
        ax.text(maxlen * 1.01, y, f"{al}  {ln:,} bp  AC {ac}", va="center", fontsize=5 if n > 60 else 7,
                color=INK)
    ax.set_ylim(n - 0.5, -0.5)
    ax.set_xlim(0, maxlen * 1.18)
    ax.set_yticks([])
    ax.set_xlabel("Position in allele (bp, reference + strand)", fontsize=9, color=INK)
    for sp in ["top", "right", "left"]:
        ax.spines[sp].set_visible(False)
    ax.tick_params(colors=INK2, labelsize=8)
    handles = [Patch(color=color[c], label=f"{names[c]} ({int(d[f'bp_{names[c]}'].gt(0).sum())} alleles)") for c in motifs]
    handles += [Patch(color=OTHER_REP, label="other repeat"), Patch(color=NON_REP, label="non-repeat")]
    ax.legend(handles=handles, frameon=False, fontsize=7.5, labelcolor=INK, loc="upper center",
              bbox_to_anchor=(0.45, -0.04 * 12 / h), ncol=min(5, len(handles)))
    ax.set_title(f"{args.label + ' ' if args.label else ''}{args.trv_id}: repeat-unit structure of {n} alleles (REF first; ALT "
                 f"grouped by dominant motif, sorted by length; AN {rec['AN']})", fontsize=10, color=INK)
    fig.tight_layout()
    fig.savefig(f"{args.out_prefix}.pdf")
    plt.close(fig)


if __name__ == "__main__":
    main()
