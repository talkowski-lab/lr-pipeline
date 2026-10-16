#!/usr/bin/env python3
"""Link methylation outlier haplotypes (RegionHaplotypeMethylation outlier_calls) to rare variants within a window.

Subcommands:
  windows   write a merged BED of [start - window, end + window) around every outlier region on one contig (input for
            `bcftools view -R`).
  link      given the genotypes of those windows (`bcftools query -f '%CHROM\\t%POS\\t%END\\t%ID\\t%FILTER\\t
            %INFO/allele_type\\t%INFO/allele_length\\t%INFO/cadd_phred[\\t%GT:%PS]\\n'` plus a sample list), report
            every rare allele the outlier sample carries within the window, and per outlier call how unusual its rare-
            variant burden is compared with all other retained samples.

Retained samples: in the VCF and the methylation matrix and not QC-excluded; they define allele frequencies and the
burden background. Rare allele: allele count among the retained samples' called haplotypes / allele number <= --max-af
(each ALT allele of a multiallelic record separately).

Outlier haplotype: for HyperASM / HypoASM the haplotype whose region mean is farther from the cohort median of all
retained haplotype means in that region; for HyperMethylated / HypoMethylated both haplotypes. The VCF phase
(GT a|b, PS block) and the methylation hap1/hap2 (read HP tags) may not share one coordinate system, so
alt_on_outlier_haplotype is reported for phased genotypes only and must be validated (e.g. against variants inside the
region) before use; the phase-agnostic columns (carried / zygosity / burden ranks) do not depend on it.

link outputs (<prefix>.):
  outlier_variants.tsv.gz   one row per outlier call x rare allele carried by the outlier sample within the window
  outlier_burden.tsv.gz     one row per outlier call: hap means, outlier haplotype, and for each distance in
                            --distances the number of rare alleles carried (all, heterozygous) and an empirical
                            p-value = (1 + #other retained samples with >= that count) / (1 + #other retained samples)
"""
import argparse
import gzip
import re

import numpy as np
import pandas as pd

GT_RE = re.compile(r"^([0-9]+|\.)(?:([|/])([0-9]+|\.))?$")
VARIANT_COLS = ["region_index", "#chrom", "start", "end", "sample", "label", "outlier_haplotype", "var_pos", "var_end",
                "var_id", "allele", "filter", "allele_type", "allele_length", "cadd_phred", "distance", "in_region", "ac",
                "an", "af", "zygosity", "alt_haplotype", "ps", "alt_on_outlier_haplotype", "n_carrier_samples",
                "n_outlier_samples_in_region", "n_outlier_samples_carrying", "carried_only_by_outliers"]
CALL_COLS = ["region_index", "#chrom", "start", "end", "sample", "label", "hap1_mean", "hap2_mean", "cohort_hap_median",
             "outlier_haplotype"]
INFO_COLS = ["chrom", "pos", "end", "id", "filter", "allele_type", "allele_length", "cadd_phred"]


def read_samples(path):
    with open(path) as fh:
        return [line.split("\t")[0].strip() for line in fh if line.strip()]


def load_calls(path, contig):
    calls = pd.read_csv(path, sep="\t")
    calls = calls[calls["#chrom"] == contig].reset_index(drop=True)
    return calls


def windows(args):
    calls = load_calls(args.outlier_calls, args.contig)
    iv = sorted({(max(0, int(s) - args.window), int(e) + args.window) for s, e in zip(calls["start"], calls["end"])})
    merged = []
    for s, e in iv:
        if merged and s <= merged[-1][1]:
            merged[-1][1] = max(merged[-1][1], e)
        else:
            merged.append([s, e])
    with open(args.out_bed, "w") as out:
        for s, e in merged:
            out.write(f"{args.contig}\t{s}\t{e}\n")
    print(f"{len(calls)} outlier calls on {args.contig}; {len(merged)} merged windows, "
          f"{sum(e - s for s, e in merged):,} bp")


def parse_gt_block(cells):
    """cells: 2-D array of 'GT:PS' strings -> allele0, allele1 (int16, -1 missing), phased (bool), ps (int64, -1)."""
    flat = pd.Series(cells.ravel())
    parts = flat.str.split(":", n=1, expand=True)
    gt = parts[0]
    ps = pd.to_numeric(parts[1], errors="coerce").fillna(-1).astype(np.int64) if parts.shape[1] > 1 \
        else pd.Series(np.full(len(flat), -1, np.int64))
    m = gt.str.extract(GT_RE)
    a0 = pd.to_numeric(m[0], errors="coerce").fillna(-1).astype(np.int16).to_numpy()
    a1 = pd.to_numeric(m[2], errors="coerce").fillna(-1).astype(np.int16).to_numpy()
    phased = (m[1] == "|").to_numpy()
    shape = cells.shape
    return a0.reshape(shape), a1.reshape(shape), phased.reshape(shape), ps.to_numpy().reshape(shape)


def outlier_haplotypes(calls, mean_matrix, retained):
    """Add hap1_mean, hap2_mean, cohort_hap_median, outlier_haplotype to calls."""
    idx = set(calls["region_index"])
    chunks = []
    for chunk in pd.read_csv(mean_matrix, sep="\t", chunksize=20000):
        chunks.append(chunk[chunk["region_index"].isin(idx)])
    m = pd.concat(chunks).set_index("region_index")
    hap_cols = [c for c in m.columns[3:] if c.rsplit("_hap", 1)[0] in retained]
    med = m[hap_cols].median(axis=1)
    h1 = np.array([m.at[r, f"{s}_hap1"] for r, s in zip(calls["region_index"], calls["sample"])], float)
    h2 = np.array([m.at[r, f"{s}_hap2"] for r, s in zip(calls["region_index"], calls["sample"])], float)
    c = med.loc[calls["region_index"]].to_numpy(float)
    calls = calls.copy()
    calls["hap1_mean"], calls["hap2_mean"], calls["cohort_hap_median"] = h1.round(2), h2.round(2), c.round(2)
    asm = calls["label"].str.endswith("ASM").to_numpy()
    pick = np.where(np.abs(h1 - c) >= np.abs(h2 - c), "1", "2")
    calls["outlier_haplotype"] = np.where(asm, pick, "both")
    return calls


def link(args):
    calls = load_calls(args.outlier_calls, args.contig)
    with gzip.open(args.genotypes, "rt") as fh:
        empty = not fh.readline()
    if calls.empty or empty:
        write_outputs(args, [], [])
        print(f"no outlier calls or no genotypes on {args.contig}; wrote empty outputs")
        return
    excluded = set(read_samples(args.excluded_samples)) if args.excluded_samples else set()
    vcf_samples = read_samples(args.vcf_samples)
    header = pd.read_csv(args.mean_matrix, sep="\t", nrows=0).columns[4:]
    meth_samples = {c.rsplit("_hap", 1)[0] for c in header}
    samples = [s for s in vcf_samples if s in meth_samples and s not in excluded]
    col = {s: i for i, s in enumerate(samples)}
    missing = sorted(set(calls["sample"]) - set(col))
    if missing:
        print(f"outlier samples absent from the VCF (dropped): {missing}")
        calls = calls[calls["sample"].isin(col)].reset_index(drop=True)
    calls = outlier_haplotypes(calls, args.mean_matrix, set(samples))
    out_cols = sorted({col[s] for s in calls["sample"]})
    distances = sorted(args.distances)

    info, dosage, code, ps_out = [], [], [], []
    n_records = 0
    reader = pd.read_csv(args.genotypes, sep="\t", header=None, dtype=str, chunksize=args.chunk_size,
                         names=INFO_COLS + vcf_samples, usecols=INFO_COLS + samples)
    for chunk in reader:
        n_records += len(chunk)
        a0, a1, phased, ps = parse_gt_block(chunk[samples].to_numpy())
        max_allele = int(max(a0.max(), a1.max(), 0))
        an = (a0 >= 0).sum(1) + (a1 >= 0).sum(1)
        for k in range(1, max_allele + 1):
            h0, h1 = a0 == k, a1 == k
            ac = h0.sum(1) + h1.sum(1)
            with np.errstate(invalid="ignore", divide="ignore"):
                af = np.where(an > 0, ac / an, np.nan)
            rare = (ac > 0) & (af <= args.max_af)
            if not rare.any():
                continue
            d = (h0[rare].astype(np.int8) + h1[rare].astype(np.int8))
            # 0 none, 1 ALT on first GT allele only (phased), 2 on second only (phased), 3 het unphased, 4 hom / both
            c = np.zeros(d.shape, np.int8)
            ph = phased[rare]
            c[(d == 1) & ph & h0[rare]] = 1
            c[(d == 1) & ph & h1[rare]] = 2
            c[(d == 1) & ~ph] = 3
            c[d == 2] = 4
            sub = chunk.loc[rare, INFO_COLS].copy()
            sub["allele"], sub["ac"], sub["an"], sub["af"] = k, ac[rare], an[rare], np.round(af[rare], 5)
            info.append(sub)
            dosage.append(d)
            code.append(c)
            ps_out.append(ps[rare][:, out_cols])
    if not info:
        write_outputs(args, [], [])
        print("no rare alleles in the windows; wrote empty outputs")
        return
    var = pd.concat(info, ignore_index=True)
    dos = np.vstack(dosage)
    code = np.vstack(code)
    psm = np.vstack(ps_out)
    ps_col = {j: i for i, j in enumerate(out_cols)}
    var["pos"] = var["pos"].astype(np.int64)
    var["end"] = pd.to_numeric(var["end"], errors="coerce").fillna(var["pos"]).astype(np.int64)
    order = np.argsort(var["pos"].to_numpy(), kind="stable")
    var, dos, code, psm = var.iloc[order].reset_index(drop=True), dos[order], code[order], psm[order]
    var["id"] = var["id"].str.slice(0, 120)
    pos = var["pos"].to_numpy()
    vend = var["end"].to_numpy()
    carriers = (dos > 0).sum(1)
    print(f"{n_records:,} records in windows; {len(var):,} rare alleles (AF <= {args.max_af}); "
          f"{len(samples)} retained samples; {len(calls)} outlier calls")

    outl_by_region = calls.groupby("region_index")["sample"].apply(lambda x: [col[s] for s in x]).to_dict()
    max_d = max(distances + [args.window])
    rows, burden = [], []
    for _, r in calls.iterrows():
        s, start, end = col[r["sample"]], int(r["start"]), int(r["end"])
        lo, hi = np.searchsorted(pos, start - max_d - args.max_sv_span), np.searchsorted(pos, end + max_d, side="right")
        p, e = pos[lo:hi], vend[lo:hi]
        dist = np.where(e < start + 1, start + 1 - e, np.where(p > end, p - end, 0))
        d_all = dos[lo:hi]
        het = d_all == 1
        b = {k: r[k] for k in ("region_index", "#chrom", "start", "end", "sample", "label", "hap1_mean", "hap2_mean",
                               "cohort_hap_median", "outlier_haplotype")}
        others = np.ones(len(samples), bool)
        others[s] = False
        for dd in distances:
            w = dist <= dd
            cnt = (d_all[w] > 0).sum(0)
            cnt_het = het[w].sum(0)
            b[f"n_rare_{dd}"] = int(cnt[s])
            b[f"p_rare_{dd}"] = round((1 + (cnt[others] >= cnt[s]).sum()) / (1 + others.sum()), 4)
            b[f"n_rare_het_{dd}"] = int(cnt_het[s])
            b[f"p_rare_het_{dd}"] = round((1 + (cnt_het[others] >= cnt_het[s]).sum()) / (1 + others.sum()), 4)
        burden.append(b)

        sel = np.where((dist <= args.window) & (d_all[:, s] > 0))[0]
        if not len(sel):
            continue
        co = outl_by_region[r["region_index"]]
        for i in sel:
            g = lo + i
            cd = int(code[g, s])
            alt_hap = {1: "1", 2: "2", 3: "unphased", 4: "both"}[cd]
            if r["outlier_haplotype"] == "both" or cd == 4:
                on = "yes"
            elif cd in (1, 2):
                on = "yes" if alt_hap == r["outlier_haplotype"] else "no"
            else:
                on = "NA"
            n_out = int((dos[g, co] > 0).sum())
            rows.append({"region_index": r["region_index"], "#chrom": r["#chrom"], "start": start, "end": end,
                         "sample": r["sample"], "label": r["label"], "outlier_haplotype": r["outlier_haplotype"],
                         "var_pos": int(pos[g]), "var_end": int(vend[g]), "var_id": var.at[g, "id"],
                         "allele": int(var.at[g, "allele"]), "filter": var.at[g, "filter"],
                         "allele_type": var.at[g, "allele_type"], "allele_length": var.at[g, "allele_length"],
                         "cadd_phred": var.at[g, "cadd_phred"], "distance": int(dist[i]), "in_region": dist[i] == 0,
                         "ac": int(var.at[g, "ac"]), "an": int(var.at[g, "an"]), "af": var.at[g, "af"],
                         "zygosity": "hom" if cd == 4 else "het", "alt_haplotype": alt_hap,
                         "ps": int(psm[g, ps_col[s]]) if psm[g, ps_col[s]] >= 0 else "NA",
                         "alt_on_outlier_haplotype": on, "n_carrier_samples": int(carriers[g]),
                         "n_outlier_samples_in_region": len(co), "n_outlier_samples_carrying": n_out,
                         "carried_only_by_outliers": n_out == int(carriers[g])})

    write_outputs(args, rows, burden)
    print(f"{len(rows):,} outlier-call x rare-allele rows")


def burden_cols(distances):
    cols = list(CALL_COLS)
    for dd in sorted(distances):
        cols += [f"n_rare_{dd}", f"p_rare_{dd}", f"n_rare_het_{dd}", f"p_rare_het_{dd}"]
    return cols


def write_outputs(args, rows, burden):
    with gzip.open(f"{args.prefix}.outlier_variants.tsv.gz", "wt") as fh:
        pd.DataFrame(rows, columns=VARIANT_COLS).to_csv(fh, sep="\t", index=False)
    with gzip.open(f"{args.prefix}.outlier_burden.tsv.gz", "wt") as fh:
        pd.DataFrame(burden, columns=burden_cols(args.distances)).to_csv(fh, sep="\t", index=False)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    w = sub.add_parser("windows")
    w.add_argument("--outlier-calls", required=True)
    w.add_argument("--contig", required=True)
    w.add_argument("--window", type=int, default=1000000)
    w.add_argument("--out-bed", required=True)
    w.set_defaults(func=windows)
    k = sub.add_parser("link")
    k.add_argument("--outlier-calls", required=True)
    k.add_argument("--mean-matrix", required=True, help="region x haplotype mean matrix (RegionHaplotypeMethylation)")
    k.add_argument("--genotypes", required=True, help="bcftools query output (see module docstring), gzip ok")
    k.add_argument("--vcf-samples", required=True, help="sample order of the genotype columns (bcftools query -l)")
    k.add_argument("--excluded-samples", help="QC-excluded samples (sample<TAB>anything)")
    k.add_argument("--contig", required=True)
    k.add_argument("--window", type=int, default=1000000)
    k.add_argument("--distances", type=int, nargs="+", default=[10000, 100000, 1000000])
    k.add_argument("--max-af", type=float, default=0.01)
    k.add_argument("--max-sv-span", type=int, default=1000000,
                   help="look this far upstream of the window for long variants whose span reaches it")
    k.add_argument("--chunk-size", type=int, default=20000)
    k.add_argument("--prefix", required=True)
    k.set_defaults(func=link)
    args = ap.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
