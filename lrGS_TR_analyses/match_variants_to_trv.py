#!/usr/bin/env python3
"""Match external indel calls inside a TR locus to the TRV record and in-repeat routine indels of gnomAD-LR VCFs.

Each external variant (whitespace-separated: chrom pos id ref alt; VCF-style, 1-based, anchored) is applied to the TRV REF
sequence to build the full haplotype of the repeat region. Matching is done on haplotype sequence, so indels placed at
different positions inside the repeat still match when they produce the same sequence:
  routine_match      PASS-or-not routine del/ins record in the VCF (TRID = the TRV) giving the identical haplotype
  trv_exact_allele   a TRV ALT allele whose sequence equals that haplotype (allele index, AC)
  same_length        TRV ALT alleles (AC > 0) with the same length change: count, summed AC, and the closest one by
                     approximate differing bases (difflib; len(longer) - matched bases)
Usage:
  match_variants_to_trv.py --variants variants.tsv --trv-id chr6-30949309-TRV-2636 \
      --vcf hprc_hgsvc=hprc.vcf.gz AoU_I=aou.vcf.gz --out matches.tsv
"""
import argparse
import difflib
import subprocess

import pandas as pd


def query(vcf, expr, fmt):
    out = subprocess.run(["bcftools", "query", "-i", expr, "-f", fmt, vcf], check=True, capture_output=True, text=True).stdout
    return [line.split("\t") for line in out.strip().split("\n") if line]


def apply(ref_seq, ref_start, pos, ref, alt):
    off = pos - ref_start
    if off < 0 or ref_seq[off:off + len(ref)] != ref:
        return None
    return ref_seq[:off] + alt + ref_seq[off + len(ref):]


def diff_bases(a, b):
    m = difflib.SequenceMatcher(None, a, b, autojunk=False)
    return max(len(a), len(b)) - sum(bl.size for bl in m.get_matching_blocks())


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--variants", required=True)
    p.add_argument("--trv-id", required=True)
    p.add_argument("--vcf", nargs="+", required=True, help="LABEL=VCF")
    p.add_argument("--out", required=True)
    args = p.parse_args()

    ext = pd.read_csv(args.variants, sep=r"\s+", header=None, names=["chrom", "pos", "id", "ref", "alt"], dtype={"pos": int})
    rows = {i: {"ext_id": r.id, "pos": r.pos, "delta": len(r.alt) - len(r.ref)} for i, r in ext.iterrows()}
    for item in args.vcf:
        label, vcf = item.split("=", 1)
        trv = query(vcf, f'ID="{args.trv_id}"', "%POS\t%REF\t%ALT\t%INFO/AC\n")[0]
        tpos, tref, talts = int(trv[0]), trv[1], trv[2].split(",")
        tac = [int(x) for x in trv[3].split(",")]
        tend = tpos + len(tref) - 1
        routine = query(vcf, 'INFO/allele_type="del" || INFO/allele_type="ins"',
                        "%POS\t%ID\t%REF\t%ALT\t%FILTER\t%INFO/AC\n")
        routine = [r for r in routine if tpos <= int(r[0]) and int(r[0]) + len(r[2]) - 1 <= tend]
        routine_hap = {}
        for r in routine:
            h = apply(tref, tpos, int(r[0]), r[2], r[3])
            if h is not None:
                routine_hap.setdefault(h, []).append(f"{r[1]}(FILTER={r[4]},AC={r[5]})")
        allele_idx = {a: i + 1 for i, a in enumerate(talts)}
        for i, r in ext.iterrows():
            h = apply(tref, tpos, r.pos, r.ref, r.alt)
            row = rows[i]
            if h is None:
                row[f"{label}_note"] = "REF does not match TRV reference"
                continue
            row[f"{label}_routine_match"] = ";".join(routine_hap.get(h, [])) or "."
            k = allele_idx.get(h)
            row[f"{label}_trv_exact_allele"] = f"ALT{k}(AC={tac[k - 1]})" if k else "."
            same = [(j, a, tac[j]) for j, a in enumerate(talts) if len(a) == len(h) and tac[j] > 0]
            row[f"{label}_same_length_alleles"] = len(same)
            row[f"{label}_same_length_AC"] = sum(c for _, _, c in same)
            if same:
                best = min(same, key=lambda x: diff_bases(x[1], h))
                row[f"{label}_closest_same_length"] = f"ALT{best[0] + 1}(AC={best[2]},~{diff_bases(best[1], h)} bp diff)"
            else:
                row[f"{label}_closest_same_length"] = "."
        print(f"{label}: TRV {args.trv_id} REF {len(tref) - 1} bp, {len(talts)} ALT alleles; "
              f"routine del/ins inside: {len(routine)}")
    out = pd.DataFrame(rows.values())
    out.to_csv(args.out, sep="\t", index=False)
    with pd.option_context("display.width", 250, "display.max_columns", 30, "display.max_colwidth", 60):
        print(out.to_string(index=False))


if __name__ == "__main__":
    main()
