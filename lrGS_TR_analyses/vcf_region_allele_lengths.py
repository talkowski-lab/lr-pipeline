#!/usr/bin/env python3
"""Per-sample allele lengths of a tandem-repeat region reconstructed from the variants a VCF calls inside it.

1. Extract records lying inside --region (1-based, inclusive; an indel's anchor base may sit one base before the start):
   POS >= start - 1 and record end <= end (record end = INFO/END for symbolic alleles when defined, else POS + len(REF) - 1).
   TRV summary records (ID containing "TRV", or INFO/allele_type == trv when that field exists) are skipped unless
   --keep-trv-records, so they are not double counted with the routine indels.
2. Each ALT allele's type and size are inferred from REF/ALT (no INFO/allele_type needed):
     snv  len(REF) == len(ALT) == 1        mnv  len(REF) == len(ALT) > 1         (size 0)
     ins  len(ALT) > len(REF)              del  len(REF) > len(ALT)              (size = len(ALT) - len(REF))
     "*"  spanning deletion, size 0        <SYMBOLIC>  size from INFO/SVLEN if defined (sign from <DEL>/<INS>/<DUP>),
                                           otherwise size 0 and counted in hapN_n_unsized
3. For each sample and haplotype (GT allele order), sum the sizes of the ALT alleles carried. Output one row per sample:
     hapN_length       region reference length + summed length change (bp)
     hapN_units        hapN_length / --motif-len (estimated repeat-unit count)
     hapN_sv_size      hapN_length - region reference length (bp; SV size relative to REF)
     hapN_variants     ALT alleles carried on that haplotype as ID:type+size (e.g. 6_58705:del-57; comma-separated)
     hapN_n_snv / hapN_n_ins / hapN_n_del   counts of carried ALT alleles by inferred type (mnv counted with snv)
     hapN_n_unsized    symbolic alleles with no SVLEN (size unknown, counted as 0)
     hapN_n_missing    records with a missing call on that haplotype
   plus phased (yes: all records with a call are phased; no: at least one unphased GT, so the haplotype split is GT order
   only) and n_records_in_region.
Requires bcftools. Uses -r (indexed VCF) when an index exists, otherwise -t (streams the file). Optional INFO fields
(allele_type, END, SVLEN) are read only when the header defines them; "chr6" vs "6" naming is resolved from the header;
bcftools errors are reported.

Usage:
  vcf_region_allele_lengths.py --vcf input.vcf.gz --region chr6:30949310-30951945 --motif-len 57 --out alleles.tsv
"""
import argparse
import os
import re
import subprocess

import pandas as pd


def bcftools(args):
    r = subprocess.run(["bcftools"] + args, capture_output=True, text=True)
    if r.returncode != 0:
        raise SystemExit(f"bcftools {' '.join(args[:3])} ... failed (exit {r.returncode}):\n{r.stderr.strip()}")
    return r.stdout


def header_info(vcf):
    """INFO IDs and contig names declared in the VCF header."""
    header = bcftools(["view", "-h", vcf])
    info = set(re.findall(r"^##INFO=<ID=([^,>]+)", header, flags=re.M))
    contigs = set(re.findall(r"^##contig=<ID=([^,>]+)", header, flags=re.M))
    return info, contigs


def classify(ref, alt, svlen):
    """Inferred type and size (bp, ALT minus REF) of one ALT allele; size None when it cannot be determined."""
    if alt == "*":
        return "span_del", 0
    if alt.startswith("<"):
        kind = alt.strip("<>").split(":")[0].upper()
        if svlen is None:
            return kind.lower(), None
        size = abs(svlen)
        return kind.lower(), -size if kind == "DEL" else size
    if len(ref) == len(alt):
        return ("snv" if len(ref) == 1 else "mnv"), 0
    return ("ins" if len(alt) > len(ref) else "del"), len(alt) - len(ref)


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--vcf", required=True)
    p.add_argument("--region", required=True, help="chrom:start-end, 1-based inclusive")
    p.add_argument("--motif-len", type=float, required=True, help="Repeat unit length (bp) for the unit-count estimate")
    p.add_argument("--pass-only", action="store_true", help="Use only FILTER PASS records")
    p.add_argument("--keep-trv-records", action="store_true")
    p.add_argument("--out", required=True)
    args = p.parse_args()

    chrom, rng = args.region.split(":")
    start, end = (int(x.replace(",", "")) for x in rng.split("-"))
    ref_len = end - start + 1
    indexed = any(os.path.exists(args.vcf + ext) for ext in (".tbi", ".csi"))
    info, contigs = header_info(args.vcf)
    if contigs and chrom not in contigs:  # chr6 vs 6 naming
        alt_name = chrom[3:] if chrom.startswith("chr") else f"chr{chrom}"
        if alt_name in contigs:
            print(f"note: VCF names the contig {alt_name}, not {chrom}; using {alt_name}")
            chrom = alt_name
        else:
            raise SystemExit(f"contig {chrom} not in VCF header (e.g. {sorted(contigs)[:5]})")
    q_region = f"{chrom}:{max(start - 1, 1)}-{end}"
    samples = bcftools(["query", "-l", args.vcf]).split()
    opt = {k: (f"%INFO/{k}" if k in info else ".") for k in ("allele_type", "END", "SVLEN")}
    fmt = f"%CHROM\t%POS\t%ID\t%REF\t%ALT\t%FILTER\t{opt['allele_type']}\t{opt['END']}\t{opt['SVLEN']}[\t%GT]\n"
    raw = bcftools(["query", "-r" if indexed else "-t", q_region, "-f", fmt, args.vcf])

    records = []
    for line in raw.strip().split("\n"):
        if not line:
            continue
        t = line.split("\t")
        pos, vid, ref, filt, atype, rec_end, svlens = int(t[1]), t[2], t[3], t[5], t[6], t[7], t[8]
        alts, gts = t[4].split(","), t[9:]
        stop = int(rec_end) if rec_end not in (".", "") and any(a.startswith("<") for a in alts) else pos + len(ref) - 1
        if t[0] != chrom or pos < start - 1 or stop > end:
            continue
        if not args.keep_trv_records and ("TRV" in vid or atype == "trv"):
            continue
        if args.pass_only and filt != "PASS":
            continue
        sv = [None if x in (".", "") else int(float(x)) for x in svlens.split(",")]
        sv = sv + [sv[-1]] * (len(alts) - len(sv)) if sv else [None] * len(alts)
        if vid == ".":
            vid = f"{t[0]}-{pos}-{ref[:10]}-{alts[0][:10]}"
        records.append((vid, [classify(ref, a, sv[i]) for i, a in enumerate(alts)], gts))
    kinds = pd.Series([k for _, alleles, _ in records for k, _ in alleles]).value_counts().to_dict()
    print(f"{args.region}: reference length {ref_len} bp; records used {len(records)}; samples {len(samples)}; "
          f"{'indexed (-r)' if indexed else 'streamed (-t)'}; inferred ALT allele types {kinds}")

    rows = []
    for si, sample in enumerate(samples):
        hap = {h: {"delta": 0, "ids": [], "miss": 0, "snv": 0, "ins": 0, "del": 0, "unsized": 0} for h in (1, 2)}
        phased, ploidy = True, 2
        for vid, alleles, gts in records:
            gt = gts[si]
            idx = re.split(r"[|/]", gt)
            ploidy = len(idx)
            if "/" in gt and len(idx) > 1:
                phased = False
            for h, a in enumerate(idx[:2], start=1):
                if a == ".":
                    hap[h]["miss"] += 1
                elif a != "0":
                    kind, size = alleles[int(a) - 1]
                    if size is None:
                        hap[h]["unsized"] += 1
                        size = 0
                    hap[h]["delta"] += size
                    if kind in ("snv", "mnv"):
                        hap[h]["snv"] += 1
                    elif kind in ("ins", "del"):
                        hap[h][kind] += 1
                    elif kind in ("dup", "cnv") or (kind not in ("span_del",) and size > 0):
                        hap[h]["ins"] += 1
                    elif size < 0:
                        hap[h]["del"] += 1
                    hap[h]["ids"].append(f"{vid}:{kind}{size:+d}" if size else f"{vid}:{kind}")
        row = {"sample": sample, "phased": "yes" if phased else "no", "n_records_in_region": len(records)}
        for h in (1, 2):
            if h > ploidy:
                for k in ("length", "units", "sv_size", "variants", "n_snv", "n_ins", "n_del", "n_unsized", "n_missing"):
                    row[f"hap{h}_{k}"] = "."
                continue
            length = ref_len + hap[h]["delta"]
            row[f"hap{h}_length"] = length
            row[f"hap{h}_units"] = round(length / args.motif_len, 2)
            row[f"hap{h}_sv_size"] = hap[h]["delta"]
            row[f"hap{h}_variants"] = ",".join(hap[h]["ids"]) or "."
            for k in ("snv", "ins", "del", "unsized"):
                row[f"hap{h}_n_{k}"] = hap[h][k]
            row[f"hap{h}_n_missing"] = hap[h]["miss"]
        rows.append(row)
    cols = ["sample", "hap1_length", "hap2_length", "hap1_units", "hap2_units", "hap1_sv_size", "hap2_sv_size",
            "hap1_variants", "hap2_variants", "hap1_n_snv", "hap2_n_snv", "hap1_n_ins", "hap2_n_ins", "hap1_n_del", "hap2_n_del",
            "hap1_n_unsized", "hap2_n_unsized", "hap1_n_missing", "hap2_n_missing", "phased", "n_records_in_region"]
    out = pd.DataFrame(rows)[cols]
    out.to_csv(args.out, sep="\t", index=False)
    sv = pd.concat([out["hap1_sv_size"], out["hap2_sv_size"]])
    sv = pd.to_numeric(sv, errors="coerce").dropna().astype(int)
    print("haplotype SV size distribution (bp):")
    print(sv.value_counts().sort_index().to_string())


if __name__ == "__main__":
    main()
