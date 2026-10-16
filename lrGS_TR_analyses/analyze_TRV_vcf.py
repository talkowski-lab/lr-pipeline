#!/usr/bin/env python3
"""Per-contig tandem-repeat variant (TRV) analysis of gnomAD-LR VCFs. Standard library only.

Input is the TRV `bcftools query` table of one contig (made by run_analyze_TRV_vcf.sh or the AnalyzeTRVariants WDL):
  bcftools view -i 'INFO/allele_type="trv"' <vcf> -Oz -o <prefix>.TRV.vcf.gz
  bcftools query -l <prefix>.TRV.vcf.gz > <prefix>.samples.txt
  bcftools query -f "$QUERY_FORMAT" <prefix>.TRV.vcf.gz | gzip > <prefix>.TRV.query.tsv.gz
with QUERY_FORMAT (tab-separated, one line per site):
  %CHROM %POS %ID %REF %ALT %FILTER %INFO/TRID %INFO/MOTIFS %INFO/AC %INFO/AF %INFO/AN %INFO/MC_allele [%GT per sample]

Subcommands
  analyze          one contig -> <prefix>.TRV.sites.bed.gz, <prefix>.TRV.sample_size_diff.tsv.gz,
                   <prefix>.TRV.sample_motif_count_diff.tsv.gz, <prefix>.TRV.per_sample_summary.tsv
  merge-summaries  sum per-contig per_sample_summary tables -> one genome-wide table

Definitions
  Repeat span      TRID start..end (0-based half-open; min start / max end across TRID components). Allele repeat
                   sequences are REF/ALT with the VCF padding base and any flank outside the TRID span trimmed.
  size diff        len(allele) - len(REF)  (identical for full and trimmed sequences).
  motif count      INFO/MC_allele when present (multi-allelic sites, the pipeline's definition); otherwise a greedy
                   exact-match count of MOTIFS along the repeat sequence (longest motif first). Column mc_source says which.
  size diff / motif length uses the shortest motif.
  Genic context    >=1 bp overlap with CDS -> coding; else exon -> UTR; else transcript -> intronic; else intergenic.
                   Intergenic TRs get the distance (bp gap) to, and gene of, the closest 5'UTR and closest 3'UTR;
                   UTRs are split into 5'/3' per transcript by position relative to its CDS and strand.
  Per-sample cells "d1,d2" follow GT allele order (GT is mostly unphased); "." = missing allele; haploid calls give one value.
  Per-sample summary counts PASS sites only unless --all-filters.

Usage
  analyze_TRV_vcf.py analyze --query-tsv q.tsv.gz --samples samples.txt --gtf gencode.gtf.gz --prefix out/chr22 [--all-filters]
  analyze_TRV_vcf.py merge-summaries --summaries a.tsv b.tsv ... --out TRV.per_sample_summary.tsv
"""
import argparse
import bisect
import collections
import gzip
import re

CONTEXTS = ["coding", "UTR", "intronic", "intergenic"]
SUMMARY_FIELDS = ["n_sites_called", "n_nonref_sites", "n_nonref_alleles", "n_expansion_alleles", "n_contraction_alleles",
                  "n_same_length_nonref_alleles", "n_motif_count_changed_alleles", "sum_abs_size_diff", "sum_size_diff"]
BED_HEADER = ["#chrom", "start", "end", "ID", "FILTER", "TRID", "motifs", "ref_seq", "ref_len", "ref_motif_count", "mc_source",
              "n_alt", "alt_len", "alt_size_diff", "alt_size_diff_per_motif", "alt_motif_count", "alt_motif_count_diff",
              "AC", "AF", "AN", "genic_context", "genes", "dist_closest_5UTR", "gene_closest_5UTR", "dist_closest_3UTR",
              "gene_closest_3UTR"]


class IntervalIndex:
    """Per-chromosome sorted intervals (0-based half-open) with a label; overlap and nearest queries."""

    def __init__(self, rows):
        by_chrom = collections.defaultdict(list)
        for chrom, s, e, label in rows:
            by_chrom[chrom].append((s, e, label))
        self.data = {}
        for chrom, iv in by_chrom.items():
            iv.sort()
            starts = [x[0] for x in iv]
            maxlen = max(x[1] - x[0] for x in iv)
            by_end = sorted((x[1], x[2]) for x in iv)
            self.data[chrom] = (iv, starts, maxlen, [x[0] for x in by_end], [x[1] for x in by_end])

    def overlaps(self, chrom, s, e):
        if chrom not in self.data:
            return set()
        iv, starts, maxlen, _, _ = self.data[chrom]
        out = set()
        j = bisect.bisect_left(starts, e) - 1
        while j >= 0 and starts[j] > s - maxlen:
            if iv[j][1] > s:
                out.add(iv[j][2])
            j -= 1
        return out

    def nearest(self, chrom, s, e):
        """Gap (bp) to the closest non-overlapping interval and its label(s); (None, '.') if chrom absent."""
        if chrom not in self.data:
            return None, "."
        iv, starts, _, ends, end_labels = self.data[chrom]
        best, labels = None, set()
        k = bisect.bisect_left(starts, e)
        if k < len(iv):
            best, labels = iv[k][0] - e, {x[2] for x in iv[k:] if x[0] == iv[k][0]}
        k = bisect.bisect_right(ends, s) - 1
        if k >= 0:
            d = s - ends[k]
            same = {end_labels[i] for i in range(k, -1, -1) if ends[i] == ends[k]}
            if best is None or d < best:
                best, labels = d, same
            elif d == best:
                labels |= same
        return best, ",".join(sorted(labels)) if labels else "."


def load_gtf(gtf):
    feats = {"CDS": [], "exon": [], "transcript": [], "UTR": []}
    cds_span, strand = {}, {}
    with gzip.open(gtf, "rt") as f:
        for line in f:
            if line.startswith("#"):
                continue
            t = line.rstrip("\n").split("\t")
            if t[2] not in feats:
                continue
            gene = re.search(r'gene_name "([^"]+)"', t[8]).group(1)
            tx = re.search(r'transcript_id "([^"]+)"', t[8]).group(1)
            s, e = int(t[3]) - 1, int(t[4])
            feats[t[2]].append((t[0], s, e, gene, tx))
            strand[tx] = t[6]
            if t[2] == "CDS":
                lo, hi = cds_span.get(tx, (s, e))
                cds_span[tx] = (min(lo, s), max(hi, e))
    utr5, utr3 = [], []
    for chrom, s, e, gene, tx in feats["UTR"]:
        lo, hi = cds_span[tx]
        upstream = e <= lo if strand[tx] == "+" else s >= hi
        (utr5 if upstream else utr3).append((chrom, s, e, gene))
    idx = {k: IntervalIndex([(c, s, e, g) for c, s, e, g, _ in v]) for k, v in feats.items() if k != "UTR"}
    idx["UTR5"], idx["UTR3"] = IntervalIndex(utr5), IntervalIndex(utr3)
    return idx


def genic_annotation(idx, chrom, s, e):
    for ctx, feat in [("coding", "CDS"), ("UTR", "exon"), ("intronic", "transcript")]:
        genes = idx[feat].overlaps(chrom, s, e)
        if genes:
            return ctx, ",".join(sorted(genes)), ".", ".", ".", "."
    d5, g5 = idx["UTR5"].nearest(chrom, s, e)
    d3, g3 = idx["UTR3"].nearest(chrom, s, e)
    return "intergenic", ".", "." if d5 is None else str(d5), g5, "." if d3 is None else str(d3), g3


def greedy_motif_count(seq, motifs):
    i = n = 0
    while i < len(seq):
        for m in motifs:
            if seq.startswith(m, i):
                n += 1
                i += len(m)
                break
        else:
            i += 1
    return n


def analyze(args):
    idx = load_gtf(args.gtf)
    samples = [x.strip() for x in open(args.samples) if x.strip()]
    summary = collections.defaultdict(collections.Counter)

    bed = gzip.open(f"{args.prefix}.TRV.sites.bed.gz", "wt")
    size_mat = gzip.open(f"{args.prefix}.TRV.sample_size_diff.tsv.gz", "wt")
    mc_mat = gzip.open(f"{args.prefix}.TRV.sample_motif_count_diff.tsv.gz", "wt")
    bed.write("\t".join(BED_HEADER) + "\n")
    size_mat.write("\t".join(["ID"] + samples) + "\n")
    mc_mat.write("\t".join(["ID"] + samples) + "\n")

    n = 0
    with gzip.open(args.query_tsv, "rt") as q:
        for line in q:
            t = line.rstrip("\n").split("\t")
            chrom, pos, vid, ref, alt, filt, trid, motifs, ac, af, an, mca = t[:12]
            gts = t[12:]
            pos0 = int(pos) - 1
            comps = [c.split("-") for c in trid.split(",")]
            s0 = min(int(c[1]) for c in comps)
            e0 = max(int(c[2]) for c in comps)
            left = max(s0 - pos0, 0)
            right = max(pos0 + len(ref) - e0, 0)
            alleles = [ref] + alt.split(",")
            rep = [a[left:len(a) - right] for a in alleles]
            mot = sorted(set(motifs.split(",")), key=lambda m: (-len(m), m))
            min_motif = min(len(m) for m in mot)
            if mca != ".":
                mc = [int(x) for x in mca.split(",")]
                mc_source = "INFO_MC_allele"
            else:
                mc = [greedy_motif_count(r, mot) for r in rep]
                mc_source = "greedy_exact"
            size_d = [len(a) - len(ref) for a in alleles]
            mc_d = [m - mc[0] for m in mc]
            ctx, genes, d5, g5, d3, g3 = genic_annotation(idx, chrom, s0, e0)

            bed.write("\t".join([chrom, str(s0), str(e0), vid, filt, trid, ",".join(mot), rep[0], str(len(rep[0])), str(mc[0]),
                                 mc_source, str(len(alleles) - 1), ",".join(str(len(r)) for r in rep[1:]),
                                 ",".join(map(str, size_d[1:])), ",".join(f"{d / min_motif:.2f}" for d in size_d[1:]),
                                 ",".join(map(str, mc[1:])), ",".join(map(str, mc_d[1:])), ac, af, an,
                                 ctx, genes, d5, g5, d3, g3]) + "\n")

            count = args.all_filters or filt == "PASS"
            srow, mrow = [vid], [vid]
            for smp, gt in zip(samples, gts):
                idxs = re.split(r"[|/]", gt)
                srow.append(",".join("." if i == "." else str(size_d[int(i)]) for i in idxs))
                mrow.append(",".join("." if i == "." else str(mc_d[int(i)]) for i in idxs))
                if not count:
                    continue
                called = [int(i) for i in idxs if i != "."]
                if not called:
                    continue
                c = summary[(smp, ctx)]
                c["n_sites_called"] += 1
                nonref = [i for i in called if i != 0]
                if nonref:
                    c["n_nonref_sites"] += 1
                for i in nonref:
                    c["n_nonref_alleles"] += 1
                    d = size_d[i]
                    if d > 0:
                        c["n_expansion_alleles"] += 1
                    elif d < 0:
                        c["n_contraction_alleles"] += 1
                    else:
                        c["n_same_length_nonref_alleles"] += 1
                    c["n_motif_count_changed_alleles"] += mc_d[i] != 0
                    c["sum_abs_size_diff"] += abs(d)
                    c["sum_size_diff"] += d
            size_mat.write("\t".join(srow) + "\n")
            mc_mat.write("\t".join(mrow) + "\n")
            n += 1
    for fh in (bed, size_mat, mc_mat):
        fh.close()
    write_summary(f"{args.prefix}.TRV.per_sample_summary.tsv", samples, summary)
    print(f"{args.prefix}: {n:,} TRV sites, {len(samples)} samples")


def write_summary(path, samples, summary):
    with open(path, "w") as f:
        f.write("\t".join(["sample", "genic_context"] + SUMMARY_FIELDS) + "\n")
        for smp in samples:
            for ctx in CONTEXTS + ["all"]:
                if ctx == "all":
                    c = sum((summary[(smp, x)] for x in CONTEXTS), collections.Counter())
                else:
                    c = summary[(smp, ctx)]
                f.write("\t".join([smp, ctx] + [str(c[k]) for k in SUMMARY_FIELDS]) + "\n")


def merge_summaries(args):
    samples, summary = [], collections.defaultdict(collections.Counter)
    for path in args.summaries:
        with open(path) as f:
            header = f.readline().rstrip("\n").split("\t")
            for line in f:
                row = dict(zip(header, line.rstrip("\n").split("\t")))
                if row["genic_context"] == "all":
                    continue
                if row["sample"] not in samples:
                    samples.append(row["sample"])
                c = summary[(row["sample"], row["genic_context"])]
                for k in SUMMARY_FIELDS:
                    c[k] += int(row[k])
    write_summary(args.out, samples, summary)


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = p.add_subparsers(dest="cmd", required=True)
    a = sub.add_parser("analyze", help="Analyze one contig's TRV query table")
    a.add_argument("--query-tsv", required=True, help="Gzipped bcftools query table (format in module docstring)")
    a.add_argument("--samples", required=True, help="Sample list in VCF order (bcftools query -l)")
    a.add_argument("--gtf", required=True, help="GENCODE GTF (gz) with CDS/exon/transcript/UTR features and gene_name")
    a.add_argument("--prefix", required=True, help="Output path prefix")
    a.add_argument("--all-filters", action="store_true", help="Per-sample summary over all sites (default: PASS only)")
    m = sub.add_parser("merge-summaries", help="Sum per-contig per-sample summary tables")
    m.add_argument("--summaries", nargs="+", required=True)
    m.add_argument("--out", required=True)
    args = p.parse_args()
    analyze(args) if args.cmd == "analyze" else merge_summaries(args)


if __name__ == "__main__":
    main()
