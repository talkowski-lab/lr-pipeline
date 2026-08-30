#!/usr/bin/env python

import pysam
import argparse
import csv
from collections import defaultdict
from intervaltree import IntervalTree
from Levenshtein import distance as levenshtein_distance

VAMOS_REF_SAMPLE = "GRCh38_to_GRCh38/GRCh38_to_GRCh38"


def load_vcf_records(vcf_path, contigs=None):
    with pysam.VariantFile(vcf_path) as vcf:
        samples = list(vcf.header.samples)
        records = {}
        tree = defaultdict(IntervalTree)
        iterator = (vcf.fetch(c) for c in contigs) if contigs else [vcf.fetch()]
        for chunk in iterator:
            for record in chunk:
                chrom, start, end = (
                    record.chrom,
                    record.pos,
                    max(record.stop, record.pos + 1),
                )
                key = (chrom, start, end)
                records[key] = record
                tree[chrom][start : end + 1] = key
    return records, tree, samples


def parse_vamos_gt(record, sample_name):
    sample_idx = list(record.samples.keys()).index(sample_name)
    gt_str = str(record).split("\t")[9 + sample_idx].split(":")[0].replace("\n", "")
    return tuple(int(x) if x != "." else None for x in gt_str.split("/"))


def build_vamos_sequence(record, sample_name, ref_fa=None):
    motifs = record.info.get("RU")
    motifs = motifs if isinstance(motifs, tuple) else motifs.split(",")
    alt_haps = record.info.get("ALTANNO")
    alt_haps = alt_haps if isinstance(alt_haps, tuple) else alt_haps.split(",")
    sequences = []
    for allele_idx in parse_vamos_gt(record, sample_name):
        if allele_idx is None:
            continue
        try:
            motif_idx = [int(x) for x in alt_haps[allele_idx - 1].split("-")]
            sequences.append("".join(motifs[i] for i in motif_idx).upper())
        except BaseException:
            pass
    return tuple(sequences)


def build_trgt_sequence(record, sample_name, ref_fa=None):
    sequences = []
    for allele_idx in record.samples[sample_name]["GT"]:
        if allele_idx is None:
            continue
        seq = str(record.ref) if allele_idx == 0 else str(record.alts[allele_idx - 1])
        sequences.append(seq.upper())
    return tuple(sequences)


def get_vamos_motif_size(record):
    motifs = record.info.get("RU")
    motifs = motifs if isinstance(motifs, tuple) else motifs.split(",")
    return min(len(m) for m in motifs)


def get_trgt_motif_size(record):
    motifs = record.info.get("MOTIFS")
    motifs = motifs if isinstance(motifs, tuple) else motifs.split(",")
    return min(len(m) for m in motifs)


def merge_sequences(
    keys, records, sample, window_start, window_end, chrom, ref_fa, build_seq_func
):
    sorted_keys = sorted(keys, key=lambda k: k[1])
    hap_seqs = {0: [], 1: []}
    last_pos = {0: window_start - 1, 1: window_start - 1}
    for key in sorted_keys:
        record = records[key]
        variant_seqs = build_seq_func(record, sample, ref_fa)
        var_start_orig, var_end = record.pos, record.stop
        for hap_idx, var_seq in enumerate(variant_seqs[:2]):
            var_start = var_start_orig
            if var_start <= last_pos[hap_idx]:
                overlap_len = last_pos[hap_idx] - var_start + 1
                var_seq = var_seq[overlap_len:]
                var_start = last_pos[hap_idx] + 1
            if last_pos[hap_idx] < var_start - 1:
                gap_seq = ref_fa.fetch(
                    region=f"{chrom}:{last_pos[hap_idx]+1}-{var_start-1}"
                ).upper()
                hap_seqs[hap_idx].append(gap_seq)
            if var_seq:
                hap_seqs[hap_idx].append(var_seq)
                last_pos[hap_idx] = var_end
    result = []
    for hap_idx in [0, 1]:
        if last_pos[hap_idx] < window_end:
            gap_seq = ref_fa.fetch(
                region=f"{chrom}:{last_pos[hap_idx]+1}-{window_end}"
            ).upper()
            hap_seqs[hap_idx].append(gap_seq)
        result.append("".join(hap_seqs[hap_idx]))
    return tuple(result)


def sequence_similarity(seq1, seq2):
    if not seq1 and not seq2:
        return 1.0
    if not seq1 or not seq2:
        return 0.0
    return 1.0 - levenshtein_distance(seq1, seq2) / max(len(seq1), len(seq2))


def check_allele_match(vamos_seqs, trgt_seqs, threshold=0.9):
    if len(vamos_seqs) != 2 or len(trgt_seqs) != 2:
        return "DIVERGENT", 0.0, 0
    sim_direct = (
        sequence_similarity(vamos_seqs[0], trgt_seqs[0])
        + sequence_similarity(vamos_seqs[1], trgt_seqs[1])
    ) / 2
    edit_direct = levenshtein_distance(
        vamos_seqs[0], trgt_seqs[0]
    ) + levenshtein_distance(vamos_seqs[1], trgt_seqs[1])
    sim_swap = (
        sequence_similarity(vamos_seqs[0], trgt_seqs[1])
        + sequence_similarity(vamos_seqs[1], trgt_seqs[0])
    ) / 2
    edit_swap = levenshtein_distance(
        vamos_seqs[0], trgt_seqs[1]
    ) + levenshtein_distance(vamos_seqs[1], trgt_seqs[0])
    if sim_swap > sim_direct:
        similarity, edit_dist = sim_swap, edit_swap
    else:
        similarity, edit_dist = sim_direct, edit_direct
    if similarity == 1.0:
        category = "EXACT"
    elif similarity >= threshold:
        category = "SIMILAR"
    else:
        category = "DIVERGENT"
    return category, similarity, edit_dist


def find_overlaps(vamos_records, vamos_tree, trgt_records, trgt_tree):
    overlap_map = {"vamos": defaultdict(list), "trgt": defaultdict(list)}
    overlapped = {"vamos": set(), "trgt": set()}
    for caller, records, other_tree, other_name in [
        ("vamos", vamos_records, trgt_tree, "trgt"),
        ("trgt", trgt_records, vamos_tree, "vamos"),
    ]:
        for key in records:
            chrom, pos, end = key
            for interval in other_tree[chrom][pos : end + 1]:
                overlapped[caller].add(key)
                overlapped[other_name].add(interval.data)
                overlap_map[caller][key].append(interval.data)
    all_keys = overlapped["vamos"] | overlapped["trgt"]
    parent = {k: k for k in all_keys}

    def find(x):
        if parent[x] != x:
            parent[x] = find(parent[x])
        return parent[x]

    def union(x, y):
        px, py = find(x), find(y)
        if px != py:
            parent[px] = py

    for caller in ["vamos", "trgt"]:
        for key, refs in overlap_map[caller].items():
            for ref_key in refs:
                union(key, ref_key)
    components = defaultdict(set)
    for k in all_keys:
        components[find(k)].add(k)
    regions = []
    for comp_keys in components.values():
        chrom = next(iter(comp_keys))[0]
        vamos_keys = tuple(sorted(k for k in comp_keys if k in vamos_records))
        trgt_keys = tuple(sorted(k for k in comp_keys if k in trgt_records))
        if vamos_keys and trgt_keys:
            regions.append((chrom, vamos_keys, trgt_keys))
    return regions, overlapped


def get_non_ref_counts(
    vamos_keys, trgt_keys, vamos_records, trgt_records, trgt_sample, vamos_target_sample
):
    trgt_h0_nonref = False
    trgt_h1_nonref = False
    for k in trgt_keys:
        rec = trgt_records[k]
        gt = rec.samples[trgt_sample]["GT"]
        if len(gt) > 0 and gt[0] is not None and gt[0] > 0:
            trgt_h0_nonref = True
        if len(gt) > 1 and gt[1] is not None and gt[1] > 0:
            trgt_h1_nonref = True
    trgt_count = int(trgt_h0_nonref) + int(trgt_h1_nonref)

    vamos_h0_nonref = False
    vamos_h1_nonref = False
    for k in vamos_keys:
        rec = vamos_records[k]
        tgt_gt = parse_vamos_gt(rec, vamos_target_sample)
        ref_gt = parse_vamos_gt(rec, VAMOS_REF_SAMPLE)
        t0 = tgt_gt[0] if len(tgt_gt) > 0 else None
        r0 = ref_gt[0] if len(ref_gt) > 0 else None
        if t0 is not None and t0 != r0:
            vamos_h0_nonref = True
        t1 = tgt_gt[1] if len(tgt_gt) > 1 else None
        r1 = ref_gt[1] if len(ref_gt) > 1 else None
        if t1 is not None and t1 != r1:
            vamos_h1_nonref = True
    vamos_count = int(vamos_h0_nonref) + int(vamos_h1_nonref)
    return trgt_count, vamos_count


def compare_region(
    chrom,
    vamos_keys,
    trgt_keys,
    vamos_records,
    trgt_records,
    vamos_sample,
    trgt_sample,
    ref_fa,
):
    all_starts = [vamos_records[k].pos for k in vamos_keys] + [
        trgt_records[k].pos for k in trgt_keys
    ]
    all_ends = [vamos_records[k].stop for k in vamos_keys] + [
        trgt_records[k].stop for k in trgt_keys
    ]
    window_start, window_end = min(all_starts), max(all_ends)
    vamos_merged = merge_sequences(
        vamos_keys,
        vamos_records,
        vamos_sample,
        window_start,
        window_end,
        chrom,
        ref_fa,
        build_vamos_sequence,
    )
    trgt_merged = merge_sequences(
        trgt_keys,
        trgt_records,
        trgt_sample,
        window_start,
        window_end,
        chrom,
        ref_fa,
        build_trgt_sequence,
    )
    category, similarity, edit_dist = check_allele_match(vamos_merged, trgt_merged)
    ref_seq = ref_fa.fetch(region=f"{chrom}:{window_start}-{window_end}").upper()
    vamos_ref_dist = sum(levenshtein_distance(seq, ref_seq) for seq in vamos_merged)
    trgt_ref_dist = sum(levenshtein_distance(seq, ref_seq) for seq in trgt_merged)
    length_diff = sum(len(s) for s in vamos_merged) - sum(len(s) for s in trgt_merged)
    locus_size = window_end - window_start + 1
    vamos_motif_size = min(get_vamos_motif_size(vamos_records[k]) for k in vamos_keys)
    trgt_motif_size = min(get_trgt_motif_size(trgt_records[k]) for k in trgt_keys)
    trgt_nr_count, vamos_nr_count = get_non_ref_counts(
        vamos_keys, trgt_keys, vamos_records, trgt_records, trgt_sample, vamos_sample
    )
    return {
        "vamos_keys": vamos_keys,
        "trgt_keys": trgt_keys,
        "vamos_seqs": vamos_merged,
        "trgt_seqs": trgt_merged,
        "window_start": window_start,
        "window_end": window_end,
        "category": category,
        "similarity": similarity,
        "edit_dist": edit_dist,
        "vamos_ref_dist": vamos_ref_dist,
        "trgt_ref_dist": trgt_ref_dist,
        "length_diff": length_diff,
        "locus_size": locus_size,
        "vamos_motif_size": vamos_motif_size,
        "trgt_motif_size": trgt_motif_size,
        "trgt_non_ref_count": trgt_nr_count,
        "vamos_non_ref_count": vamos_nr_count,
    }


def format_key(key):
    return f"{key[0]}:{key[1]}-{key[2]}"


def print_summary(vamos_records, trgt_records, overlapped, results):
    print("\n--- Site-Level ---")
    for name, records, ovl in [
        ("Vamos", vamos_records, overlapped["vamos"]),
        ("TRGT", trgt_records, overlapped["trgt"]),
    ]:
        total = len(records)
        print(f"{name} Overlapped: {len(ovl)} / {total} ({len(ovl)/total:.1%})")
        print(
            f"{name} Exclusive: {total - len(ovl)} / {total} ({(total-len(ovl))/total:.1%})\n"
        )
    print("--- Genotype Level ---")
    total = len(results)
    by_cat = defaultdict(list)
    for r in results:
        by_cat[r["category"]].append(r)
    print(f"Total regions: {total}")
    for cat in ["EXACT", "SIMILAR", "DIVERGENT"]:
        print(f" {cat}: {len(by_cat[cat])} ({len(by_cat[cat])/total:.1%})")


def write_stats_tsv(output_path, results, sample_id):
    with open(output_path, "w", newline="") as f:
        writer = csv.writer(f, delimiter="\t")
        writer.writerow(
            [
                "sample_id",
                "category",
                "similarity",
                "edit_dist",
                "vamos_ref_dist",
                "trgt_ref_dist",
                "length_diff",
                "locus_size",
                "vamos_motif_size",
                "trgt_motif_size",
                "trgt_non_ref_count",
                "vamos_non_ref_count",
            ]
        )
        for r in results:
            writer.writerow(
                [
                    sample_id,
                    r["category"],
                    r["similarity"],
                    r["edit_dist"],
                    r["vamos_ref_dist"],
                    r["trgt_ref_dist"],
                    r["length_diff"],
                    r["locus_size"],
                    r["vamos_motif_size"],
                    r["trgt_motif_size"],
                    r["trgt_non_ref_count"],
                    r["vamos_non_ref_count"],
                ]
            )


def run_benchmark(
    vamos_vcf, trgt_vcf, ref_fa_path, contigs=None, output_stats=None, sample_id=None
):
    vamos_records, vamos_tree, vamos_samples = load_vcf_records(vamos_vcf, contigs)
    trgt_records, trgt_tree, trgt_samples = load_vcf_records(trgt_vcf, contigs)
    ref_fa = pysam.FastaFile(ref_fa_path)

    trgt_sample = trgt_samples[0]

    regions, overlapped = find_overlaps(
        vamos_records, vamos_tree, trgt_records, trgt_tree
    )
    results = [
        compare_region(
            chrom,
            vamos_keys,
            trgt_keys,
            vamos_records,
            trgt_records,
            trgt_sample,
            trgt_sample,
            ref_fa,
        )
        for chrom, vamos_keys, trgt_keys in regions
    ]
    print_summary(vamos_records, trgt_records, overlapped, results)

    if output_stats:
        write_stats_tsv(output_stats, results, sample_id or "unknown")


def main():
    parser = argparse.ArgumentParser(description="Benchmark VAMOS and TRGT VCF files")
    parser.add_argument("--vamos-vcf", required=True)
    parser.add_argument("--trgt-vcf", required=True)
    parser.add_argument("--ref-fa", required=True)
    parser.add_argument("--contigs", help="Comma-separated list of contigs")
    parser.add_argument(
        "--output-stats", default=None, help="Output TSV file for per-sample statistics"
    )
    parser.add_argument(
        "--sample-id", default=None, help="Sample ID for statistics output"
    )
    args = parser.parse_args()
    contigs = [c.strip() for c in args.contigs.split(",")] if args.contigs else None
    run_benchmark(
        args.vamos_vcf,
        args.trgt_vcf,
        args.ref_fa,
        contigs,
        args.output_stats,
        args.sample_id,
    )


if __name__ == "__main__":
    main()
