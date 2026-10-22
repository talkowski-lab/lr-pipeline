#!/usr/bin/env python3


import argparse

import edlib
from pysam import VariantFile


def parse_intersection_line(line: str):
    fields = line.strip().split("\t")
    palmer_alts_raw = fields[13]
    palmer_alts = palmer_alts_raw.split(",") if "," in palmer_alts_raw else [palmer_alts_raw]
    return {
        "chrom": fields[0],
        "sv_positions": fields[3].split(","),
        "sv_refs": fields[4].split(","),
        "sv_lengths": fields[5].split(","),
        "mei_start_0based": int(fields[7]),
        "mei_pos_1based": int(fields[8]),
        "mei_id": fields[9],
        "mei_length": abs(int(fields[10])),
        "mei_samples": fields[11].rstrip(",").split(","),
        "palmer_ref": fields[12],
        "palmer_alts": palmer_alts,
    }


def compute_size_similarity(target_length: int, mei_length: int) -> float:
    if target_length <= 0 or mei_length <= 0:
        return 0.0
    return min(target_length, mei_length) / max(target_length, mei_length)


def compute_reciprocal_overlap(target_start: int, target_length: int, mei_start: int, mei_length: int) -> float:
    target_end = target_start + target_length
    mei_end = mei_start + mei_length
    overlap = min(target_end, mei_end) - max(target_start, mei_start)
    if overlap <= 0:
        return 0.0
    return min(overlap / target_length, overlap / mei_length)


def compute_sequence_similarity(seq_a: str, seq_b: str) -> float:
    seq_a = seq_a.upper()
    seq_b = seq_b.upper()
    if not seq_a and not seq_b:
        return 1.0
    if not seq_a or not seq_b:
        return 0.0
    if seq_a == seq_b:
        return 1.0
    result = edlib.align(seq_a, seq_b, task="distance")
    max_len = max(len(seq_a), len(seq_b))
    if max_len == 0:
        return 1.0
    return 1.0 - (result["editDistance"] / max_len)


def is_non_ref_genotype(gt: tuple) -> bool:
    if gt is None:
        return False
    return any(allele != 0 for allele in gt)


def record_passes_filters(
    record,
    ref_allele: str,
    sv_length: int,
    mei_data: dict,
    palmer_alts: list[str],
    samples: list[str],
    args: argparse.Namespace,
) -> bool:
    allele_length = record.info.get("allele_length")
    if isinstance(allele_length, (list, tuple)):
        allele_length = abs(allele_length[0])
    else:
        allele_length = abs(allele_length)
    target_length = abs(sv_length)
    if allele_length != target_length - 1:
        return False

    allele_type = record.info.get("allele_type")
    if isinstance(allele_type, (list, tuple)):
        allele_type = allele_type[0]
    if allele_type != "ins":
        return False

    if record.ref != ref_allele:
        return False

    if not record.alts:
        return False

    best_sequence_similarity = 0.0
    for target_alt in record.alts:
        for palmer_alt in palmer_alts:
            sim = compute_sequence_similarity(target_alt, palmer_alt)
            best_sequence_similarity = max(best_sequence_similarity, sim)
    if best_sequence_similarity < args.sequence_similarity:
        return False

    size_similarity = compute_size_similarity(target_length, mei_data["mei_length"])
    if size_similarity < args.size_similarity:
        return False

    reciprocal_overlap = compute_reciprocal_overlap(
        record.pos - 1,
        target_length,
        mei_data["mei_start_0based"],
        mei_data["mei_length"],
    )
    if reciprocal_overlap < args.reciprocal_overlap:
        return False

    breakpoint_distance = abs(record.pos - mei_data["mei_pos_1based"])
    if breakpoint_distance > args.breakpoint_window:
        return False

    sv_samples = [sample for sample in samples if is_non_ref_genotype(record.samples[sample]["GT"])]
    shared = set(sv_samples) & set(mei_data["mei_samples"])
    if len(shared) < args.min_shared_samples:
        return False

    return True


def main() -> None:
    parser = argparse.ArgumentParser(description="Transfer PALMER MEI annotations onto a short-read VCF")
    parser.add_argument("--intersection", required=True, help="Path to bedtools intersect output")
    parser.add_argument("--target-vcf", required=True, help="Target VCF to annotate")
    parser.add_argument("--me-type", required=True, help="MEI class under evaluation")
    parser.add_argument(
        "--output",
        required=True,
        help="Destination TSV containing transferable annotations",
    )
    parser.add_argument(
        "--reciprocal-overlap",
        type=float,
        required=True,
        help="Minimum reciprocal overlap between PALMER and short-read events (0-1)",
    )
    parser.add_argument(
        "--size-similarity",
        type=float,
        required=True,
        help="Minimum size similarity between PALMER and short-read events (0-1)",
    )
    parser.add_argument(
        "--sequence-similarity",
        type=float,
        required=True,
        help="Minimum sequence similarity between PALMER and short-read events (0-1)",
    )
    parser.add_argument(
        "--breakpoint-window",
        type=int,
        required=True,
        help="Maximum distance (bp) between event breakpoints",
    )
    parser.add_argument(
        "--min-shared-samples",
        type=int,
        required=True,
        help="Minimum number of shared samples required for annotation",
    )
    args = parser.parse_args()

    target_vcf = VariantFile(args.target_vcf)
    samples = list(target_vcf.header.samples)

    with open(args.intersection, "r", encoding="utf-8") as source, open(args.output, "w", encoding="utf-8") as output_handle:
        for raw_line in source:
            mei_data = parse_intersection_line(raw_line)
            chrom = mei_data["chrom"]
            for ref, pos, length in zip(mei_data["sv_refs"], mei_data["sv_positions"], mei_data["sv_lengths"]):
                sv_pos = int(pos)
                sv_length = int(length)
                ref_allele = ref.split("_")[0]
                palmer_alts = mei_data["palmer_alts"]

                for record in target_vcf.fetch(chrom, sv_pos - 1, sv_pos):
                    if "ME_TYPE" in record.info:
                        continue

                    if not record_passes_filters(
                        record,
                        ref_allele,
                        sv_length,
                        mei_data,
                        palmer_alts,
                        samples,
                        args,
                    ):
                        continue

                    all_alts = ",".join(record.alts)
                    output_handle.write(f"{chrom}\t{sv_pos}\t{record.ref}\t{all_alts}\t{record.id}\t")
                    output_handle.write(f"{args.me_type}\t{mei_data['mei_length']}\t{mei_data['mei_id']}\n")

    target_vcf.close()


if __name__ == "__main__":
    main()
