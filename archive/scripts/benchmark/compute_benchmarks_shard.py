#!/usr/bin/env python3

import argparse
from typing import Dict, Tuple, List

import pandas as pd
from collections import defaultdict, Counter


def is_af_field(key: str) -> bool:
    lower_key = key.lower()
    return lower_key == "af" or lower_key.startswith("af_") or lower_key.endswith("_af")


def normalize_af_value(value):
    if isinstance(value, tuple) or isinstance(value, list):
        return value[0] if value else None
    return value


def normalize_af_field(field: str) -> frozenset:
    parts = field.lower().replace("af_", "").replace("_af", "").split("_")
    normalized_parts = set()
    for part in parts:
        if part == "male":
            normalized_parts.add("xy")
        elif part == "female":
            normalized_parts.add("xx")
        else:
            normalized_parts.add(part)
    return frozenset(normalized_parts)


def parse_truth_info_string(info_str: str) -> Dict[str, str]:
    items = info_str.strip().split(";") if info_str else []
    kv = {}
    for it in items:
        if "=" in it:
            k, v = it.split("=", 1)
            kv[k] = v
        else:
            kv[it] = ""
    return kv


def parse_vep_header_line(header_line: str) -> Tuple[str, List[str]]:
    line = header_line.strip()
    lower = line.lower()
    if "id=csq" in lower:
        id_key = "CSQ"
    elif "id=vep" in lower:
        id_key = "VEP"
    else:
        raise ValueError("VEP/CSQ ID not found in header line")
    fmt_part = line.split("Format:")[-1].strip().strip('">')
    fmt_fields = [f.strip().lower() for f in fmt_part.split("|")]
    return id_key, fmt_fields


def get_case_insensitive(info_dict: Dict[str, object], key: str):
    target = key.lower()
    for k, v in info_dict.items():
        try:
            if str(k).lower() == target:
                return v
        except Exception:
            continue
    return None


def extract_vep_annotations(
    info_dict: Dict[str, object], vep_key: str, indices: Dict[int, str]
) -> Dict[str, str]:
    vep_string = ""
    v = get_case_insensitive(info_dict, vep_key)
    if v is not None:
        if isinstance(v, tuple) or isinstance(v, list):
            vep_string = v[0] if v else ""
        else:
            vep_string = v
    fields = str(vep_string).lower().split("|") if vep_string else []
    annos = {}
    for idx, cat in indices.items():
        if idx < len(fields) and fields[idx]:
            annos[cat] = fields[idx]
        else:
            annos[cat] = "N/A"
    return annos


def load_enriched_shard(path: str) -> List[Tuple[str, str, str, str]]:
    rows = []
    with open(path, "rt") as f:
        for line in f:
            if not line.strip():
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 4:
                continue
            eval_id, truth_id, eval_info, truth_info = parts[0], parts[1], parts[2], parts[3]
            rows.append((eval_id, truth_id, eval_info, truth_info))
    return rows


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--prefix", required=True)
    ap.add_argument("--contig", required=True)
    ap.add_argument("--matched_shard_tsv", required=True)
    ap.add_argument("--eval_vep_header", required=True)
    ap.add_argument("--truth_vep_header", required=True)
    ap.add_argument("--shard_label", required=True)
    ap.add_argument(
        "--skip_vep_categories", help="Comma-separated list of VEP categories to skip"
    )
    args = ap.parse_args()

    af_out_path = f"{args.prefix}.shard_{args.shard_label}.af_pairs.tsv"
    vep_out_path = f"{args.prefix}.shard_{args.shard_label}.vep_pairs.tsv"

    enriched_rows = load_enriched_shard(args.matched_shard_tsv)

    with open(args.eval_vep_header, "r") as f:
        eval_header_line = f.readline()
    eval_vep_key, eval_vep_fields = parse_vep_header_line(eval_header_line)

    with open(args.truth_vep_header, "r") as f:
        truth_header_line = f.readline()
    truth_vep_key, truth_vep_fields = parse_vep_header_line(truth_header_line)

    # Parse skip_vep_categories
    skip_categories = set()
    if args.skip_vep_categories:
        skip_categories = {
            cat.strip().lower() for cat in args.skip_vep_categories.split(",")
        }
        print(f"Skipping VEP categories: {', '.join(sorted(skip_categories))}")

    # Filter out skipped categories
    common_categories = set(eval_vep_fields) & set(truth_vep_fields)
    common_categories = {cat for cat in common_categories if cat not in skip_categories}

    eval_idx = {
        i: cat for i, cat in enumerate(eval_vep_fields) if cat in common_categories
    }
    truth_idx = {
        i: cat for i, cat in enumerate(truth_vep_fields) if cat in common_categories
    }
    af_rows = []
    vep_counts: Dict[str, Counter] = defaultdict(Counter)

    for _, _, eval_info_str, truth_info_str in enriched_rows:
        eval_info = parse_truth_info_string(eval_info_str)
        truth_info = parse_truth_info_string(truth_info_str)

        eval_af_pairs = {}
        for k, v in eval_info.items():
            if is_af_field(k):
                try:
                    val = v.split(",")[0] if "," in v else v
                    eval_af_pairs[normalize_af_field(k)] = float(val)
                except (ValueError, AttributeError):
                    continue

        truth_af_pairs = {}
        for k, v in truth_info.items():
            if is_af_field(k):
                try:
                    val = v.split(",")[0] if "," in v else v
                    truth_af_pairs[normalize_af_field(k)] = float(val)
                except (ValueError, AttributeError):
                    continue

        for key_set, eval_val in eval_af_pairs.items():
            if key_set in truth_af_pairs and eval_val is not None:
                e = float(eval_val)
                t = float(truth_af_pairs[key_set])
                if e > 0 and t > 0:
                    af_rows.append(
                        {
                            "af_key": "_".join(sorted(list(key_set))),
                            "eval_af": e,
                            "truth_af": t,
                        }
                    )

        # Process VEP annotations
        eval_ann = extract_vep_annotations(eval_info, eval_vep_key, eval_idx)
        truth_ann = extract_vep_annotations(truth_info, truth_vep_key, truth_idx)
        for cat in common_categories:
            eval_val = eval_ann.get(cat, "N/A")
            truth_val = truth_ann.get(cat, "N/A")
            vep_counts[cat][(eval_val, truth_val)] += 1

    if af_rows:
        df_af = pd.DataFrame(af_rows)
        with open(af_out_path, "wt") as f:
            df_af.to_csv(f, sep="\t", index=False)
    else:
        with open(af_out_path, "wt") as f:
            f.write("af_key\teval_af\ttruth_af\n")

    rows = []
    for cat, ctr in vep_counts.items():
        for (e, t), c in ctr.items():
            rows.append({"category": cat, "eval": e, "truth": t, "count": c})
    df_vep = pd.DataFrame(rows, columns=["category", "eval", "truth", "count"])
    with open(vep_out_path, "wt") as f:
        df_vep.to_csv(f, sep="\t", index=False)


if __name__ == "__main__":
    main()
