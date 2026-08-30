#!/usr/bin/env python3

import argparse
import pysam
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from scipy.stats import pearsonr
import os
import tarfile
import subprocess
from collections import defaultdict
from collections import Counter


def is_af_field(key):
    lower_key = key.lower()
    return lower_key == "af" or lower_key.startswith("af_") or lower_key.endswith("_af")


def normalize_af_value(value):
    return value[0] if isinstance(value, tuple) else value


def normalize_af_field(field):
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


def plot_af_correlation(df, pop, output_dir, log_scale):
    plt.figure(figsize=(10, 10))
    sns.scatterplot(data=df, x="eval_af", y="truth_af")

    if log_scale:
        plt.xscale("log")
        plt.yscale("log")
    plt.xlabel(f"Eval VCF AF{' (log scale)' if log_scale else ''}")
    plt.ylabel(f"Truth VCF AF{' (log scale)' if log_scale else ''}")
    plt.title(f"{pop}")
    plt.grid(True, which="both", ls="--")

    r, _ = pearsonr(df["eval_af"], df["truth_af"])
    r_squared = r**2
    plt.text(
        0.05,
        0.95,
        f"$R^2 = {r_squared:.4f}$",
        transform=plt.gca().transAxes,
        fontsize=12,
        verticalalignment="top",
    )

    plot_path = os.path.join(output_dir, f"{pop}.png")
    plt.savefig(plot_path)
    plt.close()


def get_vep_format(vcf_path):
    with pysam.VariantFile(vcf_path) as vcf:
        for record in vcf.header.records:
            if record.key == "INFO" and record.get("ID"):
                vep_id = record.get("ID")
                if vep_id.lower() in ["vep", "csq"]:
                    description_str = record.get("Description", "")
                    format_str = description_str.split("Format: ")[-1]
                    cleaned_format_str = format_str.strip('"').strip().lower()
                    return vep_id, cleaned_format_str


def get_vep_annotations(info_dict, vep_key, vep_idx):
    vep_string = ""
    if vep_key in info_dict:
        vep_value = info_dict.get(vep_key, "")
        if isinstance(vep_value, tuple):
            vep_string = vep_value[0] if vep_value else ""
        else:
            vep_string = vep_value

    vep_fields = vep_string.lower().split("|")
    vep_categories = {
        category: (
            vep_fields[index]
            if index < len(vep_fields) and vep_fields[index]
            else "N/A"
        )
        for index, category in vep_idx.items()
    }
    return vep_categories


def plot_vep_heatmap(eval_values, truth_values, column, output_dir):
    # Skip if no data or only N/A
    eval_value_counts = Counter(eval_values)
    if not eval_value_counts:
        print(f"No data for VEP category {column}, skipping plot")
        return
    elif len(eval_value_counts) == 1 and "N/A" in eval_value_counts:
        print(
            f"Skipping plot for VEP category '{column}' as it only "
            f"contains 'N/A' values."
        )
        return

    # Determine labels
    top_n = 10
    if len(eval_value_counts) > top_n:
        top_labels = [val for val, count in eval_value_counts.most_common(top_n)]
        col_labels = top_labels + ["Other"]
        row_labels = top_labels + ["Other"]
        print(
            f"Warning: Too many distinct values "
            f"({len(eval_value_counts)}) for '{column}'. "
            f"Plotting top {top_n} + Other."
        )
    else:
        all_unique_values = set(eval_values) | set(truth_values)
        sorted_labels = sorted(
            list(all_unique_values),
            key=lambda x: eval_value_counts.get(x, 0),
            reverse=True,
        )
        col_labels = sorted_labels
        row_labels = sorted_labels

    # Create concordance matrix
    concordance = pd.DataFrame(0, index=row_labels, columns=col_labels)
    for eval_val, truth_val in zip(eval_values, truth_values):
        col = eval_val if eval_val in col_labels and eval_val != "Other" else "Other"
        row = truth_val if truth_val in row_labels and truth_val != "Other" else "Other"
        if col in concordance.columns and row in concordance.index:
            concordance.loc[row, col] += 1

    # Create annotation and percentage matrices
    annot_matrix = pd.DataFrame("", index=row_labels, columns=col_labels)
    percent_matrix = pd.DataFrame(0.0, index=row_labels, columns=col_labels)
    zero_count_mask = pd.DataFrame(False, index=row_labels, columns=col_labels)
    column_totals = concordance.sum(axis=0)
    for col in col_labels:
        total = column_totals[col]
        if total > 0:
            for row in row_labels:
                count = concordance.loc[row, col]
                percentage = (count / total) * 100
                percent_matrix.loc[row, col] = percentage
                annot_matrix.loc[row, col] = f"{count}/{total}\n({percentage:.1f}%)"
                if count == 0:
                    zero_count_mask.loc[row, col] = True
        else:
            for row in row_labels:
                annot_matrix.loc[row, col] = "0/0\n(0.0%)"
                zero_count_mask.loc[row, col] = True

    # Create plot
    fig_height = max(10, len(row_labels) * 0.6)
    fig_width = max(12, len(col_labels) * 0.8)
    plt.figure(figsize=(fig_width, fig_height))

    masked_percent_matrix = np.ma.masked_where(
        zero_count_mask.values, percent_matrix.values
    )

    sns.heatmap(
        masked_percent_matrix,
        annot=annot_matrix,
        fmt="s",
        cmap="viridis",
        cbar=True,
        cbar_kws={"label": "Concordance (%)"},
        linewidths=0.5,
        linecolor="black",
    )

    ax = plt.gca()

    # Add white cells for zero counts
    for i, row in enumerate(row_labels):
        for j, col in enumerate(col_labels):
            if zero_count_mask.loc[row, col]:
                ax.add_patch(
                    plt.Rectangle(
                        (j, i),
                        1,
                        1,
                        fill=True,
                        facecolor="white",
                        linewidth=0.5,
                        edgecolor="black",
                    )
                )
                ax.text(
                    j + 0.5,
                    i + 0.5,
                    annot_matrix.loc[row, col],
                    ha="center",
                    va="center",
                    fontsize=10,
                )

    plt.title(f"{column}")
    plt.xlabel("Eval VCF Annotation")
    plt.ylabel("Truth VCF Annotation")
    plt.xticks(rotation=45, ha="right")
    plt.yticks(rotation=0)
    plt.tight_layout(pad=2.0)

    plot_path = os.path.join(output_dir, f"{column.replace('/', '_')}.png")
    plt.savefig(plot_path, bbox_inches="tight")
    plt.close()


def write_vep_table(eval_values, truth_values, column, output_dir):
    eval_labels = sorted(list(set(val for val in eval_values)))
    truth_labels = sorted(list(set(val for val in truth_values)))
    if not eval_labels:
        return

    concordance = pd.DataFrame(0, index=truth_labels, columns=eval_labels)
    for eval_val, truth_val in zip(eval_values, truth_values):
        if eval_val in eval_labels and truth_val in truth_labels:
            concordance.loc[truth_val, eval_val] += 1

    table_path = os.path.join(output_dir, f"{column}.tsv")
    concordance.to_csv(table_path, sep="\t")


def write_summary_stats(final_vcf_path, contig, output_path):
    total_variants = 0
    match_counts = defaultdict(int)

    with pysam.VariantFile(final_vcf_path) as vcf_in:
        for record in vcf_in:
            total_variants += 1
            if "gnomAD_V4_match" in record.info:
                match_type = record.info["gnomAD_V4_match"]
                match_counts[match_type] += 1
            else:
                match_counts["UNMATCHED"] += 1

    total_matched = total_variants - match_counts["UNMATCHED"]
    summary_stats = {
        "contig": contig,
        "total_variants": total_variants,
        "total_matched": total_matched,
        "total_unmatched": match_counts["UNMATCHED"],
        "percent_matched": (
            (total_matched / total_variants * 100) if total_variants > 0 else 0.0
        ),
        "percent_unmatched": (
            (match_counts["UNMATCHED"] / total_variants * 100)
            if total_variants > 0
            else 0.0
        ),
    }

    for match_type, count in match_counts.items():
        if match_type != "UNMATCHED":
            summary_stats[f"{match_type.lower()}_count"] = count
            summary_stats[f"{match_type.lower()}_percent"] = (
                (count / total_variants * 100) if total_variants > 0 else 0.0
            )

    with open(output_path, "w") as f:
        f.write("\t".join(summary_stats.keys()) + "\n")
        f.write("\t".join(str(v) for v in summary_stats.values()) + "\n")


def write_summary_table(final_vcf_path, truth_variants, vep_keys, output_path):
    all_variants_data = []

    eval_vep_key, truth_vep_key, eval_idx, truth_idx = vep_keys
    common_vep_categories = set(eval_idx.values())

    with pysam.VariantFile(final_vcf_path) as vcf_in:
        for record in vcf_in:
            row_data = {
                "eval_variant_id": record.id,
                "match_status": False,
                "truth_variant_id": ".",
            }

            if "gnomAD_V4_match_ID" in record.info:
                match_id = record.info["gnomAD_V4_match_ID"]
                if match_id in truth_variants:
                    row_data["match_status"] = True
                    row_data["truth_variant_id"] = match_id

                    truth_info = truth_variants[match_id]

                    eval_af_pairs = {
                        normalize_af_field(k): normalize_af_value(v)
                        for k, v in record.info.items()
                        if is_af_field(k)
                    }
                    truth_af_pairs = {
                        normalize_af_field(k): normalize_af_value(v)
                        for k, v in truth_info.items()
                        if is_af_field(k)
                    }

                    common_af_keys = set(eval_af_pairs.keys()) & set(
                        truth_af_pairs.keys()
                    )
                    for af_key_set in common_af_keys:
                        af_key_str = "_".join(sorted(list(af_key_set)))
                        row_data[f"{af_key_str}_eval"] = eval_af_pairs[af_key_set]
                        row_data[f"{af_key_str}_truth"] = truth_af_pairs[af_key_set]

                    eval_annos = get_vep_annotations(
                        record.info, eval_vep_key, eval_idx
                    )
                    truth_annos = get_vep_annotations(
                        truth_info, truth_vep_key, truth_idx
                    )
                    for category in common_vep_categories:
                        if category in eval_annos or category in truth_annos:
                            row_data[f"{category}_eval"] = eval_annos.get(
                                category, "N/A"
                            )
                            row_data[f"{category}_truth"] = truth_annos.get(
                                category, "N/A"
                            )

            all_variants_data.append(row_data)

    df = pd.DataFrame(all_variants_data)
    cols_to_drop = []
    prefixes = set()
    for col in df.columns:
        if col.endswith("_eval") or col.endswith("_truth"):
            prefixes.add(col.rsplit("_", 1)[0])

    # Drop columns that are empty in both eval and truth
    for prefix in prefixes:
        eval_col = f"{prefix}_eval"
        truth_col = f"{prefix}_truth"
        if eval_col in df.columns and truth_col in df.columns:
            eval_is_empty = (df[eval_col].isna() | (df[eval_col] == "N/A")).all()
            truth_is_empty = (df[truth_col].isna() | (df[truth_col] == "N/A")).all()
            if eval_is_empty and truth_is_empty:
                cols_to_drop.extend([eval_col, truth_col])
    df.drop(columns=cols_to_drop, inplace=True, errors="ignore")

    df.to_csv(output_path, sep="\t", index=False, na_rep=".")


def main():
    parser = argparse.ArgumentParser(
        description="Annotate VCF and/or run benchmarking."
    )
    parser.add_argument("--prefix", help="Prefix for output files.")
    parser.add_argument("--contig", help="Contig being processed.")
    parser.add_argument("--exact_matched_vcf", help="VCF of variants matched exactly.")
    parser.add_argument(
        "--truvari_matched_vcf", help="VCF of variants matched by Truvari."
    )
    parser.add_argument(
        "--truvari_too_small_vcf", help="VCF of variants too small for Truvari."
    )
    parser.add_argument(
        "--truvari_unmatched_vcf", help="Unmatched VCF from Truvari step."
    )
    parser.add_argument(
        "--closest_bed", help="BED file with closest matches from bedtools."
    )
    parser.add_argument("--vcf_truth_snv", help="Original truth VCF with SNVs.")
    parser.add_argument("--vcf_truth_sv", help="Original truth VCF with SVs.")
    parser.add_argument(
        "--log_scale", action="store_true", help="Use log scale for AF plots."
    )
    args = parser.parse_args()

    bedtools_matches = defaultdict(list)
    with open(args.closest_bed, "r") as f:
        for line in f:
            if line.startswith("query_svid"):
                continue
            fields = line.strip().split("\t")
            if len(fields) >= 2:
                query_svid = fields[0]
                ref_svid = fields[1]
                if ref_svid != ".":
                    bedtools_matches[query_svid].append(ref_svid)

    # --- Part 1: Annotate Bedtools Matches  ---
    bedtools_matched_tmp_path = f"{args.prefix}.bedtools_matched.tmp.vcf"
    final_unmatched_tmp_path = f"{args.prefix}.final_unmatched.tmp.vcf"
    bedtools_matched_out_path = f"{args.prefix}.bedtools_matched.vcf.gz"
    final_unmatched_out_path = f"{args.prefix}.final_unmatched.vcf.gz"

    vcf_in_unmatched = pysam.VariantFile(args.truvari_unmatched_vcf)
    if "gnomAD_V4_match" not in vcf_in_unmatched.header.info:
        vcf_in_unmatched.header.info.add(
            "gnomAD_V4_match", "1", "String", "Matching status against gnomAD v4."
        )
    if "gnomAD_V4_match_ID" not in vcf_in_unmatched.header.info:
        vcf_in_unmatched.header.info.add(
            "gnomAD_V4_match_ID", "1", "String", "Matching variant ID from gnomAD v4."
        )

    with open(bedtools_matched_tmp_path, "w") as matched_out, open(
        final_unmatched_tmp_path, "w"
    ) as unmatched_out:
        matched_out.write(str(vcf_in_unmatched.header))
        unmatched_out.write(str(vcf_in_unmatched.header))
        for record in vcf_in_unmatched:
            if record.id in bedtools_matches:
                record.info["gnomAD_V4_match"] = "BEDTOOLS_CLOSEST"
                record.info["gnomAD_V4_match_ID"] = bedtools_matches[record.id][0]
                matched_out.write(str(record))
            else:
                unmatched_out.write(str(record))

    subprocess.run(
        [
            "bcftools",
            "view",
            "-Oz",
            "-o",
            bedtools_matched_out_path,
            bedtools_matched_tmp_path,
        ],
        check=True,
    )
    subprocess.run(["tabix", "-p", "vcf", "-f", bedtools_matched_out_path], check=True)

    subprocess.run(
        [
            "bcftools",
            "view",
            "-Oz",
            "-o",
            final_unmatched_out_path,
            final_unmatched_tmp_path,
        ],
        check=True,
    )
    subprocess.run(["tabix", "-p", "vcf", "-f", final_unmatched_out_path], check=True)

    # --- Part 2: Concatenate VCFs ---
    vcf_files_to_concat = [
        args.exact_matched_vcf,
        args.truvari_matched_vcf,
        args.truvari_too_small_vcf,
        bedtools_matched_out_path,
        final_unmatched_out_path,
    ]
    final_vcf_path = f"{args.prefix}.final_annotated.vcf.gz"
    concat_cmd = [
        "bcftools",
        "concat",
        "-a",
        "-Oz",
        "-o",
        final_vcf_path,
    ] + vcf_files_to_concat
    subprocess.run(concat_cmd, check=True)
    subprocess.run(["tabix", "-p", "vcf", "-f", final_vcf_path], check=True)

    # --- Part 3: Benchmark Concatenated VCF ---
    truth_variants = {}
    for vcf_path in [args.vcf_truth_snv, args.vcf_truth_sv]:
        with pysam.VariantFile(vcf_path) as vcf_in:
            for record in vcf_in:
                truth_variants[record.id] = record.info

    matched_data = []
    with pysam.VariantFile(final_vcf_path) as vcf_in:
        for record in vcf_in:
            if "gnomAD_V4_match" in record.info:
                # if record.info['gnomAD_V4_match'] == 'BEDTOOLS_CLOSEST':
                #     continue

                match_id = record.info["gnomAD_V4_match_ID"]
                if match_id in truth_variants:
                    matched_data.append(
                        {"eval": record.info, "truth": truth_variants[match_id]}
                    )
                else:
                    print(
                        f"Matching variant {match_id} for {record.id} "
                        "not found in truth VCF"
                    )

    # Output directories
    output_basedir = f"{args.prefix}_benchmark_results"
    af_plot_dir = os.path.join(output_basedir, "AF_plots", args.contig)
    vep_plot_dir = os.path.join(output_basedir, "VEP_plots", args.contig)
    vep_table_dir = os.path.join(output_basedir, "VEP_tables", args.contig)
    os.makedirs(af_plot_dir, exist_ok=True)
    os.makedirs(vep_plot_dir, exist_ok=True)
    os.makedirs(vep_table_dir, exist_ok=True)

    # Data structures
    eval_vep_key, eval_vep_format = get_vep_format(final_vcf_path)
    truth_vep_key, truth_vep_format = get_vep_format(args.vcf_truth_snv)
    eval_vep_categories = eval_vep_format.split("|")
    truth_vep_categories = truth_vep_format.split("|")

    common_categories = set(eval_vep_categories) & set(truth_vep_categories)
    eval_idx = {
        i: cat for i, cat in enumerate(eval_vep_categories) if cat in common_categories
    }
    truth_idx = {
        i: cat for i, cat in enumerate(truth_vep_categories) if cat in common_categories
    }

    af_data = defaultdict(list)
    vep_data_by_category = defaultdict(lambda: {"eval": [], "truth": []})

    # Plotting
    for item in matched_data:
        # AF data
        eval_af_pairs = {
            normalize_af_field(k): normalize_af_value(v)
            for k, v in item["eval"].items()
            if is_af_field(k)
        }
        truth_af_pairs = {
            normalize_af_field(k): normalize_af_value(v)
            for k, v in item["truth"].items()
            if is_af_field(k)
        }
        for eval_key, eval_af in eval_af_pairs.items():
            if eval_key in truth_af_pairs:
                if eval_af > 0 and truth_af_pairs[eval_key] > 0:
                    af_data[eval_key].append(
                        {"eval_af": eval_af, "truth_af": truth_af_pairs[eval_key]}
                    )

        # VEP data
        eval_annos = get_vep_annotations(item["eval"], eval_vep_key, eval_idx)
        truth_annos = get_vep_annotations(item["truth"], truth_vep_key, truth_idx)
        for category in common_categories:
            if category in eval_annos or category in truth_annos:
                vep_data_by_category[category]["eval"].append(
                    eval_annos.get(category, "N/A")
                )
                vep_data_by_category[category]["truth"].append(
                    truth_annos.get(category, "N/A")
                )

    # AF plots
    for af_key, af_data in af_data.items():
        af_key_str = "_".join(sorted(list(af_key)))
        plot_af_correlation(
            pd.DataFrame(af_data), af_key_str, af_plot_dir, args.log_scale
        )

    # VEP plots
    for category, data in vep_data_by_category.items():
        plot_vep_heatmap(data["eval"], data["truth"], category, vep_plot_dir)
        write_vep_table(data["eval"], data["truth"], category, vep_table_dir)

    # Summary table
    summary_table_path = f"{args.prefix}.benchmark_summary.tsv"
    vep_keys = (eval_vep_key, truth_vep_key, eval_idx, truth_idx)
    write_summary_table(final_vcf_path, truth_variants, vep_keys, summary_table_path)

    # Summary stats
    summary_stats_path = f"{args.prefix}.summary_stats.tsv"
    write_summary_stats(final_vcf_path, args.contig, summary_stats_path)

    # Tarball
    with tarfile.open(f"{args.prefix}.benchmarks.tar.gz", "w:gz") as tar:
        tar.add(output_basedir, arcname=os.path.basename(output_basedir))


if __name__ == "__main__":
    main()
