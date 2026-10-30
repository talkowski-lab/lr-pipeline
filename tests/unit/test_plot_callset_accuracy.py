import subprocess
import sys

import pandas as pd
from conftest import REPO_ROOT

SCRIPT = "scripts/qc/plot_callset_accuracy.py"
KEYS = ["level", "sample", "allele_type", "length_bin", "region", "match_stage"]


def accuracy_tsv(path, counts):
    """A CountCallsetAccuracy TSV from {(level, sample, allele_type, length_bin, region, stage): (tp_eval, fp, tp_truth, fn)}."""
    rows = [[*key, *values, ".", ".", "."] for key, values in counts.items()]
    pd.DataFrame(rows, columns=KEYS + ["tp_eval", "fp", "tp_truth", "fn", "precision", "recall", "f1"]).to_csv(
        path, sep="\t", index=False
    )


def strata(scale):
    """Counts for SNVs and 50-1999 bp deletions in two regions, at both stages, with every marginal."""
    counts = {}
    for stage, gained in (("exact", 0), ("truvari", 2)):
        for level, samples in (("site", ["all"]), ("call", ["S1", "all"])):
            for sample in samples:
                for allele_type, length_bin, region, base in (
                    ("snv", "0", "US/RM", (90, 10, 80, 20)),
                    ("del", "50-1999", "SD/SR", (40, 10, 30, 10)),
                ):
                    tp_eval, fp, tp_truth, fn = (scale * value for value in base)
                    values = (tp_eval + gained, fp - gained, tp_truth + gained, fn - gained)
                    for key in (
                        (allele_type, length_bin, region),
                        (allele_type, length_bin, "all"),
                        (allele_type, "all", region),
                        (allele_type, "all", "all"),
                    ):
                        counts[(level, sample, *key, stage)] = values
                totals = [sum(counts[(level, sample, t, "all", "all", stage)][i] for t in ("snv", "del")) for i in range(4)]
                counts[(level, sample, "all", "all", "all", stage)] = tuple(totals)
    return counts


def test_sums_contigs_recomputes_metrics_and_plots(tmp_path):
    first, second = tmp_path / "chr1.tsv", tmp_path / "chr2.tsv"
    accuracy_tsv(first, strata(1))
    accuracy_tsv(second, strata(2))
    combined, per_contig, plots = tmp_path / "combined.tsv", tmp_path / "per_contig.tsv", tmp_path / "plots"
    subprocess.run(
        [
            sys.executable,
            str(REPO_ROOT / SCRIPT),
            "--tsvs",
            str(first),
            str(second),
            "--contigs",
            "chr1",
            "chr2",
            "--combined_tsv",
            str(combined),
            "--per_contig_tsv",
            str(per_contig),
            "--output_dir",
            str(plots),
        ],
        check=True,
    )

    table = pd.read_csv(combined, sep="\t", dtype=str, keep_default_na=False).set_index(KEYS)
    row = table.loc[("site", "all", "snv", "0", "US/RM", "truvari")]
    # 3x the single-contig counts, with the 2 truvari gains of each contig
    assert [int(row[c]) for c in ("tp_eval", "fp", "tp_truth", "fn")] == [274, 26, 244, 56]
    assert row["precision"] == f"{274 / 300:.4f}" and row["recall"] == f"{244 / 300:.4f}"

    stacked = pd.read_csv(per_contig, sep="\t", dtype=str, keep_default_na=False)
    assert list(stacked.columns[:2]) == ["contig", "level"]
    assert set(stacked["contig"]) == {"chr1", "chr2"} and len(stacked) == 2 * len(table)

    for name in ("accuracy_sites.png", "accuracy_calls.png", "accuracy_by_sample.png", "accuracy_by_stage.png"):
        assert (plots / name).stat().st_size > 0
    summary = pd.read_csv(plots / "summary.tsv", sep="\t", dtype=str, keep_default_na=False)
    del_row = summary[(summary["allele_type"] == "del") & (summary["length_bin"] == "50-1999") & (summary["region"] == "all")]
    assert del_row[["site_tp_eval", "call_tp_eval"]].values.tolist() == [["124", "124"]]
