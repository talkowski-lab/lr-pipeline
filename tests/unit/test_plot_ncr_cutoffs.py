import itertools
import subprocess
import sys

import pandas as pd
from conftest import REPO_ROOT, load_script

SCRIPT = "scripts/qc/plot_ncr_cutoffs.py"
mod = load_script(SCRIPT)

CUTOFFS = [round(i * 0.1, 1) for i in range(11)]
COLUMNS = [
    "prec_numerator",
    "prec_denominator",
    "prec_skipped",
    "rec_numerator",
    "rec_denominator",
    "sites_passing",
    "sites_total",
]


def counts_tsv(path, calls):
    """An EvaluateNCRCutoffs TSV of PASS sites, with every `all` marginal, from per-stratum calls.

    `calls` maps (allele_type, size_bin, region) to (ncr, inherited, de novo, all-sample calls, sites) tuples.
    """
    by_stratum = {}
    for (allele_type, size_bin, region), values in calls.items():
        for stratum in itertools.product([allele_type, "all"], [size_bin, "all"], ["PASS", "all"], [region, "all"]):
            by_stratum.setdefault(stratum, []).extend(values)
    rows = []
    for stratum, values in by_stratum.items():
        total_calls = sum(value[3] for value in values)
        total_sites = sum(value[4] for value in values)
        for cutoff in CUTOFFS:
            kept = [value for value in values if value[0] <= cutoff + 1e-9]
            inherited = sum(value[1] for value in kept)
            rows.append(
                [
                    *stratum,
                    f"{cutoff:g}",
                    inherited,
                    inherited + sum(value[2] for value in kept),
                    0,
                    sum(value[3] for value in kept),
                    total_calls,
                    sum(value[4] for value in kept),
                    total_sites,
                ]
            )
    table = pd.DataFrame(rows, columns=["allele_type", "size_bin", "filter", "region", "cutoff", *COLUMNS])
    table.to_csv(path, sep="\t", index=False)


# De novo calls pile up at high NCR, so the knee sits below them
DENSE = [(0.0, 400, 2, 500, 400), (0.3, 50, 3, 60, 50), (0.8, 5, 20, 30, 25)]
SPARSE = [(0.0, 5, 0, 6, 5), (0.6, 1, 1, 2, 2)]


def run(tmp_path, calls, *extra):
    tsv = tmp_path / "counts.tsv"
    counts_tsv(tsv, calls)
    recommended = tmp_path / "recommended.tsv"
    subprocess.run(
        [
            sys.executable,
            str(REPO_ROOT / SCRIPT),
            "--tsv",
            str(tsv),
            "--output_dir",
            str(tmp_path / "plots"),
            "--recommended_tsv",
            str(recommended),
            *extra,
        ],
        check=True,
    )
    return pd.read_csv(recommended, sep="\t", dtype={"region": str})


def test_length_bounds():
    assert mod.length_bounds("0-0") == (0, 0)
    assert mod.length_bounds("50-499") == (50, 499)
    assert mod.length_bounds("500+") == (500, mod.MAX_LENGTH)


def test_recommendation_falls_back_from_region_to_all_regions_to_fallback(tmp_path):
    calls = {("del", "1-49", "US"): DENSE, ("del", "1-49", "RM"): SPARSE, ("del", "500+", "US"): SPARSE}
    recommended = run(tmp_path, calls, "--fallback_cutoff", "0.5")

    sources = {(row.allele_type, row.min_length, row.region): (row.source, row.cutoff) for row in recommended.itertuples()}
    assert sources[("del", 1, "US")][0] == "region"
    assert sources[("del", 1, "RM")][0] == "all_regions"
    assert sources[("del", 500, "US")] == ("fallback", 0.5)
    # Regions without sites in a stratum get no row
    assert ("del", 500, "RM") not in sources
    assert (tmp_path / "plots" / "DEL.combined.png").exists()


def test_recommendation_leaves_out_sparse_strata_without_fallback(tmp_path):
    calls = {("del", "1-49", "US"): DENSE, ("del", "500+", "US"): SPARSE}
    recommended = run(tmp_path, calls)

    assert list(recommended["min_length"]) == [1]


def test_plots_allele_types_too_sparse_for_any_histogram(tmp_path):
    recommended = run(tmp_path, {("del", "1-49", "US"): SPARSE}, "--fallback_cutoff", "0.5")

    assert list(recommended["source"]) == ["fallback"]
    assert (tmp_path / "plots" / "DEL.combined.png").exists()


def test_recommended_rows_never_overlap(tmp_path):
    calls = {
        (allele_type, size_bin, region): DENSE
        for allele_type in ("del", "ins")
        for size_bin in ("1-49", "50-499", "500+")
        for region in ("US", "RM", "SD", "SR")
    }
    recommended = run(tmp_path, calls, "--fallback_cutoff", "0.5")

    assert len(recommended) == 24
    for first, second in itertools.combinations(recommended.itertuples(), 2):
        same_stratum = first.allele_type == second.allele_type and first.region == second.region
        assert not (same_stratum and first.min_length <= second.max_length and second.min_length <= first.max_length)
