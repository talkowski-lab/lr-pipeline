import os
import subprocess
import sys

import numpy as np
import pandas as pd
import pytest
from conftest import REPO_ROOT, load_script

mod = load_script("scripts/helper/estimated_CN_denoising.py")
SCRIPT = REPO_ROOT / "scripts" / "helper" / "estimated_CN_denoising.py"


def test_calculate_median_mad_known_values():
    median, mad = mod.calculate_median_mad(np.array([1, 2, 3, 4, 100]))
    assert median == 3
    assert mad == 1


def test_calculate_median_mad_constant_has_zero_mad():
    assert mod.calculate_median_mad(np.array([2.0, 2.0, 2.0])) == (2.0, 0.0)


def test_main_writes_per_sample_and_contig_stats(tmp_path):
    binwise = tmp_path / "binwise.tsv"
    pd.DataFrame(
        {
            "#Chr": ["chr1", "chr1", "chr1", "chr2", "chr2", "chr2"],
            "Start": [0, 10, 20, 0, 10, 20],
            "End": [10, 20, 30, 10, 20, 30],
            "S1": [2, 2, 2, 1, 1, 4],
            "S2": [3, 2, 4, 2, 2, 2],
        }
    ).to_csv(binwise, sep="\t", index=False)
    estimated = tmp_path / "estimated.tsv"
    pd.DataFrame(
        {
            "sample_id": ["S1", "S2"],
            "chr1_CopyNumber": [2.0, 3.0],
            "chr2_CopyNumber": [1.5, 2.0],
        }
    ).to_csv(estimated, sep="\t", index=False)
    stats, pdf = tmp_path / "stats.tsv", tmp_path / "plots.pdf"
    env = dict(os.environ, MPLBACKEND="Agg")
    subprocess.run(
        [
            sys.executable,
            str(SCRIPT),
            "--binwise-copy-number",
            str(binwise),
            "--estimated-copy-number",
            str(estimated),
            "--output-stats",
            str(stats),
            "--output-pdf",
            str(pdf),
        ],
        check=True,
        env=env,
        capture_output=True,
    )
    result = pd.read_csv(stats, sep="\t").set_index(["Chr", "SampleID"])
    assert len(result) == 4
    assert result.loc[("chr1", "S2"), "BinwiseCopyNumber_median"] == 3
    assert result.loc[("chr1", "S2"), "BinwiseCopyNumber_mad"] == 1
    assert result.loc[("chr2", "S1"), "BinwiseCopyNumber_median"] == 1
    assert result.loc[("chr2", "S1"), "BinwiseCopyNumber_mad"] == 0
    assert result.loc[("chr2", "S1"), "CopyNumber"] == pytest.approx(1.5)
    assert pdf.read_bytes().startswith(b"%PDF")
