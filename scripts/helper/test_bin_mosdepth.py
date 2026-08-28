import gzip
import tempfile
import unittest
from pathlib import Path

from bin_mosdepth import bin_mosdepth


class BinMosdepthTest(unittest.TestCase):
    def write_input(self, directory, contents):
        path = Path(directory) / "coverage.bed.gz"
        with gzip.open(path, "wt") as output:
            output.write(contents)
        return path

    def test_full_bins_use_base_weighted_median(self):
        with tempfile.TemporaryDirectory() as directory:
            source = self.write_input(
                directory,
                "chr1\t0\t2\t10\nchr1\t2\t4\t20\n"
                "chr1\t4\t6\t30\nchr1\t6\t8\t40\n",
            )
            output = Path(directory) / "binned.bed"
            histogram = bin_mosdepth(
                source, output, 4, truncate_depth=True
            )
            self.assertEqual(output.read_text(), "chr1\t0\t4\t15\nchr1\t4\t8\t35\n")
            self.assertEqual(histogram, {15: 1, 35: 1})

    def test_one_based_contig_output_drops_partial_bin(self):
        with tempfile.TemporaryDirectory() as directory:
            source = self.write_input(
                directory,
                "chr1\t0\t3\t1.9\nchr1\t3\t5\t2.1\n"
                "chr2\t0\t4\t8\n",
            )
            output = Path(directory) / "counts.tsv"
            bin_mosdepth(
                source,
                output,
                4,
                coordinate_system="one-based",
                contig="chr1",
                output_contig="chr1",
                truncate_depth=True,
            )
            self.assertEqual(output.read_text(), "chr1\t1\t4\t1\n")

    def test_contiguous_mode_keeps_terminal_partial_and_excludes_contig(self):
        with tempfile.TemporaryDirectory() as directory:
            source = self.write_input(
                directory,
                "# comment\nchr1\t0\t2\t10\nchr1\t2\t5\t20\n"
                "chrY\t0\t4\t99\n",
            )
            output = Path(directory) / "binned.tsv.gz"
            histogram = bin_mosdepth(
                source,
                output,
                4,
                require_contiguous=True,
                allow_partial=True,
                skip_comments=True,
                exclude_contigs={"chrY"},
            )
            with gzip.open(output, "rt") as source:
                self.assertEqual(
                    source.read(), "chr1\t0\t4\t15\nchr1\t4\t5\t20\n"
                )
            self.assertEqual(histogram, {15.0: 1, 20.0: 1})

    def test_contiguous_mode_rejects_gaps(self):
        with tempfile.TemporaryDirectory() as directory:
            source = self.write_input(
                directory, "chr1\t0\t2\t10\nchr1\t3\t4\t20\n"
            )
            with self.assertRaisesRegex(ValueError, "not contiguous"):
                bin_mosdepth(
                    source,
                    Path(directory) / "binned.tsv",
                    4,
                    require_contiguous=True,
                )


if __name__ == "__main__":
    unittest.main()
