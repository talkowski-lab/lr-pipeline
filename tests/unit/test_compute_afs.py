import subprocess
import sys

import pysam
import pytest
from conftest import REPO_ROOT, load_script

SCRIPT = "scripts/annotation/compute_AFs.py"
mod = load_script(SCRIPT)


def test_resolve_lps_alleles_ignores_haplotype_order():
    # Allele 1 is fixed by a homozygous carrier, then propagates through hets listed in either order
    carriers = [((1, 1), [5, 5]), ((1, 2), [7, 5]), ((3, 2), [7, 9]), ((0, 0), [3, 3])]
    assert mod.resolve_lps_alleles(carriers, 5) == [3, 5, 7, 9, None]


def test_resolve_lps_alleles_leaves_unpinned_hets_missing():
    assert mod.resolve_lps_alleles([((1, 2), [4, 6])], 3) == [None, None, None]


@pytest.mark.parametrize("header_last", [False, True])
def test_lps_allele_written_from_table(tmp_path, header_last):
    vcf_path = tmp_path / "in.vcf"
    vcf_path.write_text(
        "##fileformat=VCFv4.2\n"
        "##contig=<ID=chr22,length=50818468>\n"
        '##INFO=<ID=TRID,Number=1,Type=String,Description="TR ID">\n'
        '##INFO=<ID=MOTIFS,Number=.,Type=String,Description="Motifs">\n'
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2\tS3\tS4\n"
        "chr22\t100\tchr22-100-TRV-8\tCATATAT\tCATAT,CATATATAT\t.\t.\tTRID=22-100-108-AT;MOTIFS=AT\tGT\t0/0\t1|2\t2/2\t1/1\n"
    )
    lps_path = tmp_path / "lps.tsv"
    header = "trid\tmotif\tS1\tS2\tS3\tS4\n"
    rows = "22-100-108-AT\tAT\t3,3\t4,2\t4,4\t.,.\n22-100-108-AT\tGC\t1,1\t1,1\t1,1\t1,.\n"
    lps_path.write_text(rows + header if header_last else header + rows)
    out_path = tmp_path / "out.vcf"
    subprocess.run([sys.executable, str(REPO_ROOT / SCRIPT), str(vcf_path), str(out_path), "-l", str(lps_path)], check=True)

    record = next(pysam.VariantFile(str(out_path)))
    assert record.info["LPS_allele"] == (3, 2, 4)
