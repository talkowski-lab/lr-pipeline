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


def test_males_count_one_allele_outside_par_in_either_genotype_form(tmp_path):
    # Males M1-M4 carry hemizygous calls written after ploidy normalization (./0, ./1) or before it (0/0, 1/1)
    vcf_path = tmp_path / "in.vcf"
    vcf_path.write_text(
        "##fileformat=VCFv4.2\n"
        "##contig=<ID=chrX,length=156040895>\n"
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tM1\tM2\tM3\tM4\tF1\n"
        "chrX\t5000000\tv1\tA\tG\t.\t.\t.\tGT\t./1\t./0\t1/1\t0/0\t0/1\n"
        "chrX\t5000100\tv2\tA\tG\t.\t.\t.\tGT\t./1\t./1\t./.\t./0\t1/1\n"
    )
    ped_path = tmp_path / "samples.ped"
    ped_path.write_text(
        "".join(f"fam\t{s}\t0\t0\t{sex}\t0\n" for s, sex in [("M1", 1), ("M2", 1), ("M3", 1), ("M4", 1), ("F1", 2)])
    )
    par_path = tmp_path / "par.bed"
    par_path.write_text("chrX\t10000\t2781479\n")
    out_path = tmp_path / "out.vcf"
    subprocess.run(
        [sys.executable, str(REPO_ROOT / SCRIPT), str(vcf_path), str(out_path), "-f", str(ped_path), "--par", str(par_path)],
        check=True,
    )

    first, second = pysam.VariantFile(str(out_path))
    assert (first.info["AN_XY"], first.info["AC_XY"], first.info["nhemiref_XY"]) == (4, (2,), 2)
    assert (first.info["AN"], first.info["AC"]) == (6, (3,))
    assert (second.info["AN_XY"], second.info["AC_XY"], second.info["nhemiref_XY"]) == (3, (2,), 1)
    assert (second.info["AN"], second.info["AC"]) == (5, (4,))
    assert second.info["NCR_XY"] == pytest.approx(0.25)


def test_no_call_rate_expects_one_allele_per_male_outside_par(tmp_path):
    # Fully called males and a cleared female carry no missing alleles on chrY; one missing male allele is one of two
    vcf_path = tmp_path / "in.vcf"
    vcf_path.write_text(
        "##fileformat=VCFv4.2\n"
        "##contig=<ID=chrX,length=156040895>\n"
        "##contig=<ID=chrY,length=57227415>\n"
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tM1\tM2\tF1\n"
        "chrY\t5000000\tv1\tA\tG\t.\t.\t.\tGT\t./1\t1/1\t./.\n"
        "chrY\t5000100\tv2\tA\tG\t.\t.\t.\tGT\t./1\t./.\t./.\n"
        "chrX\t5000000\tv3\tA\tG\t.\t.\t.\tGT\t./0\t./1\t0/.\n"
        "chrX\t1000000\tv4\tA\tG\t.\t.\t.\tGT\t0/0\t0/1\t./.\n"
    )
    ped_path = tmp_path / "samples.ped"
    ped_path.write_text("".join(f"fam\t{s}\t0\t0\t{sex}\t0\n" for s, sex in [("M1", 1), ("M2", 1), ("F1", 2)]))
    par_path = tmp_path / "par.bed"
    par_path.write_text("chrX\t10000\t2781479\n")
    out_path = tmp_path / "out.vcf"
    subprocess.run(
        [sys.executable, str(REPO_ROOT / SCRIPT), str(vcf_path), str(out_path), "-f", str(ped_path), "--par", str(par_path)],
        check=True,
    )

    y_called, y_missing, x_nonpar, x_par = pysam.VariantFile(str(out_path))
    assert (y_called.info["NCR"], y_called.info["NCR_XY"]) == (0, 0)
    assert y_missing.info["NCR"] == pytest.approx(0.5)
    assert x_nonpar.info["NCR"] == pytest.approx(0.25)
    # Inside the PAR every sample is diploid
    assert x_par.info["NCR"] == pytest.approx(2 / 6, abs=1e-6)
