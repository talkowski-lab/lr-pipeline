import subprocess
import sys

import pytest
from conftest import REPO_ROOT

SCRIPT = REPO_ROOT / "scripts/qc/vcf_qc/vcf2bed.py"

# Covers every branch of svtk vcf2bed: symbolic and sequence alleles, INS, BND by bracket ALT and by CHR2/END2, CTX,
# CPX with and without an insertion class, a CNV called through CN, a multi-value FILTER and a missing FILTER
VCF = (
    "##fileformat=VCFv4.2\n"
    "##contig=<ID=chr1,length=248956422>\n"
    "##contig=<ID=chr2,length=242193529>\n"
    '##ALT=<ID=DEL,Description="Deletion">\n'
    '##ALT=<ID=DUP,Description="Duplication">\n'
    '##ALT=<ID=INV,Description="Inversion">\n'
    '##ALT=<ID=INS,Description="Insertion">\n'
    '##ALT=<ID=CPX,Description="Complex">\n'
    '##ALT=<ID=CTX,Description="Translocation">\n'
    '##ALT=<ID=CNV,Description="CNV">\n'
    '##FILTER=<ID=LowQual,Description="low">\n'
    '##FILTER=<ID=UNRESOLVED,Description="unres">\n'
    '##INFO=<ID=END,Number=1,Type=Integer,Description="end">\n'
    '##INFO=<ID=SVTYPE,Number=1,Type=String,Description="type">\n'
    '##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="len">\n'
    '##INFO=<ID=CHR2,Number=1,Type=String,Description="chr2">\n'
    '##INFO=<ID=END2,Number=1,Type=Integer,Description="end2">\n'
    '##INFO=<ID=CPX_INTERVALS,Number=.,Type=String,Description="cpx">\n'
    '##INFO=<ID=CPX_TYPE,Number=1,Type=String,Description="cpxtype">\n'
    '##INFO=<ID=AF,Number=A,Type=Float,Description="af">\n'
    '##INFO=<ID=FLAGX,Number=0,Type=Flag,Description="flag">\n'
    '##FORMAT=<ID=GT,Number=1,Type=String,Description="gt">\n'
    '##FORMAT=<ID=CN,Number=1,Type=Integer,Description="cn">\n'
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2\tS3\n"
    "chr1\t1000\tdel1\tN\t<DEL>\t.\tPASS\tEND=2000;SVTYPE=DEL;SVLEN=-1000;AF=0.5\tGT\t0/1\t0/0\t./.\n"
    "chr1\t1500\tdel2\tA\tACGT\t.\tLowQual\tEND=1500;SVTYPE=INS;SVLEN=3;FLAGX\tGT\t1/1\t0/1\t0/0\n"
    "chr1\t3000\tdup1\tN\t<DUP>\t.\tPASS;LowQual\tEND=5000;SVTYPE=DUP;SVLEN=2000\tGT\t0|1\t1|0\t0/0\n"
    "chr1\t6000\tinv1\tN\t<INV>\t.\t.\tEND=6000;SVTYPE=INV;SVLEN=0\tGT\t0/1\t./.\t0/0\n"
    "chr1\t7000\tins1\tN\t<INS>\t.\tPASS\tEND=7001;SVTYPE=INS;SVLEN=300\tGT\t0/1\t0/0\t0/0\n"
    "chr1\t8000\tbnd1\tN\tN]chr2:5000]\t.\tPASS\tSVTYPE=BND;CHR2=chr2;END2=5000\tGT\t0/1\t0/0\t0/0\n"
    "chr1\t8500\tbnd2\tN\t<BND>\t.\tPASS\tEND=8500;SVTYPE=BND;CHR2=chr2;END2=9000\tGT\t0/1\t0/0\t0/0\n"
    "chr1\t9000\tctx1\tN\t<CTX>\t.\tPASS\tEND=12000;SVTYPE=CTX;CHR2=chr2\tGT\t0/1\t0/0\t0/0\n"
    "chr1\t10000\tcpx1\tN\t<CPX>\t.\tPASS\tEND=11000;SVTYPE=CPX;CPX_TYPE=dDUP;CPX_INTERVALS=DUP_chr2:100-300,INV_chr2:300-400\tGT\t0/1\t0/0\t0/0\n"
    "chr1\t10500\tcpx2\tN\t<CPX>\t.\tPASS\tEND=11500;SVTYPE=CPX;CPX_TYPE=delINV;CPX_INTERVALS=DEL_chr1:10500-10700,INV_chr1:10700-11500\tGT\t0/1\t0/0\t0/0\n"
    "chr1\t13000\tcnv1\tN\t<CNV>\t.\tPASS\tEND=15000;SVTYPE=CNV;SVLEN=2000\tGT:CN\t./.:3\t./.:2\t./.:1\n"
    "chr1\t16000\tsnv1\tA\tG\t.\tPASS\tEND=16000\tGT\t0/1\t0/0\t0/0\n"
)

# Expected outputs were produced by `svtk vcf2bed` in the GATK-SV sv-pipeline image with the same arguments
EXPECTED = {
    ("-i", "SVTYPE", "-i", "SVLEN", "--split-cpx"): (
        "#chrom\tstart\tend\tname\tsvtype\tsamples\tSVTYPE\tSVLEN\n"
        "chr1\t999\t2000\tdel1\tDEL\tS1\tDEL\t-1000\n"
        "chr1\t1499\t1500\tdel2\tACGT\tS1,S2\tINS\t3\n"
        "chr1\t2999\t5000\tdup1\tDUP\tS1,S2\tDUP\t2000\n"
        "chr1\t5999\t6000\tinv1\tINV\tS1\tINV\t0\n"
        "chr1\t6999\t7001\tins1\tINS\tS1\tINS\t300\n"
        "chr1\t7999\t8000\tbnd1\tN]chr2:5000]\tS1\tBND\tNA\n"
        "chr1\t8499\t8500\tbnd2\tBND\tS1\tBND\tNA\n"
        "chr1\t8999\t9000\tctx1\tCTX\tS1\tCTX\tNA\n"
        "chr2\t99\t300\tcpx1\tDUP\tS1\tCPX\tNA\n"
        "chr2\t299\t400\tcpx1\tINV\tS1\tCPX\tNA\n"
        "chr1\t9999\t10000\tcpx1\tDEL\tS1\tCPX\tNA\n"
        "chr1\t10499\t10700\tcpx2\tDEL\tS1\tCPX\tNA\n"
        "chr1\t10699\t11500\tcpx2\tINV\tS1\tCPX\tNA\n"
        "chr1\t12999\t15000\tcnv1\tCNV\tS1,S3\tCNV\t2000\n"
        "chr1\t15999\t16000\tsnv1\tG\tS1\tNA\tNA\n"
    ),
    ("--info", "SVLEN", "--info", "END", "--no-samples"): (
        "#chrom\tstart\tend\tname\tsvtype\tSVLEN\tEND\n"
        "chr1\t999\t2000\tdel1\tDEL\t-1000\t2000\n"
        "chr1\t1499\t1500\tdel2\tACGT\t3\t1500\n"
        "chr1\t2999\t5000\tdup1\tDUP\t2000\t5000\n"
        "chr1\t5999\t6000\tinv1\tINV\t0\t6000\n"
        "chr1\t6999\t7001\tins1\tINS\t300\t7001\n"
        "chr1\t7999\t8000\tbnd1\tN]chr2:5000]\tNA\t8000\n"
        "chr1\t8499\t8500\tbnd2\tBND\tNA\t8500\n"
        "chr1\t8999\t9000\tctx1\tCTX\tNA\t12000\n"
        "chr1\t9999\t11000\tcpx1\tCPX\tNA\t11000\n"
        "chr1\t10499\t11500\tcpx2\tCPX\tNA\t11500\n"
        "chr1\t12999\t15000\tcnv1\tCNV\t2000\t15000\n"
        "chr1\t15999\t16000\tsnv1\tG\tNA\t16000\n"
    ),
    ("--info", "SVTYPE", "--include-filters", "--split-bnd", "--split-cpx"): (
        "#chrom\tstart\tend\tname\tsvtype\tsamples\tSVTYPE\tFILTER\n"
        "chr1\t999\t2000\tdel1\tDEL\tS1\tDEL\tPASS\n"
        "chr1\t1499\t1500\tdel2\tACGT\tS1,S2\tINS\tLowQual\n"
        "chr1\t2999\t5000\tdup1\tDUP\tS1,S2\tDUP\tPASS,LowQual\n"
        "chr1\t5999\t6000\tinv1\tINV\tS1\tINV\t\n"
        "chr1\t6999\t7001\tins1\tINS\tS1\tINS\tPASS\n"
        "chr1\t7999\t8000\tbnd1\tN]chr2:5000]\tS1\tBND\tPASS\n"
        "chr2\t4999\t5000\tbnd1\tN]chr2:5000]\tS1\tBND\tPASS\n"
        "chr1\t8499\t8500\tbnd2\tBND\tS1\tBND\tPASS\n"
        "chr2\t8999\t9000\tbnd2\tBND\tS1\tBND\tPASS\n"
        "chr1\t8999\t9000\tctx1\tCTX\tS1\tCTX\tPASS\n"
        "chr2\t11999\t12000\tctx1\tCTX\tS1\tCTX\tPASS\n"
        "chr2\t99\t300\tcpx1\tDUP\tS1\tCPX\tPASS\n"
        "chr2\t299\t400\tcpx1\tINV\tS1\tCPX\tPASS\n"
        "chr1\t9999\t10000\tcpx1\tDEL\tS1\tCPX\tPASS\n"
        "chr1\t10499\t10700\tcpx2\tDEL\tS1\tCPX\tPASS\n"
        "chr1\t10699\t11500\tcpx2\tINV\tS1\tCPX\tPASS\n"
        "chr1\t12999\t15000\tcnv1\tCNV\tS1,S3\tCNV\tPASS\n"
        "chr1\t15999\t16000\tsnv1\tG\tS1\tNA\tPASS\n"
    ),
}


@pytest.mark.parametrize("args", list(EXPECTED))
def test_vcf2bed_matches_svtk(tmp_path, args):
    vcf_path = tmp_path / "in.vcf"
    vcf_path.write_text(VCF)
    bed_path = tmp_path / "out.bed"
    subprocess.run([sys.executable, str(SCRIPT), *args, str(vcf_path), str(bed_path)], check=True)
    assert bed_path.read_text() == EXPECTED[args]


def test_vcf2bed_info_all_lists_header_fields(tmp_path):
    vcf_path = tmp_path / "in.vcf"
    vcf_path.write_text(VCF)
    bed_path = tmp_path / "out.bed"
    subprocess.run([sys.executable, str(SCRIPT), "--info", "ALL", str(vcf_path), str(bed_path)], check=True)
    header = bed_path.read_text().splitlines()[0].split("\t")
    assert header[6:] == ["END", "SVTYPE", "SVLEN", "CHR2", "END2", "CPX_INTERVALS", "CPX_TYPE", "AF", "FLAGX"]
