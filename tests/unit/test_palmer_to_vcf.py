import pysam
import pytest
from conftest import load_script

mod = load_script("scripts/mei/PALMER_to_vcf.py")

TSD_HEADER = "cluster_id\tc1\tc2\tc3\tc4\tc5\tseq\n"


def call_line(cid, pos=(10, 10, 12, 12), tsd_5=(5, 6), tsd_3=(30, 31), conf=3, ori="+", polya=20, inv=0, extra=("0", "0")):
    fields = [
        cid,
        "chr1",
        *map(str, pos),
        *map(str, tsd_5),
        *map(str, tsd_3),
        str(conf),
        "x",
        "x",
        ori,
        str(polya),
        "7",
        "8",
        "9",
        str(inv),
    ]
    fields += list(extra)
    return "\t".join(fields) + "\n"


@pytest.fixture
def ref_fasta(tmp_path):
    fa = tmp_path / "ref.fa"
    fa.write_text(">chr1\n" + "ACGTACGTAC" * 3 + "\n")
    pysam.faidx(str(fa))
    handle = pysam.FastaFile(str(fa))
    yield handle
    handle.close()


def test_parse_tsd_reads_keeps_first_sequence_per_cluster_and_skips_header(tmp_path):
    path = tmp_path / "tsd.txt"
    path.write_text(TSD_HEADER + "c1\ta\tb\tc\td\te\tAAAA\nc1\ta\tb\tc\td\te\tCCCC\nc2\ta\tb\tc\td\te\tGGGG\n")
    assert mod.parse_tsd_reads(str(path)) == {"c1": "AAAA", "c2": "GGGG"}


def test_parse_palmer_calls_builds_ref_alt_and_length(tmp_path, ref_fasta):
    path = tmp_path / "calls.txt"
    path.write_text("cluster_id\tchrom\n" + call_line("c1", pos=(10, 10, 12, 12), tsd_5=(5, 7), tsd_3=(30, 36)))
    calls = mod.parse_palmer_calls(str(path), {"c1": "tttt"}, ref_fasta)
    call = calls["c1"]
    assert call["pos"] == 11
    assert call["ref"] == "A"
    assert call["alt"] == "ATTTT"
    assert call["length"] == 33 - 6 + 9
    assert call["5inv"] == 0


def test_parse_palmer_calls_adds_inversion_span_to_length(tmp_path, ref_fasta):
    path = tmp_path / "calls.txt"
    path.write_text(call_line("c1", tsd_5=(5, 5), tsd_3=(25, 25), inv=1, extra=("100", "80")))
    assert mod.parse_palmer_calls(str(path), {}, ref_fasta)["c1"]["length"] == (25 - 5) + 9 + (100 - 80)


def test_parse_palmer_calls_filters_low_confidence_and_duplicates(tmp_path, ref_fasta):
    path = tmp_path / "calls.txt"
    path.write_text(call_line("low", conf=0) + call_line("c1") + call_line("dup"))
    calls = mod.parse_palmer_calls(str(path), {}, ref_fasta)
    assert list(calls) == ["c1"]


def test_parse_palmer_calls_missing_tsd_sequence_gives_ref_only_alt(tmp_path, ref_fasta):
    path = tmp_path / "calls.txt"
    path.write_text(call_line("c1"))
    call = mod.parse_palmer_calls(str(path), {}, ref_fasta)["c1"]
    assert call["alt"] == call["ref"]


def test_write_vcf_header_lists_contigs_and_sample(tmp_path, capsys):
    fai = tmp_path / "ref.fa.fai"
    fai.write_text("chr1\t100\t0\t60\t61\nchr2\t50\t0\t60\t61\n")
    mod.write_vcf_header(str(fai), "SAMPLE1")
    lines = capsys.readouterr().out.splitlines()
    assert lines[0] == "##fileformat=VCFv4.2"
    assert "##contig=<ID=chr1,length=100>" in lines
    assert "##contig=<ID=chr2,length=50>" in lines
    assert lines[-1].endswith("FORMAT\tSAMPLE1")
    assert any(line.startswith("##INFO=<ID=ME_TYPE") for line in lines)


def test_write_vcf_records_round_trips_through_pysam(tmp_path, capsys, ref_fasta):
    calls_path = tmp_path / "calls.txt"
    calls_path.write_text(call_line("c1", inv=1, extra=("100", "80")))
    calls = mod.parse_palmer_calls(str(calls_path), {"c1": "ACGT"}, ref_fasta)
    fai = tmp_path / "ref.fa.fai"
    fai.write_text("chr1\t30\t0\t60\t61\n")
    mod.write_vcf_header(str(fai), "S1")
    mod.write_vcf_records(calls, "ALU", "1|0")
    vcf = tmp_path / "out.vcf"
    vcf.write_text(capsys.readouterr().out)
    record = next(iter(pysam.VariantFile(str(vcf))))
    assert record.id == "c1"
    assert record.info["ME_TYPE"] == ("ALU",)
    assert record.info["allele_type"] == ("ins",)
    assert record.info["INVERSION_5PRIME"] is True
    assert record.samples["S1"]["GT"] == (1, 0)
    assert record.samples["S1"].phased
