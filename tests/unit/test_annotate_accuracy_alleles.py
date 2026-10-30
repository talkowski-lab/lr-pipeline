import subprocess
import sys

import pysam
from conftest import REPO_ROOT

SCRIPT = "scripts/qc/annotate_accuracy_alleles.py"

HEADER = (
    "##fileformat=VCFv4.2\n"
    "##contig=<ID=chr22,length=50818468>\n"
    '##INFO=<ID=SVTYPE,Number=1,Type=String,Description="SV type">\n'
    '##INFO=<ID=SVLEN,Number=.,Type=Integer,Description="SV length">\n'
    '##INFO=<ID=END,Number=1,Type=Integer,Description="End position">\n'
    '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\n"
)


def annotate(tmp_path, records, *extra):
    vcf_path = tmp_path / "in.vcf"
    vcf_path.write_text(
        HEADER + "".join(f"chr22\t{pos}\t{vid}\t{ref}\t{alt}\t.\t.\t{info}\tGT\t0/1\n" for pos, vid, ref, alt, info in records)
    )
    out_path = tmp_path / "out.vcf.gz"
    subprocess.run(
        [
            sys.executable,
            str(REPO_ROOT / SCRIPT),
            "--vcf",
            str(vcf_path),
            "--output",
            str(out_path),
            "--id_prefix",
            "eval",
            *extra,
        ],
        check=True,
    )
    assert (tmp_path / "out.vcf.gz.tbi").exists()
    with pysam.VariantFile(str(out_path)) as vcf:
        return [(r.id, r.info["allele_type"], r.info["allele_length"], r.info.get("SVTYPE")) for r in vcf]


def test_types_and_sizes_from_ref_and_alt_or_sv_fields(tmp_path):
    records = [
        (100, "snv", "A", "G", "."),
        (200, "mnv", "AC", "GT", "."),
        (300, "ins", "A", "ACGT", "."),
        (400, "del", "ACGT", "A", "."),
        # A duplication written as an insertion counts as one; a symbolic one stays a duplication
        (500, "seqdup", "A", "A" + "C" * 60, "SVTYPE=DUP;SVLEN=60"),
        (600, "symdup", "A", "<DUP>", "SVTYPE=DUP;SVLEN=1000;END=1600"),
        (700, "symdel", "A", "<DEL>", "SVTYPE=DEL;SVLEN=-300;END=1000"),
    ]
    assert annotate(tmp_path, records) == [
        ("snv", "snv", 0, None),
        ("mnv", "mnv", 0, None),
        ("ins", "ins", 3, None),
        ("del", "del", -3, None),
        ("seqdup", "ins", 60, "INS"),
        ("symdup", "dup", 1000, "DUP"),
        ("symdel", "del", -300, "DEL"),
    ]


def test_keeps_existing_type_and_length(tmp_path):
    vcf_path_records = [(100, "trv", "AT", "ATATAT", "allele_type=trv;allele_length=4")]
    header = HEADER.replace(
        "##FORMAT",
        '##INFO=<ID=allele_type,Number=1,Type=String,Description="Allele type">\n'
        '##INFO=<ID=allele_length,Number=1,Type=Integer,Description="Allele length">\n##FORMAT',
        1,
    )
    vcf_path = tmp_path / "in.vcf"
    vcf_path.write_text(
        header + "".join(f"chr22\t{p}\t{i}\t{r}\t{a}\t.\t.\t{info}\tGT\t0/1\n" for p, i, r, a, info in vcf_path_records)
    )
    out_path = tmp_path / "out.vcf.gz"
    subprocess.run(
        [sys.executable, str(REPO_ROOT / SCRIPT), "--vcf", str(vcf_path), "--output", str(out_path), "--id_prefix", "eval"],
        check=True,
    )
    record = next(pysam.VariantFile(str(out_path)))
    assert (record.info["allele_type"], record.info["allele_length"]) == ("trv", 4)


def test_renames_only_missing_and_repeated_ids_unless_told_to_replace(tmp_path):
    records = [(100, ".", "A", "G", "."), (200, "v1", "A", "G", "."), (300, "v1", "A", "G", "."), (400, "v2", "A", "G", ".")]
    assert [r[0] for r in annotate(tmp_path, records)] == ["eval_1", "v1", "eval_3", "v2"]
    assert [r[0] for r in annotate(tmp_path, records, "--replace_ids")] == ["eval_1", "eval_2", "eval_3", "eval_4"]
