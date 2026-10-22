import argparse
from types import SimpleNamespace

import pytest
from conftest import load_script

mod = load_script("scripts/mei/PALMER_transfer_annotations.py")


def intersection_line(overrides=None):
    fields = {
        0: "chr1",
        3: "100,200",
        4: "A_x,C_y",
        5: "300,-310",
        7: "99",
        8: "100",
        9: "mei1",
        10: "-300",
        11: "S1,S2,",
        12: "A",
        13: "AACC,AATT",
    }
    fields.update(overrides or {})
    return "\t".join(fields.get(i, ".") for i in range(14)) + "\n"


def test_parse_intersection_line_fields():
    parsed = mod.parse_intersection_line(intersection_line())
    assert parsed["chrom"] == "chr1"
    assert parsed["sv_positions"] == ["100", "200"]
    assert parsed["sv_refs"] == ["A_x", "C_y"]
    assert parsed["sv_lengths"] == ["300", "-310"]
    assert parsed["mei_start_0based"] == 99
    assert parsed["mei_pos_1based"] == 100
    assert parsed["mei_id"] == "mei1"
    assert parsed["palmer_ref"] == "A"
    assert parsed["palmer_alts"] == ["AACC", "AATT"]


def test_parse_intersection_line_takes_absolute_length_and_strips_trailing_sample_comma():
    parsed = mod.parse_intersection_line(intersection_line())
    assert parsed["mei_length"] == 300
    assert parsed["mei_samples"] == ["S1", "S2"]


def test_parse_intersection_line_single_alt():
    assert mod.parse_intersection_line(intersection_line({13: "AACC"}))["palmer_alts"] == ["AACC"]


@pytest.mark.parametrize(
    "target,mei,expected",
    [(100, 100, 1.0), (50, 100, 0.5), (100, 50, 0.5), (0, 10, 0.0), (10, 0, 0.0), (-5, 10, 0.0)],
)
def test_compute_size_similarity(target, mei, expected):
    assert mod.compute_size_similarity(target, mei) == expected


@pytest.mark.parametrize(
    "args,expected",
    [
        ((100, 100, 100, 100), 1.0),
        ((100, 100, 150, 100), 0.5),
        ((100, 100, 200, 100), 0.0),
        ((100, 100, 300, 100), 0.0),
        ((100, 100, 110, 50), 0.5),
        ((100, 50, 100, 100), 0.5),
    ],
)
def test_compute_reciprocal_overlap(args, expected):
    assert mod.compute_reciprocal_overlap(*args) == pytest.approx(expected)


def test_compute_sequence_similarity_identical_and_case_insensitive():
    assert mod.compute_sequence_similarity("ACGT", "acgt") == 1.0


def test_compute_sequence_similarity_empty():
    assert mod.compute_sequence_similarity("", "") == 1.0
    assert mod.compute_sequence_similarity("ACGT", "") == 0.0


def test_compute_sequence_similarity_edit_distance():
    assert mod.compute_sequence_similarity("ACGT", "ACGA") == pytest.approx(0.75)
    assert mod.compute_sequence_similarity("AAAA", "AA") == pytest.approx(0.5)


@pytest.mark.parametrize(
    "gt,expected",
    [(None, False), ((0, 0), False), ((0, 1), True), ((1, 1), True), ((None, 0), True), ((0,), False)],
)
def test_is_non_ref_genotype(gt, expected):
    assert mod.is_non_ref_genotype(gt) is expected


class FakeRecord:
    def __init__(self, pos=100, ref="A", alts=("AACC",), allele_length=3, allele_type="ins", genotypes=None):
        self.pos = pos
        self.ref = ref
        self.alts = alts
        self.info = {"allele_length": allele_length, "allele_type": allele_type}
        self.samples = {s: {"GT": gt} for s, gt in (genotypes or {"S1": (0, 1), "S2": (0, 0)}).items()}


MEI = {"mei_length": 4, "mei_start_0based": 99, "mei_pos_1based": 100, "mei_samples": ["S1"]}
ARGS = argparse.Namespace(
    reciprocal_overlap=0.5, size_similarity=0.5, sequence_similarity=0.5, breakpoint_window=10, min_shared_samples=1
)


def passes(record, args=ARGS, mei=MEI, sv_length=4):
    return mod.record_passes_filters(record, "A", sv_length, mei, ["AACC"], ["S1", "S2"], args)


def test_record_passes_filters_accepts_matching_insertion():
    assert passes(FakeRecord())


def test_record_passes_filters_rejects_wrong_allele_length():
    assert not passes(FakeRecord(allele_length=10))


def test_record_passes_filters_rejects_non_insertion():
    assert not passes(FakeRecord(allele_type="del"))


def test_record_passes_filters_rejects_ref_mismatch():
    assert not passes(FakeRecord(ref="G"))


def test_record_passes_filters_rejects_distant_breakpoint():
    assert not passes(FakeRecord(pos=500))


def test_record_passes_filters_requires_shared_samples():
    assert not passes(FakeRecord(genotypes={"S1": (0, 0), "S2": (0, 1)}))


def test_record_passes_filters_rejects_dissimilar_sequence():
    assert not passes(FakeRecord(alts=("TTTT",)), args=SimpleNamespace(**{**vars(ARGS), "sequence_similarity": 0.9}))
