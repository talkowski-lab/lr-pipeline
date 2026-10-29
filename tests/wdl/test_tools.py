import pysam
import pytest
from conftest import CONTIG, FAI, SAMPLES, TOOLS, UTILS, WINDOW_FASTA, load_records, rewrite_vcf, site_key, write_vcf
from pysam import bcftools

pytestmark = pytest.mark.wdl


def subset_samples(src, dst, samples):
    """Write src restricted to samples as a bgzipped, indexed VCF and return its path."""
    bcftools.view("-s", ",".join(samples), "-Oz", "-o", str(dst), str(src), catch_stdout=False)
    pysam.tabix_index(str(dst), preset="vcf", force=True)
    return dst


def is_phased_het(record, sample):
    gt = record["gts"][sample]
    return record["phased"][sample] and None not in gt and len(gt) == 2 and gt[0] != gt[1]


def is_snv_allele(ref, alt):
    """True when ref and alt differ by one base once their shared prefix and suffix are trimmed."""
    while len(ref) > 1 and len(alt) > 1 and ref[-1] == alt[-1]:
        ref, alt = ref[:-1], alt[:-1]
    while len(ref) > 1 and len(alt) > 1 and ref[0] == alt[0]:
        ref, alt = ref[1:], alt[1:]
    return len(ref) == 1 and len(alt) == 1


def trim_shared_suffix(record):
    """REF and ALTs with the trailing bases they all share removed, leaving at least one base in each."""
    alleles = [record["ref"], *record["alts"]]
    while min(len(a) for a in alleles) > 1 and len({a[-1] for a in alleles}) == 1:
        alleles = [a[:-1] for a in alleles]
    return alleles[0], tuple(alleles[1:])


def test_backbone_phase_flips_phase_sets_back_to_the_base_and_reports_unmatched_samples(
    run_workflow, run_dir, subset_vcf, utils_docker, source
):
    sample = "HG00097"
    samples, records = source

    def swap(record):
        call = record.samples[sample]
        gt = call["GT"]
        if call.phased and None not in gt and len(gt) == 2 and gt[0] != gt[1]:
            call["GT"] = gt[::-1]
            call.phased = True

    target = rewrite_vcf(subset_vcf, run_dir / "target.vcf.gz", [], swap)
    base = subset_samples(subset_vcf, run_dir / "base.vcf.gz", SAMPLES[:3])
    outputs = run_workflow(
        f"{TOOLS}/BackbonePhase.wdl",
        {
            "vcf": str(target),
            "vcf_idx": f"{target}.tbi",
            "base_vcfs": [str(base)],
            "base_vcf_idxs": [f"{base}.tbi"],
            "contig": CONTIG,
            "prefix": "backbone",
            "utils_docker": utils_docker,
        },
    )

    # The base holds the original phasing, so every SNV in a swapped phase set disagrees with it and the set flips back
    evidenced = {
        r["format"][sample]["PS"]
        for r in records
        if is_phased_het(r, sample)
        and r["format"][sample].get("PS") is not None
        and any(a != 0 and is_snv_allele(r["ref"], r["alts"][a - 1]) for a in r["gts"][sample])
    }
    _, got = load_records(outputs["BackbonePhase.transferred_vcf"])
    assert [site_key(r) for r in got] == [site_key(r) for r in records]
    restored = left_swapped = 0
    for original, record in zip(records, got):
        for other in samples:
            if other != sample:
                assert record["gts"][other] == original["gts"][other]
        if not is_phased_het(original, sample):
            assert record["gts"][sample] == original["gts"][sample]
        elif original["format"][sample].get("PS") in evidenced:
            assert record["gts"][sample] == original["gts"][sample], original["id"]
            restored += 1
        else:
            assert record["gts"][sample] == original["gts"][sample][::-1], original["id"]
            left_swapped += 1
    assert restored and left_swapped
    with open(outputs["BackbonePhase.missing_samples"]) as handle:
        assert handle.read().split() == sorted(SAMPLES[3:])


def test_merge_hiphase_callsets_merges_disjoint_sample_sets_without_trgt(run_workflow, run_dir, subset_vcf, utils_docker, source):
    parts = [
        subset_samples(subset_vcf, run_dir / "first.vcf.gz", SAMPLES[:3]),
        subset_samples(subset_vcf, run_dir / "second.vcf.gz", SAMPLES[3:]),
    ]
    outputs = run_workflow(
        f"{UTILS}/MergeHiPhaseCallsets.wdl",
        {
            "phased_vcfs": [str(p) for p in parts],
            "phased_vcf_idxs": [f"{p}.tbi" for p in parts],
            "ref_fa": str(WINDOW_FASTA),
            "ref_fai": str(FAI),
            "contigs": [CONTIG],
            "prefix": "hiphase",
            "merge_trgt": False,
            "trgt_docker": utils_docker,
            "utils_docker": utils_docker,
        },
    )

    assert outputs.get("MergeHiPhaseCallsets.hiphase_merged_trgt_vcf") is None
    samples, records = source
    got_samples, got = load_records(outputs["MergeHiPhaseCallsets.hiphase_merged_integrated_vcf"])
    assert got_samples == samples
    # bcftools merge reorders records that share a position and trims the bases REF and every ALT end with, so a
    # trimmed record can share its site with another; records are matched by their unique ID instead
    merged = {r["id"]: r for r in got}
    assert sorted(merged) == sorted(r["id"] for r in records)
    for original in records:
        record = merged[original["id"]]
        assert (record["pos"], record["ref"], record["alts"]) == (original["pos"], *trim_shared_suffix(original))
        assert record["gts"] == original["gts"]
        assert record["phased"] == original["phased"]


def test_merge_palmer_callsets_merges_disjoint_samples_filling_absent_sites_with_ref(
    run_workflow, run_dir, subset_vcf, utils_docker, source
):
    me_types = ["ALU", "LINE1", "SVA"]
    insertions = [r for r in source[1] if r["info"]["allele_type"] == "ins" and len(r["alts"]) == 1]
    me_type = {site_key(r): me_types[i % len(me_types)] for i, r in enumerate(insertions)}
    first, second = insertions[: len(insertions) * 2 // 3], insertions[len(insertions) // 3 :]
    assert len({site_key(r) for r in first} & {site_key(r) for r in second}) > 10

    # PALMER callsets carry INFO/ME_TYPE, which the merge builds the record IDs from
    with pysam.VariantFile(str(subset_vcf)) as vcf:
        header = vcf.header.copy()
    header.info.add("ME_TYPE", 1, "String", "Mobile element type")
    tagged = run_dir / "tagged.vcf.gz"
    plain = tagged.with_suffix("")
    with pysam.VariantFile(str(subset_vcf)) as vcf_in, pysam.VariantFile(str(plain), "w", header=header) as vcf_out:
        for record in vcf_in:
            key = (record.contig, record.pos, record.ref, record.alts)
            if key in me_type:
                record.translate(header)
                record.info["ME_TYPE"] = me_type[key]
                vcf_out.write(record)
    pysam.tabix_index(str(plain), preset="vcf", force=True)

    callsets = {"A": (first, SAMPLES[:3]), "B": (second, SAMPLES[3:])}
    paths = []
    for name, (records, samples) in callsets.items():
        keep = {site_key(r) for r in records}
        full = write_vcf(tagged, run_dir / f"{name}.all.vcf.gz", lambda r, keep=keep: (r.contig, r.pos, r.ref, r.alts) in keep)
        paths.append(subset_samples(full, run_dir / f"{name}.vcf.gz", samples))
    outputs = run_workflow(
        f"{TOOLS}/MergePALMERCallsets.wdl",
        {
            "vcfs": [str(p) for p in paths],
            "vcf_idxs": [f"{p}.tbi" for p in paths],
            "contigs": [CONTIG],
            "prefix": "palmer",
            "utils_docker": utils_docker,
        },
    )

    got_samples, got = load_records(outputs["MergePALMERCallsets.palmer_merged_vcf"])
    assert got_samples == SAMPLES
    # bcftools merge may reorder records that share a position, so records are compared by site
    merged = {site_key(r): r for r in got}
    assert sorted(merged) == sorted({site_key(r) for r in first + second})
    for original in insertions:
        record = merged[site_key(original)]
        assert record["id"] == f"{me_type[site_key(original)]}_{CONTIG}_{original['pos']}_{original['info']['allele_length']}"
        for records, samples in callsets.values():
            present = site_key(original) in {site_key(r) for r in records}
            for sample in samples:
                expected = original["gts"][sample] if present else (0, 0)
                assert record["gts"][sample] == expected, (original["id"], sample)
