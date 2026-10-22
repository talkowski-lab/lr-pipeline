import collections
import csv
import gzip

import pysam
import pytest
from conftest import REPO_ROOT, WorkflowError

pytestmark = pytest.mark.wdl

CONTIG = "chr22"
CONTIG_LENGTH = 50818468
SAMPLES = ["HG00096", "HG00097", "HG00099", "HG00126", "HG00128"]
FAI = REPO_ROOT / "tests" / "data" / "GRCh38.chr22.fa.fai"
UTILS = "wdl/annotation_utils"


def load_records(path):
    """Read a VCF into plain dicts so assertions do not depend on pysam record lifetimes."""
    with pysam.VariantFile(str(path)) as vcf:
        samples = list(vcf.header.samples)
        records = [
            {
                "contig": r.contig,
                "pos": r.pos,
                "id": r.id,
                "ref": r.ref,
                "alts": r.alts,
                "filter": tuple(r.filter.keys()),
                "info": dict(r.info.items()),
                "gts": {s: r.samples[s]["GT"] for s in samples},
            }
            for r in vcf
        ]
    return samples, records


def raw_info_end(path):
    """Read INFO/END straight from the VCF text, since pysam derives END from the record when the tag is absent."""
    ends = []
    with gzip.open(path, "rt") as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            fields = line.split("\t")
            tags = dict(tag.partition("=")[::2] for tag in fields[7].split(";"))
            ends.append((int(fields[1]), len(fields[3]), tags.get("END")))
    return ends


def site_key(record):
    return (record["contig"], record["pos"], record["ref"], record["alts"])


def is_carrier(record, sample):
    return any(allele not in (0, None) for allele in record["gts"][sample])


def write_vcf(src, dst, keep):
    """Write the records of a VCF accepted by keep() to a bgzipped, indexed VCF and return its path."""
    plain = dst.with_suffix("")
    with pysam.VariantFile(str(src)) as vcf_in, pysam.VariantFile(str(plain), "w", header=vcf_in.header) as vcf_out:
        for record in vcf_in:
            if keep(record):
                vcf_out.write(record)
    pysam.tabix_index(str(plain), preset="vcf", force=True)
    return dst


@pytest.fixture(scope="module")
def source():
    return load_records(REPO_ROOT / "tests" / "data" / "hprc_hgsvc_v10.chr22.subset.vcf.gz")


@pytest.fixture
def vcf_inputs(subset_vcf, utils_docker):
    return {"vcf": str(subset_vcf), "vcf_idx": f"{subset_vcf}.tbi", "utils_docker": utils_docker}


def test_source_fixture_matches_documented_composition(source):
    samples, records = source
    assert samples == SAMPLES
    assert len(records) == 385
    assert all(any(r["info"]["AC"]) for r in records)


@pytest.mark.parametrize("records_per_shard", [None, 100])
def test_strip_genotypes_drops_samples_and_keeps_every_record(run_workflow, vcf_inputs, source, records_per_shard):
    inputs = {**vcf_inputs, "prefix": "stripped"}
    if records_per_shard:
        inputs["records_per_shard"] = records_per_shard
    outputs = run_workflow(f"{UTILS}/StripGenotypes.wdl", inputs)

    samples, records = load_records(outputs["StripGenotypes.dropped_vcf"])
    assert samples == []
    assert [site_key(r) for r in records] == [site_key(r) for r in source[1]]
    assert [r["info"]["AC"] for r in records] == [r["info"]["AC"] for r in source[1]]


def test_split_vcf_per_contig_splits_adds_no_geno_copy_and_rewrites_snv_ids(run_workflow, vcf_inputs, source):
    outputs = run_workflow(
        f"{UTILS}/SplitVcfPerContig.wdl",
        {**vcf_inputs, "prefix": "split", "contigs": [CONTIG], "create_no_geno": True, "modify_snv_ids": True},
    )

    samples, records = load_records(outputs["SplitVcfPerContig.contig_vcfs"][0])
    assert samples == SAMPLES
    assert [site_key(r) for r in records] == [site_key(r) for r in source[1]]
    for got, original in zip(records, source[1]):
        if original["info"]["allele_type"] == "snv":
            assert got["id"] == f"{original['contig']}-{original['pos']}-{original['ref']}-{original['alts'][0]}"
        else:
            assert got["id"] == original["id"]

    no_geno_samples, no_geno_records = load_records(outputs["SplitVcfPerContig.contig_no_geno_vcfs"][0])
    assert no_geno_samples == []
    assert len(no_geno_records) == len(source[1])


def test_split_vcf_per_contig_adds_missing_info_header_lines(run_workflow, vcf_inputs):
    outputs = run_workflow(
        f"{UTILS}/SplitVcfPerContig.wdl",
        {**vcf_inputs, "prefix": "hdr", "contigs": [CONTIG], "missing_info_header_fields": ["NEW_TAG_A", "NEW_TAG_B"]},
    )
    with pysam.VariantFile(outputs["SplitVcfPerContig.contig_vcfs"][0]) as vcf:
        assert {"NEW_TAG_A", "NEW_TAG_B"} <= set(vcf.header.info.keys())
    assert outputs["SplitVcfPerContig.contig_no_geno_vcfs"] == []


def test_concatenate_vcfs_across_contigs_joins_parts_in_order(run_workflow, run_dir, subset_vcf, utils_docker, source):
    split_at = source[1][len(source[1]) // 2]["pos"]
    low = write_vcf(subset_vcf, run_dir / "low.vcf.gz", lambda r: r.pos < split_at)
    high = write_vcf(subset_vcf, run_dir / "high.vcf.gz", lambda r: r.pos >= split_at)

    outputs = run_workflow(
        f"{UTILS}/ConcatenateVcfsAcrossContigs.wdl",
        {
            "vcfs": [str(low), str(high)],
            "vcf_idxs": [f"{low}.tbi", f"{high}.tbi"],
            "contigs": ["part_low", "part_high"],
            "prefix": "concat",
            "drop_genotypes": True,
            "utils_docker": utils_docker,
        },
    )

    samples, records = load_records(outputs["ConcatenateVcfsAcrossContigs.concat_vcf"])
    assert samples == []
    assert [site_key(r) for r in records] == [site_key(r) for r in source[1]]


def test_concatenate_vcfs_across_contigs_rejects_duplicate_contigs(run_workflow, vcf_inputs):
    with pytest.raises(WorkflowError):
        run_workflow(
            f"{UTILS}/ConcatenateVcfsAcrossContigs.wdl",
            {
                "vcfs": [vcf_inputs["vcf"]] * 2,
                "vcf_idxs": [vcf_inputs["vcf_idx"]] * 2,
                "contigs": [CONTIG, CONTIG],
                "prefix": "dup",
                "utils_docker": vcf_inputs["utils_docker"],
            },
        )


def test_extract_sample_vcfs_routes_carried_variants_by_length_and_source(run_workflow, vcf_inputs, source):
    samples_to_extract = ["HG00096", "HG00126"]
    min_sv_length = 50
    outputs = run_workflow(
        f"{UTILS}/ExtractSampleVcfs.wdl",
        {
            "cohort_vcf": vcf_inputs["vcf"],
            "cohort_vcf_idx": vcf_inputs["vcf_idx"],
            "contigs": [CONTIG],
            "prefix": "per_sample",
            "sample_ids": samples_to_extract,
            "min_sv_length": min_sv_length,
            "utils_docker": vcf_inputs["utils_docker"],
        },
    )

    for i, sample in enumerate(samples_to_extract):
        carried = [r for r in source[1] if is_carrier(r, sample)]
        # Records without allele_length, such as tandem repeats, satisfy neither bcftools expression
        sized = [r for r in carried if "allele_length" in r["info"]]
        expected_snv = [
            site_key(r) for r in sized if abs(r["info"]["allele_length"]) < min_sv_length and r["info"]["SOURCE"] == "DeepVariant"
        ]
        expected_sv = [
            site_key(r)
            for r in sized
            if abs(r["info"]["allele_length"]) >= min_sv_length and r["info"]["SOURCE"] == "HPRC_SV_Integration"
        ]

        snv_samples, snv_records = load_records(outputs["ExtractSampleVcfs.snv_indel_vcfs"][i])
        sv_samples, sv_records = load_records(outputs["ExtractSampleVcfs.sv_vcfs"][i])
        assert snv_samples == sv_samples == [sample]
        assert [site_key(r) for r in snv_records] == expected_snv
        assert [site_key(r) for r in sv_records] == expected_sv
        assert expected_snv and expected_sv
        assert all(is_carrier(r, sample) for r in snv_records + sv_records)


def test_annotate_tr_end_tags_sets_end_from_ref_length(run_workflow, vcf_inputs, source):
    outputs = run_workflow(f"{UTILS}/AnnotateTREndTags.wdl", {**vcf_inputs, "prefix": "end", "contigs": [CONTIG]})

    samples, records = load_records(outputs["AnnotateTREndTags.vcf_with_end"])
    assert samples == SAMPLES
    assert [site_key(r) for r in records] == [site_key(r) for r in source[1]]
    ends = raw_info_end(outputs["AnnotateTREndTags.vcf_with_end"])
    assert len(ends) == len(records)
    assert all(end == str(pos + ref_length - 1) for pos, ref_length, end in ends)
    assert [r["gts"] for r in records] == [r["gts"] for r in source[1]]


def test_extract_trv_length_deltas_matches_independent_tabulation(run_workflow, vcf_inputs, source):
    outputs = run_workflow(f"{UTILS}/ExtractTRVLengthDeltas.wdl", {**vcf_inputs, "prefix": "trv"})

    expected_rows = []
    for r in source[1]:
        if r["info"]["allele_type"] != "trv" or min(len(m) for m in r["info"]["MOTIFS"]) != 3:
            continue
        for alt, ac in zip(r["alts"], r["info"]["AC"]):
            expected_rows.append(
                [CONTIG, str(r["pos"]), r["info"]["TRID"], ",".join(r["info"]["MOTIFS"])]
                + [str(len(r["ref"])), str(len(alt)), str(len(alt) - len(r["ref"])), str(ac)]
            )
    assert expected_rows

    with gzip.open(outputs["ExtractTRVLengthDeltas.allele_deltas_tsv"], "rt") as handle:
        rows = list(csv.reader(handle, delimiter="\t"))
    assert rows[0] == ["contig", "pos", "trid", "motifs", "ref_length", "alt_length", "length_delta", "allele_count"]
    assert rows[1:] == expected_rows

    counts = collections.defaultdict(lambda: [0, 0])
    for row in expected_rows:
        counts[int(row[6])][0] += 1
        counts[int(row[6])][1] += int(row[7])
    with gzip.open(outputs["ExtractTRVLengthDeltas.length_deltas_tsv"], "rt") as handle:
        count_rows = list(csv.reader(handle, delimiter="\t"))
    assert count_rows[0] == ["contig", "length_delta", "n_alleles", "allele_count"]
    assert count_rows[1:] == [[CONTIG, str(d), str(n), str(a)] for d, (n, a) in sorted(counts.items())]


def read_base_table(path):
    with open(path) as handle:
        return {row["category"]: row for row in csv.DictReader(handle, delimiter="\t")}


def test_summarize_variant_bases_matches_independent_base_counts(run_workflow, vcf_inputs, source):
    outputs = run_workflow(
        f"{UTILS}/SummarizeVariantBases.wdl",
        {
            "vcfs": [vcf_inputs["vcf"]],
            "vcf_idxs": [vcf_inputs["vcf_idx"]],
            "vcf_names": ["subset"],
            "contig": CONTIG,
            "prefix": "bases",
            "ref_fai": str(FAI),
            "shard_bin_size": 15000000,
            "utils_docker": vcf_inputs["utils_docker"],
        },
    )

    def size_bin(length):
        return "1-49" if length < 50 else "50-499" if length < 500 else "500+"

    altered = collections.defaultdict(set)
    allele_bases = collections.defaultdict(int)
    for r in source[1]:
        allele_type, length = r["info"]["allele_type"], abs(r["info"].get("allele_length", 0))
        if allele_type == "snv":
            category, length = "SNV", 1
        elif allele_type in ("del", "ins") or "dup" in allele_type:
            category = f"{'DEL' if allele_type == 'del' else 'INS'} {size_bin(length)}"
        else:
            continue
        ac = sum(r["info"]["AC"])
        if ac == 0:
            continue
        allele_bases[category] += ac * length
        if allele_type == "snv":
            altered[category].add(r["pos"])
        elif allele_type == "del":
            altered[category].update(range(r["pos"] + 1, r["pos"] + length + 1))
    assert altered["DEL 1-49"] and altered["DEL 50-499"] and allele_bases["INS 50-499"]

    site = read_base_table(outputs["SummarizeVariantBases.site_bases_tsv"])
    sample = read_base_table(outputs["SummarizeVariantBases.sample_bases_tsv"])
    for category in set(allele_bases) | set(altered):
        assert float(site[category]["subset_bases"]) == len(altered[category]), category
        assert float(sample[category]["subset_bases"]) == pytest.approx(allele_bases[category] / len(SAMPLES), abs=1e-6)
    assert float(site["Total"]["subset_bases"]) == sum(len(v) for v in altered.values())
    assert float(site["Total"]["subset_proportion"]) == pytest.approx(
        sum(len(v) for v in altered.values()) / CONTIG_LENGTH, abs=1e-6
    )
    assert float(sample["Total"]["subset_bases"]) == pytest.approx(sum(allele_bases.values()) / len(SAMPLES), abs=1e-5)
    assert {row["contig_bases"] for row in site.values()} == {str(CONTIG_LENGTH)}
    assert all(float(site[f"INS {b}"]["subset_bases"]) == 0 for b in ("1-49", "50-499", "500+"))


def read_pairs(path):
    with open(path) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        return reader.fieldnames, list(reader)


def test_extract_random_calls_is_reproducible_and_respects_allele_type_filter(run_workflow, vcf_inputs, source):
    inputs = {
        "vcfs": [vcf_inputs["vcf"]],
        "vcf_idxs": [vcf_inputs["vcf_idx"]],
        "prefix": "random",
        "count": 12,
        "random_seed": 7,
        "allele_types": ["del"],
        "utils_docker": vcf_inputs["utils_docker"],
    }
    outputs = run_workflow(f"{UTILS}/ExtractRandomCalls.wdl", inputs)

    _, rows = read_pairs(outputs["ExtractRandomCalls.variant_sample_pairs"])
    by_id = {r["id"]: r for r in source[1]}
    assert len(rows) == 12
    assert len({(row["ID"], row["samples"]) for row in rows}) == 12
    for row in rows:
        record = by_id[row["ID"]]
        assert record["info"]["allele_type"] == "del" == row["allele_type"]
        assert is_carrier(record, row["samples"])

    candidates = sum(is_carrier(r, s) for r in source[1] if r["info"]["allele_type"] == "del" for s in SAMPLES)
    summary = dict(line.split("\t") for line in open(outputs["ExtractRandomCalls.candidate_summary"]).read().splitlines())
    assert summary == {"requested": "12", "total_candidates_found": str(candidates), "pairs_written": "12"}

    _, selected_vcf = load_records(outputs["ExtractRandomCalls.variant_vcf"])
    assert {r["id"] for r in selected_vcf} == {row["ID"] for row in rows}

    first = open(outputs["ExtractRandomCalls.variant_sample_pairs"]).read()
    again = run_workflow(f"{UTILS}/ExtractRandomCalls.wdl", inputs)
    assert open(again["ExtractRandomCalls.variant_sample_pairs"]).read() == first
