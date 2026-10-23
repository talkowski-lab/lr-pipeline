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


def rewrite_vcf(src, dst, drop_info, edit):
    """Write a VCF without the drop_info INFO fields, applying edit() to each record, and return its bgzipped, indexed path."""
    plain = dst.with_suffix("")
    with pysam.VariantFile(str(src)) as vcf_in:
        header = vcf_in.header.copy()
        for key in drop_info:
            header.info.remove_header(key)
        with pysam.VariantFile(str(plain), "w", header=header) as vcf_out:
            for record in vcf_in:
                edit(record)
                for key in drop_info:
                    if key in record.info:
                        del record.info[key]
                record.translate(header)
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


@pytest.mark.parametrize("sample", ["HG00096", "HG00126"])
def test_extract_sample_vcfs_routes_carried_variants_by_length_and_source(run_workflow, vcf_inputs, source, sample):
    min_sv_length = 50
    outputs = run_workflow(
        f"{UTILS}/ExtractSampleVcfs.wdl",
        {
            "cohort_vcf": vcf_inputs["vcf"],
            "cohort_vcf_idx": vcf_inputs["vcf_idx"],
            "contigs": [CONTIG],
            "prefix": sample,
            "sample_id": sample,
            "min_sv_length": min_sv_length,
            "utils_docker": vcf_inputs["utils_docker"],
        },
    )

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

    snv_samples, snv_records = load_records(outputs["ExtractSampleVcfs.snv_indel_vcf"])
    sv_samples, sv_records = load_records(outputs["ExtractSampleVcfs.sv_vcf"])
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


@pytest.mark.parametrize("compute_allele_counts", [True, False])
def test_postprocess_callset_recomputes_allele_counts_for_singleton_filter(
    run_workflow, run_dir, subset_vcf, utils_docker, source, compute_allele_counts
):
    stripped = rewrite_vcf(subset_vcf, run_dir / "no_ac.vcf.gz", ["AC", "AN"], lambda r: r.filter.clear())
    run_flags = [
        "clean_vcf_header",
        "decrement_trv_ids",
        "drop_filters",
        "filter_single_read_singletons",
        "flag_homopolymer_trvs",
        "normalize_ploidy",
        "prune_meis",
        "reassign_suffixes",
        "sorting",
        "transfer_genotypes",
        "unphase_samples",
    ]
    inputs = {
        "vcf": str(stripped),
        "vcf_idx": f"{stripped}.tbi",
        "contigs": [CONTIG],
        "prefix": "post",
        **{f"run_{flag}": False for flag in run_flags},
        "run_compute_allele_counts": compute_allele_counts,
        "run_filter_assembly_only_singletons": True,
        "utils_docker": utils_docker,
    }
    # Without an AC header definition the assembly-only singleton filter cannot read AC and fails
    if not compute_allele_counts:
        with pytest.raises(WorkflowError):
            run_workflow(f"{UTILS}/PostprocessCallset.wdl", inputs)
        return

    outputs = run_workflow(f"{UTILS}/PostprocessCallset.wdl", inputs)
    _, records = load_records(outputs["PostprocessCallset.post_processed_vcf"])
    assert [site_key(r) for r in records] == [site_key(r) for r in source[1]]
    # The fixture's AC and AN were recomputed from these genotypes when it was subset
    assert [r["info"]["AC"] for r in records] == [r["info"]["AC"] for r in source[1]]
    assert [r["info"]["AN"] for r in records] == [r["info"]["AN"] for r in source[1]]
    flagged = {r["id"] for r in records if "ASSEMBLY_ONLY_SINGLETON" in r["filter"]}
    expected = {r["id"] for r in source[1] if "ASSEMBLY_ONLY_SINGLETON" in r["filter"]}
    assert expected and expected <= flagged


def test_normalize_allele_types_reproduces_released_types(run_workflow, run_dir, subset_vcf, utils_docker, source):
    def restore_raw_type(record):
        subtype = record.info.get("allele_subtype")
        if subtype:
            record.info["allele_type"] = "dup" if subtype == "tandem_dup" else subtype

    raw = rewrite_vcf(subset_vcf, run_dir / "raw_types.vcf.gz", ["allele_subtype"], restore_raw_type)
    outputs = run_workflow(
        f"{UTILS}/NormalizeAlleleTypes.wdl",
        {"vcf": str(raw), "vcf_idx": f"{raw}.tbi", "prefix": "normalized", "utils_docker": utils_docker},
    )

    _, records = load_records(outputs["NormalizeAlleleTypes.transformed_vcf"])
    assert [site_key(r) for r in records] == [site_key(r) for r in source[1]]

    def types(r):
        return r["info"]["allele_type"], r["info"].get("allele_subtype")

    assert [types(r) for r in records] == [types(r) for r in source[1]]
    assert {"dup_interspersed", "alu_ins", "alu_del"} <= {r["info"].get("allele_subtype") for r in records}


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


def test_extract_random_calls_keeps_records_sharing_a_span(run_workflow, vcf_inputs, source):
    # The subset holds distinct insertions with identical POS and END, which must each reach the VCF
    inputs = {
        "vcfs": [vcf_inputs["vcf"]],
        "vcf_idxs": [vcf_inputs["vcf_idx"]],
        "prefix": "random",
        "count": 10000,
        "allele_types": ["ins"],
        "utils_docker": vcf_inputs["utils_docker"],
    }
    outputs = run_workflow(f"{UTILS}/ExtractRandomCalls.wdl", inputs)

    _, rows = read_pairs(outputs["ExtractRandomCalls.variant_sample_pairs"])
    carried = {r["id"] for r in source[1] if r["info"]["allele_type"] == "ins" and any(is_carrier(r, s) for s in SAMPLES)}
    assert {row["ID"] for row in rows} == carried

    _, selected_vcf = load_records(outputs["ExtractRandomCalls.variant_vcf"])
    assert sorted(r["id"] for r in selected_vcf) == sorted(carried)
