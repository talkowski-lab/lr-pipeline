import collections
import csv
import gzip
import itertools
import math
import random
import statistics
import subprocess
import sys
import tarfile
from pathlib import Path

import edlib
import pysam
import pytest
from conftest import (
    CONTIG,
    CONTIG_LENGTH,
    FAI,
    REPO_ROOT,
    SAMPLES,
    UTILS,
    WINDOW_FASTA,
    WorkflowError,
    index_table,
    is_carrier,
    load_records,
    read_dicts,
    read_table,
    rewrite_vcf,
    site_key,
    write_mosdepth_bed,
    write_table,
    write_vcf,
)
from pysam import bcftools

pytestmark = pytest.mark.wdl


def key_columns(record):
    return [record["contig"], str(record["pos"]), record["ref"], ",".join(record["alts"]), record["id"] or "."]


def histogram(values, bins):
    """Count values into [bins[i], bins[i+1]) intervals, returning the counts and the number at or above the last edge."""
    counts, larger = [0] * (len(bins) - 1), 0
    for value in values:
        if value < bins[0]:
            continue
        if value >= bins[-1]:
            larger += 1
            continue
        counts[max(i for i in range(len(bins) - 1) if bins[i] <= value)] += 1
    return counts, larger


def het_allele_balance(record):
    balances = []
    for sample in SAMPLES:
        ad = record["format"][sample].get("AD")
        if None in record["gts"][sample] or sorted(record["gts"][sample]) != [0, 1] or ad is None or None in ad[:2]:
            continue
        if ad[0] + ad[1] > 0:
            balances.append(ad[0] / (ad[0] + ad[1]))
    return balances


def format_values(record, field):
    return [(record["gts"][s], record["format"][s].get(field)) for s in SAMPLES]


def scalar(value):
    return value[0] if isinstance(value, tuple) else value


def test_annotate_gq_metrics_bins_format_values_and_allele_balance(run_workflow, vcf_inputs, source):
    gq_bins, dp_bins, ab_bins = [0, 20, 40, 60], [0, 10, 20, 30], [0, 0.25, 0.5, 0.75, 1.0]
    outputs = run_workflow(
        "wdl/annotation/AnnotateGQMetrics.wdl",
        {
            **vcf_inputs,
            "contigs": [CONTIG],
            "prefix": "gq",
            "records_per_shard": 200,
            "gq_fields": ["GQ", "DP"],
            "gq_bins": [gq_bins, dp_bins],
            "gq_variant_filters": ["", ""],
            "gq_larger_field": [True, False],
            "ab_annotation": True,
            "ab_bins": ab_bins,
        },
    )

    rows = read_table(outputs["AnnotateGQMetrics.annotations_tsv_gq"])
    assert len(rows) == len(source[1])
    for row, record in zip(rows, source[1]):
        assert row[:5] == key_columns(record)
        values = [(gt, v) for gt, v in format_values(record, "GQ") if v is not None]
        all_counts, all_larger = histogram([v for _, v in values], gq_bins)
        alt_counts, alt_larger = histogram([v for gt, v in values if any(a for a in gt if a is not None and a > 0)], gq_bins)
        assert row[5] == "|".join(map(str, all_counts))
        assert row[6] == "|".join(map(str, alt_counts))
        assert row[7:9] == [str(all_larger), str(alt_larger)]
        assert row[-1] == "|".join(map(str, histogram(het_allele_balance(record), ab_bins)[0]))
    assert any(int(row[7]) > 0 for row in rows)


def hwe_p(n_ref, n_het, n_alt):
    n = n_ref + n_het + n_alt
    p = (2 * n_ref + n_het) / (2 * n)
    expected = [p * p * n, 2 * p * (1 - p) * n, (1 - p) * (1 - p) * n]
    chisq = sum((o - e) ** 2 / e for o, e in zip([n_ref, n_het, n_alt], expected) if e > 0)
    return round(math.erfc(math.sqrt(chisq / 2)), 6)


def test_annotate_sq_metrics_computes_site_quality_from_genotype_fields(run_workflow, vcf_inputs, source):
    outputs = run_workflow(
        "wdl/annotation/AnnotateSQMetrics.wdl", {**vcf_inputs, "contigs": [CONTIG], "prefix": "sq", "records_per_shard": 200}
    )

    rows = read_table(outputs["AnnotateSQMetrics.annotations_tsv_sq"])
    assert len(rows) == len(source[1])
    checked = 0
    for row, record in zip(rows, source[1]):
        assert row[:5] == key_columns(record)
        called = [s for s in SAMPLES if None not in record["gts"][s]]
        n_non_ref = {s: sum(a > 0 for a in record["gts"][s]) for s in called}
        counts = [sum(v == k for v in n_non_ref.values()) for k in (0, 1)]
        counts.append(len(called) - sum(counts))
        assert row[10] == (str(hwe_p(*counts)) if called else ".")
        if len(record["alts"]) != 1:
            continue
        # Biallelic: depth and phred-scaled hom-ref likelihood are summed over the carriers of the alt allele
        carriers = [s for s in called if record["gts"][s].count(1) > 0]
        depth = [record["format"][s].get("DP") for s in carriers]
        pl = [record["format"][s].get("PL") for s in carriers]
        if carriers and None not in depth:
            assert row[9] == str(sum(depth))
        if carriers and all(p is not None and p[0] is not None for p in pl):
            assert row[7] == str(sum(p[0] for p in pl))
            checked += 1
        n_het = sum(record["gts"][s].count(1) == 1 for s in called)
        n_hom = sum(record["gts"][s].count(1) == 2 for s in called)
        n_ref = len(called) - n_het - n_hom
        if called:
            p = (2 * n_ref + n_het) / (2 * len(called))
            expected_het = 2 * p * (1 - p) * len(called)
            assert row[5] == (str(round(1 - n_het / expected_het, 6)) if expected_het > 0 else "0.0")
    assert checked > 100


def test_annotate_vcf_applies_filtered_and_column_subset_tsvs_in_order(run_workflow, run_dir, vcf_inputs, source):
    records = source[1]
    key = lambda r: [r["contig"], r["pos"], r["ref"], ",".join(r["alts"]), r["id"]]  # noqa: E731
    first = write_table(
        run_dir / "first.tsv",
        [key(r) + [r["pos"] % 100, f"label_{r['pos']}", "unused"] for r in reversed(records)],
    )
    # The awk condition keeps only even positions of the second table
    second = write_table(run_dir / "second.tsv", [key(r) + [r["pos"] % 7 + 0.5] for r in records])
    outputs = run_workflow(
        f"{UTILS}/AnnotateVcf.wdl",
        {
            **vcf_inputs,
            "prefix": "annotated",
            "contigs": [CONTIG],
            "records_per_shard": 100,
            "annotations_tsvs": [str(first), str(second)],
            "sort_tsvs": [True, False],
            "subset_tsv_columns": [[6, 7], []],
            "awk_tsv_conditions": ["", "$2 % 2 == 0"],
            "include_expressions": ["", 'INFO/allele_type="snv"'],
            "strip_info_fields_per_tsv": [False, True],
            "info_names": [["SCORE", "LABEL"], ["phylop"]],
            "info_descriptions": [["Score", "Label"], ["Replaced phylop"]],
            "info_types": [["Integer", "String"], ["Float"]],
            "info_numbers": [["1", "1"], ["1"]],
        },
    )

    # The second table annotates only SNVs, and the stripped field is restored only where an SNV has a row
    _, got = load_records(outputs["AnnotateVcf.annotated_vcf"])
    assert [site_key(r) for r in got] == [site_key(r) for r in records]
    for new, old in zip(got, records):
        assert new["info"]["SCORE"] == old["pos"] % 100
        assert new["info"]["LABEL"] == f"label_{old['pos']}"
        assert "unused" not in new["info"]
        if old["info"]["allele_type"] == "snv" and old["pos"] % 2 == 0:
            assert new["info"]["phylop"] == pytest.approx(old["pos"] % 7 + 0.5)
        else:
            assert "phylop" not in new["info"]
        assert new["gts"] == old["gts"]
    assert any("phylop" in r["info"] for r in got) and any("phylop" not in r["info"] for r in got)
    assert any(r["info"]["allele_type"] != "snv" for r in got)


def trv_records(records, source_name):
    return [r for r in records if r["info"]["allele_type"] == "trv" and r["info"]["SOURCE"] == source_name]


def overlap_ends(record):
    return record["pos"], record["pos"] + len(record["ref"])


def test_combine_trs_dedups_each_caller_then_adds_non_overlapping_lower_priority_loci(
    run_workflow, run_dir, subset_vcf, utils_docker, source
):
    records = source[1]
    vcfs = [
        write_vcf(
            subset_vcf,
            run_dir / f"{name}.vcf.gz",
            lambda r, name=name: r.info["SOURCE"] == name and r.info["allele_type"] == "trv",
        )
        for name in ("TRExplorer", "Vamos")
    ]
    outputs = run_workflow(
        f"{UTILS}/CombineTRs.wdl",
        {
            "tr_vcfs": [str(v) for v in vcfs],
            "tr_vcf_idxs": [f"{v}.tbi" for v in vcfs],
            "tr_callers": ["TRExplorer", "Vamos"],
            "contigs": [CONTIG],
            "sample_ids": SAMPLES,
            "prefix": "combined",
            "utils_docker": utils_docker,
        },
    )

    def dedup(group):
        ids = {}
        for r in group:
            ids.setdefault(f"{r['contig']}-{r['pos']}-TRV-{len(r['ref']) - 1}", []).append(r)
        named = []
        for new_id, same in ids.items():
            for n, r in enumerate(same, 1):
                named.append((new_id if len(same) == 1 else f"{new_id}_{n}", r))
        keep = []
        for rid, r in named:
            start, end = overlap_ends(r)
            beaten = any(
                other_id != rid
                and overlap_ends(o)[0] < end
                and start < overlap_ends(o)[1]
                and (len(r["ref"]) > len(o["ref"]) or (len(r["ref"]) == len(o["ref"]) and rid > other_id))
                for other_id, o in named
            )
            if not beaten:
                keep.append((rid, r))
        return keep

    expected = dedup(trv_records(records, "TRExplorer"))
    merged_spans = [overlap_ends(r) for _, r in expected]
    expected += [
        (rid, r)
        for rid, r in dedup(trv_records(records, "Vamos"))
        if not any(s < overlap_ends(r)[1] and overlap_ends(r)[0] < e for s, e in merged_spans)
    ]
    expected.sort(key=lambda pair: pair[1]["pos"])

    samples, got = load_records(outputs["CombineTRs.combined_tr_vcf"])
    assert samples == SAMPLES
    assert sorted((r["pos"], r["id"], r["info"]["SOURCE"]) for r in got) == sorted(
        (r["pos"], rid, r["info"]["SOURCE"]) for rid, r in expected
    )
    assert {r["info"]["SOURCE"] for r in got} == {"TRExplorer", "Vamos"}


def test_subset_tsv_to_columns_keeps_named_columns_and_filters_rows(run_workflow, run_dir, utils_docker, source):
    records = source[1]
    names = ["alpha", "beta", "gamma"]
    rows = [
        [r["contig"], r["pos"], r["ref"], ",".join(r["alts"]), r["id"], r["pos"] % 3, f"b{r['pos'] % 5}", "g"] for r in records
    ]
    tsv = write_table(run_dir / "annotations.tsv", rows)
    header = run_dir / "header.txt"
    header.write_text("\n".join(names) + "\n")

    outputs = run_workflow(
        f"{UTILS}/SubsetTsvToColumns.wdl",
        {
            "annotations_tsv": str(tsv),
            "annotations_header": str(header),
            "subset_columns": ["gamma", "alpha"],
            "subset_column_values": [[], ["0", "2"]],
            "prefix": "subset",
            "utils_docker": utils_docker,
        },
    )

    expected = [[str(v) for v in row[:5]] + ["g", str(row[5])] for row in rows if row[5] in (0, 2)]
    assert 0 < len(expected) < len(rows)
    assert read_table(outputs["SubsetTsvToColumns.subset_tsv"]) == expected


def test_create_depth_intervals_tiles_each_contig_in_whole_bins(run_workflow, utils_docker):
    bin_size = 10000000
    outputs = run_workflow(
        f"{UTILS}/CreateDepthIntervals.wdl",
        {"ref_fai": str(FAI), "contigs": [CONTIG], "prefix": "bins", "bin_size": bin_size, "utils_docker": utils_docker},
    )

    expected = [f"{CONTIG}:{end - bin_size + 1}-{end}" for end in range(bin_size, CONTIG_LENGTH + 1, bin_size)]
    assert len(expected) == 5
    assert [line[0] for line in read_table(outputs["CreateDepthIntervals.intervals"])] == expected


def test_create_depth_intervals_rejects_contigs_missing_from_index(run_workflow, utils_docker):
    with pytest.raises(WorkflowError):
        run_workflow(
            f"{UTILS}/CreateDepthIntervals.wdl",
            {"ref_fai": str(FAI), "contigs": ["chrNotThere"], "prefix": "bins", "bin_size": 1000, "utils_docker": utils_docker},
        )


@pytest.fixture
def halves(run_dir, subset_vcf, source):
    """Split the subset VCF at its median position into two indexed VCFs."""
    split_at = source[1][len(source[1]) // 2]["pos"]
    low = write_vcf(subset_vcf, run_dir / "low.vcf.gz", lambda r: r.pos < split_at)
    high = write_vcf(subset_vcf, run_dir / "high.vcf.gz", lambda r: r.pos >= split_at)
    return [low, high]


def test_merge_trv_length_deltas_concatenates_alleles_and_sums_length_counts(run_workflow, halves, utils_docker):
    parts = [
        run_workflow(
            f"{UTILS}/ExtractTRVLengthDeltas.wdl",
            {"vcf": str(v), "vcf_idx": f"{v}.tbi", "prefix": v.name[:3], "utils_docker": utils_docker},
        )
        for v in halves
    ]
    allele_tsvs = [p["ExtractTRVLengthDeltas.allele_deltas_tsv"] for p in parts]
    length_tsvs = [p["ExtractTRVLengthDeltas.length_deltas_tsv"] for p in parts]

    outputs = run_workflow(
        f"{UTILS}/MergeTRVLengthDeltas.wdl",
        {"allele_deltas_tsvs": allele_tsvs, "length_deltas_tsvs": length_tsvs, "prefix": "merged", "utils_docker": utils_docker},
    )

    allele_rows = [read_table(path) for path in allele_tsvs]
    merged = read_table(outputs["MergeTRVLengthDeltas.merged_allele_deltas_tsv"])
    assert merged[0] == allele_rows[0][0]
    assert merged[1:] == allele_rows[0][1:] + allele_rows[1][1:]

    n_alleles, allele_count = collections.defaultdict(int), collections.defaultdict(int)
    for path in length_tsvs:
        for row in read_dicts(path):
            n_alleles[int(row["length_delta"])] += int(row["n_alleles"])
            allele_count[int(row["length_delta"])] += int(row["allele_count"])
    merged_counts = read_table(outputs["MergeTRVLengthDeltas.merged_length_deltas_tsv"])
    assert merged_counts[0] == ["length_delta", "n_alleles", "allele_count"]
    assert merged_counts[1:] == [[str(d), str(n_alleles[d]), str(allele_count[d])] for d in sorted(n_alleles)]
    assert len(n_alleles) > 1


def test_merge_variant_base_tables_sums_bases_and_recomputes_proportions(run_workflow, halves, utils_docker):
    parts = [
        run_workflow(
            f"{UTILS}/SummarizeVariantBases.wdl",
            {
                "vcfs": [str(v)],
                "vcf_idxs": [f"{v}.tbi"],
                "vcf_names": ["subset"],
                "contig": CONTIG,
                "prefix": v.name[:3],
                "ref_fai": str(FAI),
                "shard_bin_size": 25000000,
                "utils_docker": utils_docker,
            },
        )
        for v in halves
    ]
    outputs = run_workflow(
        f"{UTILS}/MergeVariantBaseTables.wdl",
        {
            "site_bases_tsvs": [p["SummarizeVariantBases.site_bases_tsv"] for p in parts],
            "sample_bases_tsvs": [p["SummarizeVariantBases.sample_bases_tsv"] for p in parts],
            "prefix": "merged",
            "utils_docker": utils_docker,
        },
    )

    for kind in ("site", "sample"):
        tables = [{r["category"]: r for r in read_dicts(p[f"SummarizeVariantBases.{kind}_bases_tsv"])} for p in parts]
        merged = {r["category"]: r for r in read_dicts(outputs[f"MergeVariantBaseTables.{kind}_bases_tsv"])}
        assert list(merged) == list(tables[0])
        for category, row in merged.items():
            assert int(row["contig_bases"]) == 2 * CONTIG_LENGTH
            assert float(row["subset_bases"]) == pytest.approx(sum(float(t[category]["subset_bases"]) for t in tables), abs=1e-5)
            assert float(row["subset_proportion"]) == pytest.approx(
                float(row["subset_bases"]) / int(row["contig_bases"]), abs=1e-6
            )
        assert float(merged["Total"]["subset_bases"]) > 0


def trv_vcf(src, dst, source_name=None, sample=None):
    """Write the tandem-repeat records of a VCF, optionally of one source and one sample, as an indexed VCF."""
    plain = dst.with_suffix("")
    with pysam.VariantFile(str(src)) as vcf_in:
        if sample:
            vcf_in.subset_samples([sample])
        with pysam.VariantFile(str(plain), "w", header=vcf_in.header) as vcf_out:
            for record in vcf_in:
                if record.info["allele_type"] == "trv" and source_name in (None, record.info["SOURCE"]):
                    vcf_out.write(record)
    pysam.tabix_index(str(plain), preset="vcf", force=True)
    return dst


def test_extract_disparate_tr_loci_lists_loci_missing_from_or_overlapping_the_other_callset(
    run_workflow, run_dir, subset_vcf, utils_docker, source
):
    a = trv_vcf(subset_vcf, run_dir / "a.vcf.gz", "TRExplorer")
    counter = itertools.count()

    def shift_every_fifth(record):
        if next(counter) % 5 == 0:
            record.pos += 3

    # The second callset repeats the first with a fifth of its loci shifted by three bases, so those overlap without matching
    b = rewrite_vcf(a, run_dir / "b.vcf.gz", [], shift_every_fifth)
    outputs = run_workflow(
        f"{UTILS}/ExtractDisparateTRLoci.wdl",
        {
            "vcf_a": str(a),
            "vcf_a_idx": f"{a}.tbi",
            "vcf_b": str(b),
            "vcf_b_idx": f"{b}.tbi",
            "contig": CONTIG,
            "prefix": "disparate",
            "utils_docker": utils_docker,
        },
    )

    def loci(records):
        return [(r["pos"] - 1, r["pos"] - 1 + len(r["ref"]), r["pos"], len(r["ref"]), r["info"]["TRID"]) for r in records]

    loci_a, loci_b = loci(trv_records(source[1], "TRExplorer")), loci(load_records(b)[1])
    expected_missing = {
        ("vcf_a", "vcf_b", CONTIG, str(p), str(n)) for s, e, p, n, _ in loci_a if not any((s, e) == (x, y) for x, y, *_ in loci_b)
    }
    expected_missing |= {
        ("vcf_b", "vcf_a", CONTIG, str(p), str(n)) for s, e, p, n, _ in loci_b if not any((s, e) == (x, y) for x, y, *_ in loci_a)
    }
    expected_overlap = {
        (CONTIG, str(pa), str(na), ta, CONTIG, str(pb), str(nb), tb)
        for sa, ea, pa, na, ta in loci_a
        for sb, eb, pb, nb, tb in loci_b
        if sa < eb and sb < ea and (pa, na) != (pb, nb)
    }

    missing = read_table(outputs["ExtractDisparateTRLoci.missing_variants_tsv"])
    assert missing[0] == ["present_in", "missing_from", "chrom", "pos", "ref_length"]
    assert {tuple(r) for r in missing[1:]} == expected_missing
    overlapping = read_table(outputs["ExtractDisparateTRLoci.overlapping_variants_tsv"])
    assert overlapping[0] == ["chrom_a", "pos_a", "ref_length_a", "trid_a", "chrom_b", "pos_b", "ref_length_b", "trid_b"]
    assert {tuple(r) for r in overlapping[1:]} == expected_overlap
    assert expected_missing and expected_overlap


def test_evaluate_overlapping_tr_loci_compares_haplotypes_over_the_union_span(
    run_workflow, run_dir, subset_vcf, utils_docker, source
):
    sample = "HG00096"
    trvs = trv_vcf(subset_vcf, run_dir / "trvs.vcf.gz", sample=sample)
    # Add a copy of every second non-reference locus shifted by three bases, so that it overlaps the locus it was copied from
    plain = run_dir / "single.vcf"
    with pysam.VariantFile(str(trvs)) as vcf_in:
        records = list(vcf_in)
        copies = []
        for record in records:
            if any(a for a in record.samples[0]["GT"] if a):
                copy = record.copy()
                copy.pos += 3
                copy.info["TRID"] = f"{record.info['TRID']}_shifted"
                copies.append(copy)
        with pysam.VariantFile(str(plain), "w", header=vcf_in.header) as vcf_out:
            for record in sorted(records + copies[::2], key=lambda r: r.pos):
                vcf_out.write(record)
    pysam.tabix_index(str(plain), preset="vcf", force=True)
    vcf = run_dir / "single.vcf.gz"
    outputs = run_workflow(
        f"{UTILS}/EvaluateOverlappingTRLoci.wdl",
        {"vcf": str(vcf), "vcf_idx": f"{vcf}.tbi", "contigs": [CONTIG], "prefix": "overlap", "utils_docker": utils_docker},
    )

    loci = []
    for r in load_records(vcf)[1]:
        gt = r["gts"][sample]
        if None in gt:
            continue
        alleles = [r["ref"].upper()] + [a.upper() for a in r["alts"]]
        loci.append(
            {
                "trid": r["info"]["TRID"],
                "start": r["pos"],
                "end": r["pos"] + len(r["ref"]) - 1,
                "ref": r["ref"].upper(),
                "haps": [alleles[a] for a in gt],
                "nonref": any(a > 0 for a in gt),
            }
        )

    def padded(own, other, union_start, union_end):
        left = other["ref"][: max(own["start"] - other["start"], 0)]
        right = other["ref"][own["end"] + 1 - other["start"] : union_end - other["start"] + 1] if own["end"] < union_end else ""
        assert union_start <= own["start"]
        return [left + hap + right for hap in own["haps"]]

    expected = {}
    for i, b in enumerate(loci):
        for a in loci[:i]:
            if a["end"] < b["start"] or not (a["nonref"] or b["nonref"]):
                continue
            union_start, union_end = min(a["start"], b["start"]), max(a["end"], b["end"])
            seqs_a, seqs_b = padded(a, b, union_start, union_end), padded(b, a, union_start, union_end)
            best = min(
                sum(edlib.align(seqs_a[i], seqs_b[j], task="distance")["editDistance"] for i, j in pairing)
                for pairing in ([(0, 0), (1, 1)], [(0, 1), (1, 0)])
            )
            expected[(a["trid"], b["trid"])] = (
                [str(a["start"]), str(a["end"]), str(b["start"]), str(b["end"])],
                [str(union_start), str(union_end), str(union_end - union_start + 1)],
                best,
            )
    assert expected

    rows = read_dicts(outputs["EvaluateOverlappingTRLoci.overlapping_loci_tsv"])
    assert {(r["trid_a"], r["trid_b"]) for r in rows} == set(expected)
    for row in rows:
        spans, union, distance = expected[(row["trid_a"], row["trid_b"])]
        assert [row["locus_start_a"], row["locus_end_a"], row["locus_start_b"], row["locus_end_b"]] == spans
        assert [row["union_start"], row["union_end"], row["union_length"]] == union
        assert int(row["min_edit_distance"]) == distance
        assert 0 <= float(row["max_similarity"]) <= 1


def size_range(length):
    for limit, label in ((50, "<50"), (200, "50-199"), (500, "200-499"), (5000, "500-4999")):
        if length < limit:
            return label
    return "5000+"


def test_summarize_singleton_calls_counts_genotypes_by_stratum(run_workflow, vcf_inputs, source):
    min_length = 50
    outputs = run_workflow(
        f"{UTILS}/SummarizeSingletonCalls.wdl",
        {
            "vcfs": [vcf_inputs["vcf"]],
            "vcf_idxs": [vcf_inputs["vcf_idx"]],
            "prefix": "singletons",
            "min_length": min_length,
            "records_per_shard": 100,
            "utils_docker": vcf_inputs["utils_docker"],
        },
    )

    counts = collections.defaultdict(int)
    for r in source[1]:
        if "allele_length" not in r["info"] or abs(r["info"]["allele_length"]) < min_length:
            continue
        for sample in SAMPLES:
            if None in r["gts"][sample]:
                continue
            alt_count = sum(a > 0 for a in r["gts"][sample])
            ev = ",".join(r["format"][sample].get("EV") or ["."])
            callers = set() if ev in ("", ".") else {v.split("_(", 1)[0].lower() for v in ev.split(",")}
            count_type = "alignments" if callers - {"hapdiff", "dipcall"} else "assemblies" if callers else "other"
            singleton = "ac_0" if alt_count == 0 else "ac_1" if alt_count == 1 else "ac_>1"
            if singleton == "ac_0":
                count_type = "other"
            key = (
                r["info"]["allele_type"],
                size_range(abs(r["info"]["allele_length"])),
                r["info"].get("REGION", "."),
                count_type,
                singleton,
            )
            counts[(sample, *key)] += 1
    assert any(k[-1] == "ac_1" and k[-2] != "other" for k in counts)

    rows = read_table(outputs["SummarizeSingletonCalls.singleton_counts_tsv"])
    header, body = rows[0], rows[1:]
    assert header[0] == "sample"
    assert [row[0] for row in body] == sorted(SAMPLES)
    for row in body:
        got = {tuple(h.split(" - ")): int(v) for h, v in zip(header[1:], row[1:])}
        expected = {k[1:]: v for k, v in counts.items() if k[0] == row[0]}
        assert {k: v for k, v in got.items() if v} == expected


def test_resolve_haplotype_overlaps_clears_the_losing_call_of_each_overlapping_pair(
    run_workflow, run_dir, subset_vcf, utils_docker
):
    # Add a copy of each tenth phased, non-repeat deletion shifted by one base, so it overlaps its original on the same haplotypes
    plain = run_dir / "overlapping.vcf"
    with pysam.VariantFile(str(subset_vcf)) as vcf_in:
        records, copies = list(vcf_in), []
        for record in records:
            phased = any(record.samples[s].phased and 1 in record.samples[s]["GT"] for s in SAMPLES)
            if record.info["allele_type"] == "del" and "TR_ENVELOPED" not in record.info and phased:
                if len(copies) < 30 and record.pos % 10 < 5:
                    copy = record.copy()
                    copy.pos += 1
                    copy.id = f"{record.id}_copy"
                    copies.append(copy)
        with pysam.VariantFile(str(plain), "w", header=vcf_in.header) as vcf_out:
            for record in sorted(records + copies, key=lambda r: r.pos):
                vcf_out.write(record)
    pysam.tabix_index(str(plain), preset="vcf", force=True)
    vcf = run_dir / "overlapping.vcf.gz"
    assert copies

    outputs = run_workflow(
        f"{UTILS}/ResolveHaplotypeOverlaps.wdl",
        {
            "vcf": str(vcf),
            "vcf_idx": f"{vcf}.tbi",
            "prefix": "resolved",
            "records_per_shard": 150,
            "utils_docker": utils_docker,
        },
    )

    _, before = load_records(vcf)
    _, after = load_records(outputs["ResolveHaplotypeOverlaps.overlap_resolved_vcf"])
    by_id = {r["id"]: r for r in before}
    overlaps = read_dicts(outputs["ResolveHaplotypeOverlaps.overlap_tsv"])
    assert len(overlaps) >= len(copies)

    def span(record):
        return record["pos"], record["pos"] + len(record["ref"]) - 1

    cleared = set()
    for row in overlaps:
        hap = int(row["haplotype"]) - 1
        retained, lost = by_id[row["variant_id_retained"]], by_id[row["variant_id_cleared"]]
        for record in (retained, lost):
            assert record["gts"][row["sample"]][hap] > 0
        assert span(retained)[0] <= span(lost)[1] and span(lost)[0] <= span(retained)[1]
        assert len(retained["ref"]) >= len(lost["ref"])
        cleared.add((row["variant_id_cleared"], row["sample"]))

    assert [site_key(r) for r in after] == [site_key(r) for r in before]
    for old, new in zip(before, after):
        for sample in SAMPLES:
            if (old["id"], sample) in cleared:
                assert new["gts"][sample] == (None, None)
            else:
                assert new["gts"][sample] == old["gts"][sample]

    # No two phased calls left on the same haplotype overlap
    for sample in SAMPLES:
        for hap in (0, 1):
            spans = sorted(
                span(r)
                for r in after
                if r["phased"][sample]
                and (r["gts"][sample][hap] or 0) > 0
                and r["info"]["allele_type"] != "trv"
                and "TR_ENVELOPED" not in r["info"]
            )
            reach = 0
            for start, end in spans:
                assert start > reach
                reach = max(reach, end)


def test_fill_phased_genotypes_restores_missing_genotypes_and_format_fields_from_unphased_vcf(
    run_workflow, run_dir, subset_vcf, utils_docker, source
):
    blanked, depth_reset = "HG00096", "HG00097"

    def degrade(record):
        record.samples[blanked]["GT"] = (None, None)
        record.samples[depth_reset]["DP"] = 1
        if record.samples[depth_reset].get("AL") is not None:
            record.samples[depth_reset]["AL"] = "0"

    phased = rewrite_vcf(subset_vcf, run_dir / "phased.vcf.gz", [], degrade)
    counter = itertools.count()
    # Every tenth record has no counterpart in the unphased VCF, so its degraded values must stay
    unphased = write_vcf(subset_vcf, run_dir / "unphased.vcf.gz", lambda r: next(counter) % 10 != 0)
    missing = {site_key(r) for r in source[1]} - {site_key(r) for r in load_records(unphased)[1]}
    assert missing

    outputs = run_workflow(
        f"{UTILS}/FillPhasedGenotypes.wdl",
        {
            "phased_vcf": str(phased),
            "phased_vcf_idx": f"{phased}.tbi",
            "unphased_vcf": str(unphased),
            "unphased_vcf_idx": f"{unphased}.tbi",
            "contigs": [CONTIG],
            "prefix": "filled",
            "shard_bin_size": 20000000,
            "utils_docker": utils_docker,
        },
    )

    _, got = load_records(outputs["FillPhasedGenotypes.hiphase_phased_vcf"])
    assert [site_key(r) for r in got] == [site_key(r) for r in source[1]]
    for new, old in zip(got, source[1]):
        if site_key(old) in missing:
            assert new["gts"][blanked] == (None, None)
            assert new["format"][depth_reset].get("DP") == 1
            assert new["format"][depth_reset].get("AL") == ("0" if old["format"][depth_reset].get("AL") is not None else None)
        else:
            assert new["gts"][blanked] == old["gts"][blanked]
            # DP is compared only where the unphased VCF has one to restore
            without_dp = lambda fmt: {k: v for k, v in fmt.items() if k != "DP"}  # noqa: E731
            assert without_dp(new["format"][depth_reset]) == without_dp(old["format"][depth_reset])
        assert new["gts"]["HG00099"] == old["gts"]["HG00099"]


def test_fill_backbone_phased_genotypes_pulls_phase_only_for_unphased_hets(
    run_workflow, run_dir, subset_vcf, utils_docker, source
):
    sample, new_ps = "HG00096", 777

    def unphase(record):
        record.samples[sample].phased = False

    def retag(record):
        if record.samples[sample].phased:
            record.samples[sample]["PS"] = new_ps

    backbone = rewrite_vcf(subset_vcf, run_dir / "backbone.vcf.gz", [], unphase)
    counter = itertools.count()
    # The source for pulling is the same callset with a recognisable PS, lacking every third record
    full = rewrite_vcf(subset_vcf, run_dir / "full.vcf.gz", [], retag)
    notrgt = write_vcf(full, run_dir / "notrgt.vcf.gz", lambda r: next(counter) % 3 != 0)
    present = {site_key(r) for r in load_records(notrgt)[1]}

    inputs = {
        "backbone_phased_vcf": str(backbone),
        "backbone_phased_vcf_idx": f"{backbone}.tbi",
        "backbone_phased_notrgt_vcf": str(notrgt),
        "backbone_phased_notrgt_vcf_idx": f"{notrgt}.tbi",
        "contig": CONTIG,
        "prefix": "pulled",
        "shard_bin_size": 25000000,
        "utils_docker": utils_docker,
    }
    outputs = run_workflow(f"{UTILS}/FillBackbonePhasedGenotypes.wdl", inputs)

    n_het = n_unphased = n_post_pull = 0
    _, got = load_records(outputs["FillBackbonePhasedGenotypes.backbone_merged_vcf"])
    assert [site_key(r) for r in got] == [site_key(r) for r in source[1]]
    with pysam.VariantFile(str(outputs["FillBackbonePhasedGenotypes.backbone_merged_vcf"])) as vcf:
        phased_after = [r.samples[sample].phased for r in vcf]
    for new, old, is_phased in zip(got, source[1], phased_after):
        gt = old["gts"][sample]
        if None in gt or gt[0] == gt[1]:
            assert new["gts"][sample] == gt
            continue
        n_het += 1
        n_unphased += 1
        was_phased = old["phased"][sample]
        if was_phased and site_key(old) in present:
            assert is_phased and new["gts"][sample] == gt
            assert new["format"][sample]["PS"] == new_ps
        else:
            assert not is_phased
            n_post_pull += 1
    stats = {row["SAMPLE"]: row for row in read_dicts(outputs["FillBackbonePhasedGenotypes.backbone_merged_tsv"])}
    assert set(stats) == set(SAMPLES)
    row = stats[sample]
    assert (int(row["N_HET"]), int(row["N_HET_UNPHASED"]), int(row["N_HET_UNPHASED_POST_PULL"])) == (
        n_het,
        n_unphased,
        n_post_pull,
    )
    assert row["PROP_HET_UNPHASED_POST_PULL"] == f"{n_post_pull / n_het:.6f}"
    assert 0 < n_post_pull < n_het


def fill_format_inputs(unfilled, filled, utils_docker, **overrides):
    return {
        "unfilled_vcf": str(unfilled),
        "unfilled_vcf_idx": f"{unfilled}.tbi",
        "filled_vcf": str(filled),
        "filled_vcf_idx": f"{filled}.tbi",
        "contig": CONTIG,
        "prefix": "refilled",
        "transfer_format_fields": ["DP"],
        "drop_format_fields": [],
        "fill_alt_gts": False,
        "fill_ref_gts": False,
        "match_by_id": False,
        "unphase_gts": False,
        "add_missing_pl_via_ad": False,
        "expand_ad_across_alleles": False,
        "split_rnc_across_alleles": False,
        "normalize_unfilled_vcf": False,
        "normalize_filled_vcf": False,
        "utils_docker": utils_docker,
        **overrides,
    }


@pytest.mark.parametrize("match_by_id", [False, True])
def test_fill_format_fields_transfers_fields_and_alt_genotypes_for_the_subset(
    run_workflow, run_dir, subset_vcf, utils_docker, source, match_by_id
):
    degraded_sample = "HG00096"

    def degrade(record):
        record.samples[degraded_sample]["GT"] = (0, 0)
        record.samples[degraded_sample].phased = True
        record.samples[degraded_sample]["DP"] = 0

    unfilled = rewrite_vcf(subset_vcf, run_dir / "unfilled.vcf.gz", [], degrade)
    outputs = run_workflow(
        f"{UTILS}/FillFormatFields.wdl",
        fill_format_inputs(
            unfilled,
            subset_vcf,
            utils_docker,
            match_by_id=match_by_id,
            fill_alt_gts=True,
            unphase_gts=True,
            drop_format_fields=["PL"],
            shard_bin_size_fill=25000000,
            subset_unfilled_vcf_field="allele_type",
            subset_unfilled_vcf_value="snv",
        ),
    )

    samples, got = load_records(outputs["FillFormatFields.refilled_vcf"])
    assert [site_key(r) for r in got] == [site_key(r) for r in source[1]]
    refilled = 0
    for new, old in zip(got, source[1]):
        assert all("PL" not in fmt for fmt in new["format"].values())
        in_subset = old["info"]["allele_type"] == "snv"
        for sample in SAMPLES:
            src_gt = tuple(sorted(old["gts"][sample], key=lambda a: (a is None, a or 0)))
            if sample != degraded_sample:
                assert new["gts"][sample] == src_gt
                continue
            if in_subset and any(a for a in old["gts"][sample] if a):
                assert new["gts"][sample] == src_gt
                refilled += 1
            else:
                assert new["gts"][sample] == (0, 0)
            expected_dp = old["format"][sample].get("DP") if in_subset else 0
            assert new["format"][sample].get("DP") == expected_dp
    assert refilled > 0


def phred_likelihoods(ref_reads, alt_reads):
    if ref_reads + alt_reads == 0:
        return (0, 0, 0)
    log_likelihoods = [(alt_reads * math.log10(p) + ref_reads * math.log10(1 - p)) for p in (0.05, 0.5, 0.95)]
    return tuple(int(round(-10 * (x - max(log_likelihoods)))) for x in log_likelihoods)


def test_fill_format_fields_derives_missing_pl_from_allele_depths(run_workflow, run_dir, subset_vcf, utils_docker, source):
    plain = run_dir / "no_pl.vcf"
    with gzip.open(subset_vcf, "rt") as handle, open(plain, "w") as out:
        for line in handle:
            if line.startswith("##FORMAT=<ID=PL,"):
                continue
            if not line.startswith("#"):
                fields = line.rstrip("\n").split("\t")
                keep = [i for i, key in enumerate(fields[8].split(":")) if key != "PL"]
                for column in range(8, len(fields)):
                    values = fields[column].split(":")
                    fields[column] = ":".join(values[i] for i in keep)
                line = "\t".join(fields) + "\n"
            out.write(line)
    pysam.tabix_index(str(plain), preset="vcf", force=True)
    no_pl = run_dir / "no_pl.vcf.gz"
    assert "PL" not in pysam.VariantFile(str(no_pl)).header.formats

    outputs = run_workflow(
        f"{UTILS}/FillFormatFields.wdl",
        fill_format_inputs(no_pl, subset_vcf, utils_docker, transfer_format_fields=[], add_missing_pl_via_ad=True),
    )

    _, got = load_records(outputs["FillFormatFields.refilled_vcf"])
    derived = 0
    for new, old in zip(got, source[1]):
        if len(old["alts"]) != 1:
            continue
        for sample in SAMPLES:
            ad = old["format"][sample].get("AD")
            if ad is not None and None not in ad and len(old["gts"][sample]) == 2:
                assert new["format"][sample]["PL"] == phred_likelihoods(*ad)
                derived += 1
    assert derived > 100


def test_evaluate_backbone_phasing_counts_switch_and_flip_errors_against_base_phasing(
    run_workflow, run_dir, subset_vcf, utils_docker, source
):
    sample = "HG00096"
    trv_ids = {r["id"] for r in source[1] if r["info"]["allele_type"] == "trv"}

    def comparable(record):
        gt = record["gts"][sample]
        enveloped = "TR_ENVELOPED" in record["info"] and record["info"].get("TRID") in trv_ids
        return (
            record["info"]["allele_type"] != "trv"
            and not enveloped
            and None not in gt
            and gt[0] != gt[1]
            and record["phased"][sample]
        )

    matched = sorted((r["pos"], r["id"]) for r in source[1] if comparable(r))
    # Swapping the haplotypes of one call gives a flip error, and swapping every call from an index onward gives a switch error
    flip_at, switch_from = 10, 40
    swapped = {rid for _, rid in matched[flip_at : flip_at + 1] + matched[switch_from:]}
    assert len(matched) > switch_from + 5

    def swap(record):
        if record.id in swapped:
            record.samples[sample]["GT"] = record.samples[sample]["GT"][::-1]
            record.samples[sample].phased = True

    backbone = rewrite_vcf(subset_vcf, run_dir / "backbone.vcf.gz", [], swap)
    outputs = run_workflow(
        f"{UTILS}/EvaluateBackbonePhasing.wdl",
        {
            "backbone_phased_vcfs": [str(backbone)],
            "backbone_phased_vcf_idxs": [f"{backbone}.tbi"],
            "base_vcfs": [str(subset_vcf)],
            "base_vcf_idxs": [f"{subset_vcf}.tbi"],
            "contigs": [CONTIG],
            "prefix": "phasing",
            "utils_docker": utils_docker,
        },
    )

    rows = {r["sample"]: r for r in read_dicts(outputs["EvaluateBackbonePhasing.outside_tr_table"])}
    assert set(rows) == set(SAMPLES)
    assert (rows[sample]["matched_count"], rows[sample]["switch_error_count"], rows[sample]["flip_error_count"]) == (
        str(len(matched)),
        "1",
        "1",
    )
    for other in SAMPLES[1:]:
        assert (rows[other]["switch_error_count"], rows[other]["flip_error_count"]) == ("0", "0")
    assert read_table(outputs["EvaluateBackbonePhasing.missing_samples"]) == []

    (table,) = outputs["EvaluateBackbonePhasing.vcf_tables"]
    status = {(r["pos"], r["variant_id"]): r[sample] for r in read_dicts(table)}
    assert [status[(str(pos), rid)] for pos, rid in matched[flip_at : flip_at + 2]] == ["FL", "FL"]
    assert status[(str(matched[switch_from][0]), matched[switch_from][1])] == "SW"
    assert {s for k, s in status.items() if k in {(str(p), i) for p, i in matched}} == {"XC", "CN", "FL", "SW"}


def test_concatenate_mosdepth_joins_per_base_beds_and_indexes_the_result(run_workflow, run_dir, utils_docker):
    beds = [
        write_mosdepth_bed(run_dir / "a.regions.bed", [(0, 100, 5), (100, 250, 7)]),
        write_mosdepth_bed(run_dir / "b.regions.bed", [(250, 400, 9)]),
    ]
    outputs = run_workflow(
        f"{UTILS}/ConcatenateMosDepth.wdl",
        {"mosdepth_bed_files": [str(b) for b in beds], "prefix": "cohort", "utils_docker": utils_docker},
    )

    assert read_table(outputs["ConcatenateMosDepth.mosdepth_per_base_combined"]) == [
        [CONTIG, "0", "100", "5"],
        [CONTIG, "100", "250", "7"],
        [CONTIG, "250", "400", "9"],
    ]
    with pysam.TabixFile(
        str(outputs["ConcatenateMosDepth.mosdepth_per_base_combined"]),
        index=outputs["ConcatenateMosDepth.mosdepth_per_base_combined_idx"],
    ) as tbx:
        assert [line.split("\t")[3] for line in tbx.fetch(CONTIG, 240, 260)] == ["7", "9"]


def test_create_cohort_metadata_joins_sex_from_ped_with_ancestry(run_workflow, run_dir, utils_docker):
    ped = write_table(
        run_dir / "cohort.ped",
        [["f1", "S1", 0, 0, 1, 0], ["f1", "S2", 0, 0, 2, 0], ["f2", "S3", 0, 0, 0, 0], ["f2", "S4", 0, 0, 2, 0]],
    )
    ancestry = write_table(run_dir / "ancestry.tsv", [["S1", "afr"], ["S2", "nfe"], ["S3", "eas"]])
    outputs = run_workflow(
        f"{UTILS}/CreateCohortMetadata.wdl",
        {"ped_file": str(ped), "ancestry_file": str(ancestry), "prefix": "cohort", "utils_docker": utils_docker},
    )

    assert read_table(outputs["CreateCohortMetadata.metadata"]) == [
        ["SampleId", "Sex", "Population"],
        ["S1", "male", "afr"],
        ["S2", "female", "nfe"],
        ["S3", "unknown", "eas"],
        ["S4", "female", "unknown"],
    ]


def test_create_cohort_coverage_summary_reports_per_bin_statistics_across_samples(run_workflow, run_dir, utils_docker):
    depths = [
        [(0, 3000, 10), (3000, 7000, 30), (7000, 20000, 20)],
        [(0, 5000, 12), (5000, 20000, 40)],
        [(0, 20000, 25)],
    ]
    beds = [write_mosdepth_bed(run_dir / f"s{i}.regions.bed", runs) for i, runs in enumerate(depths)]
    fai = write_table(run_dir / "small.fai", [[CONTIG, 20000, 7, 60, 61]])
    bin_size, thresholds = 1000, [15, 30]
    outputs = run_workflow(
        f"{UTILS}/CreateCohortCoverageSummary.wdl",
        {
            "mosdepth_bed_files": [str(b) for b in beds],
            "mosdepth_bed_idx": [f"{b}.tbi" for b in beds],
            "ref_fai": str(fai),
            "contigs": [CONTIG],
            "prefix": "summary",
            "window_size": 10000,
            "bin_size": bin_size,
            "thresholds": thresholds,
            "utils_docker": utils_docker,
        },
    )

    def depth_at(runs, position):
        return next(depth for start, end, depth in runs if start <= position < end)

    expected = []
    for start in range(0, 20000, bin_size):
        values = sorted(depth_at(runs, start) for runs in depths)
        row = [CONTIG, start + 1, start + bin_size, f"{sum(values) / 3:.6f}", values[1], sum(values)]
        row += [f"{sum(v > t for v in values) / 3:.6f}" for t in thresholds]
        expected.append([str(v) for v in row])
    assert read_table(outputs["CreateCohortCoverageSummary.binned_coverage_tsv"]) == expected


def test_create_sample_read_counts_writes_gatk_counts_of_median_depth_per_full_bin(run_workflow, run_dir, utils_docker):
    # The bin from 2000 to 3000 is 700 bases at depth 10 then 300 at depth 30, and the final partial bin is dropped
    runs = [(0, 2700, 10), (2700, 4000, 30), (4000, 9500, 20)]
    bed = write_mosdepth_bed(run_dir / "sample.regions.bed", runs)
    dictionary = run_dir / "ref.dict"
    dictionary.write_text(f"@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:{CONTIG}\tLN:9500\n")

    outputs = run_workflow(
        f"{UTILS}/CreateSampleReadCounts.wdl",
        {
            "mosdepth_bed_files": [str(bed)],
            "mosdepth_bed_idx": [f"{bed}.tbi"],
            "ref_dict": str(dictionary),
            "contigs": [CONTIG],
            "prefix": "reads",
            "bin_size": 1000,
            "sample_id": "S1",
            "utils_docker": utils_docker,
        },
    )

    rows = read_table(outputs["CreateSampleReadCounts.binned_read_counts"])
    assert rows[:3] == [["@HD", "VN:1.6", "SO:coordinate"], ["@SQ", "SN:chr22", "LN:9500"], ["@RG", "ID:GATKCopyNumber", "SM:S1"]]
    assert rows[3] == ["CONTIG", "START", "END", "COUNT"]
    medians = [10, 10, 10, 30, 20, 20, 20, 20, 20]
    assert [(r[1], r[2], r[3]) for r in rows[4:]] == [
        (str(i * 1000 + 1), str((i + 1) * 1000), str(m)) for i, m in enumerate(medians)
    ]


def test_filter_low_coverage_genotypes_masks_low_depth_calls_per_sample_and_reports_them(
    run_workflow, run_dir, vcf_inputs, source
):
    depths = {s: sorted(r["format"][s]["DP"] for r in source[1] if r["format"][s].get("DP") is not None) for s in SAMPLES}
    cutoffs = {s: float(depths[s][len(depths[s]) // 2]) for s in SAMPLES}
    cutoffs_tsv = write_table(run_dir / "cutoffs.tsv", [[s, c] for s, c in cutoffs.items()], header=["sample_id", "cutoff"])
    ped = write_table(run_dir / "cohort.ped", [[f"fam{i}", s, 0, 0, 1 + i % 2, 0] for i, s in enumerate(SAMPLES)])

    outputs = run_workflow(
        f"{UTILS}/FilterLowCoverageGenotypes.wdl",
        {
            **vcf_inputs,
            "prefix": "lowcov",
            "sample_cutoffs_tsv": str(cutoffs_tsv),
            "ped": str(ped),
            "subset_unfilled_vcf_field": "allele_type",
            "subset_unfilled_vcf_value": "snv",
            "records_per_shard": 150,
        },
    )

    _, got = load_records(outputs["FilterLowCoverageGenotypes.filtered_vcf"])
    assert [site_key(r) for r in got] == [site_key(r) for r in source[1]]

    def allele_counts(record, gts):
        counts = [sum(a == i for gt in gts for a in gt if a is not None) for i in range(1, len(record["alts"]) + 1)]
        return ",".join(map(str, counts))

    expected_rows = []
    for new, old in zip(got, source[1]):
        masked = []
        for s in SAMPLES:
            dp, gt = old["format"][s].get("DP"), old["gts"][s]
            low = old["info"]["allele_type"] == "snv" and dp is not None and dp <= cutoffs[s] and any(a is not None for a in gt)
            if low:
                masked.append(s)
                assert new["gts"][s] == (None,) * len(gt) and new["format"][s]["RNC"] == ("D",) * len(gt)
            else:
                assert new["gts"][s] == gt
        if masked:
            expected_rows.append(
                [CONTIG, str(old["pos"]), old["ref"], ",".join(old["alts"]), old["id"], allele_counts(old, old["gts"].values())]
                + [allele_counts(old, [new["gts"][s] for s in SAMPLES]), str(len(masked)), ",".join(masked)]
            )
    rows = read_table(outputs["FilterLowCoverageGenotypes.filtered_genotypes_tsv"])
    assert rows[0][:5] == ["CHROM", "POS", "REF", "ALT", "ID"] and rows[0][-1] == "filtered_sample_ids"
    assert rows[1:] == expected_rows
    assert len(expected_rows) > 20


def test_filter_low_coverage_regions_flags_records_mostly_inside_low_coverage_bed(run_workflow, run_dir, vcf_inputs, source):
    cutoff = 0.8
    intervals = []
    for i, r in enumerate(source[1]):
        start, end = r["pos"] - 1, r["pos"] - 1 + len(r["ref"])
        if i % 3 == 0:
            intervals.append((start - 5, end + 5))
        elif i % 3 == 1 and end - start > 2:
            # Covers about half of the record, below the cutoff
            intervals.append((start, start + (end - start) // 2))
    bed = write_table(run_dir / "low.bed", [[CONTIG, s, e] for s, e in sorted(intervals)])

    outputs = run_workflow(
        f"{UTILS}/FilterLowCoverageRegions.wdl",
        {
            **vcf_inputs,
            "contigs": [CONTIG],
            "prefix": "regions",
            "low_coverage_regions_bed": str(bed),
            "min_region_coverage_cutoff": cutoff,
            "records_per_shard": 150,
        },
    )

    covered = set()
    for s, e in intervals:
        covered.update(range(s, e))
    _, got = load_records(outputs["FilterLowCoverageRegions.low_coverage_region_filtered_vcf"])
    assert [site_key(r) for r in got] == [site_key(r) for r in source[1]]
    flagged = 0
    for new, old in zip(got, source[1]):
        start, end = old["pos"] - 1, old["pos"] - 1 + len(old["ref"])
        if sum(b in covered for b in range(start, end)) / (end - start) >= cutoff:
            flagged += 1
            assert new["filter"] == tuple([f for f in old["filter"] if f != "PASS"] + ["LOW_COVERAGE_REGION"])
        else:
            assert new["filter"] == old["filter"]
        assert new["gts"] == old["gts"]
    assert 50 < flagged < len(source[1])


def test_identify_low_coverage_regions_flags_bins_below_per_sample_cutoffs_shared_by_enough_samples(
    run_workflow, run_dir, utils_docker
):
    bin_size, length = 10000, 50000
    runs = {
        "S1": [(0, 10000, 30), (10000, 20000, 3), (20000, 50000, 30)],
        "S2": [(0, 20000, 20), (20000, 30000, 2), (30000, 50000, 20)],
        "S3": [(0, 10000, 40), (10000, 20000, 4), (20000, 40000, 40), (40000, 43000, 40), (43000, 50000, 50)],
    }
    beds = {s: write_mosdepth_bed(run_dir / f"{s}.regions.bed", r) for s, r in runs.items()}
    # S3 has an unknown sex, so it is inferred from the (absent) chrY coverage
    ped = write_table(run_dir / "cohort.ped", [["f1", "S1", 0, 0, 1, 0], ["f2", "S2", 0, 0, 2, 0], ["f3", "S3", 0, 0, 0, 0]])

    outputs = run_workflow(
        f"{UTILS}/IdentifyLowCoverageRegions.wdl",
        {
            "mosdepth_bed_files": [str(beds[s]) for s in runs],
            "sample_ids": list(runs),
            "ped": str(ped),
            "prefix": "lowcov",
            "bin_size": bin_size,
            "utils_docker": utils_docker,
        },
    )

    # Independent model: weighted bin medians, then the median bin and a 20% cutoff per sample
    low_bins = {}
    expected_cutoffs = [["sample_id", "cutoff", "median_coverage"]]
    for sample, sample_runs in runs.items():
        per_base = [depth for start, end, depth in sample_runs for _ in range(end - start)]
        medians = [statistics.median(per_base[i : i + bin_size]) for i in range(0, length, bin_size)]
        median = sorted(medians)[math.ceil(0.5 * len(medians)) - 1]
        cutoff = math.floor(median * 0.2)
        low_bins[sample] = {i for i, m in enumerate(medians) if m <= cutoff}
        expected_cutoffs.append([sample, str(cutoff), str(int(median))])
    assert read_table(outputs["IdentifyLowCoverageRegions.sample_cutoffs_tsv"]) == expected_cutoffs
    assert [row[1:] for row in expected_cutoffs[1:]] == [["6", "30"], ["4", "20"], ["8", "40"]]

    counts = {i: sum(i in bins for bins in low_bins.values()) for i in range(length // bin_size)}
    expected_counts = [[CONTIG, str(i * bin_size), str((i + 1) * bin_size), str(n), "3"] for i, n in counts.items()]
    rows = read_table(outputs["IdentifyLowCoverageRegions.cohort_coverage_counts_tsv"])
    assert rows[0] == ["chrom", "start", "end", "low_coverage_sample_count", "eligible_sample_count"]
    assert rows[1:] == expected_counts
    failed = [[CONTIG, str(i * bin_size), str((i + 1) * bin_size)] for i, n in counts.items() if n / 3 >= 0.5]
    assert failed == [[CONTIG, "10000", "20000"]]
    assert read_table(outputs["IdentifyLowCoverageRegions.failed_bins_bed"]) == failed

    for key, expected in (
        ("sample_histograms_tar", {f"lowcov.{s}.coverage_histogram.png" for s in runs}),
        ("chromosome_low_coverage_tar", {f"{CONTIG}.low_coverage.png"}),
    ):
        with tarfile.open(outputs[f"IdentifyLowCoverageRegions.{key}"]) as tar:
            assert {Path(name).name for name in tar.getnames() if name.endswith(".png")} == expected


def test_create_cohort_methylation_file_merges_per_sample_scores_into_site_matrices(run_workflow, run_dir, utils_docker):
    samples = ["M1", "M2", "M3"]
    rng = random.Random(7)
    site_pool = [(start, start + 2) for start in range(1000, 1400, 10)]
    scores = {}  # (kind, sample) -> {site: score}
    inputs = {"combined_beds": [], "hap1_beds": [], "hap2_beds": []}
    for kind, key in (("combined", "combined_beds"), ("hap1", "hap1_beds"), ("hap2", "hap2_beds")):
        for sample in samples:
            sites = sorted(rng.sample(site_pool, 25))
            scores[kind, sample] = {site: rng.randint(0, 100) for site in sites}
            rows = [[CONTIG, s, e, score] for (s, e), score in scores[kind, sample].items()]
            inputs[key].append(str(index_table(write_table(run_dir / f"{sample}.{kind}.bed", rows), "bed")))
    for key in list(inputs):
        inputs[key.replace("beds", "bed_idxs")] = [f"{path}.tbi" for path in inputs[key]]

    outputs = run_workflow(
        f"{UTILS}/CreateCohortMethylationFile.wdl",
        {
            **inputs,
            "sample_ids": samples,
            "contigs": [CONTIG],
            "prefix": "meth",
            "samples_per_shard": 2,
            "utils_docker": utils_docker,
        },
    )

    shards = [samples[:2], samples[2:]]
    combined = [(s, "combined", s) for s in samples]
    # Haplotype columns are ordered per sample shard: all hap1 columns of the shard, then all hap2 columns
    haplotype = [(f"{s}_hap{h}", f"hap{h}", s) for shard in shards for h in (1, 2) for s in shard]
    for output, columns in (("combined_methylation_beds", combined), ("haplotype_methylation_beds", haplotype)):
        sites = sorted({site for _, kind, sample in columns for site in scores[kind, sample]})
        expected = [["#chrom", "start", "end"] + [name for name, _, _ in columns]]
        for start, end in sites:
            expected.append(
                [CONTIG, str(start), str(end)] + [str(scores[kind, s].get((start, end), ".")) for _, kind, s in columns]
            )
        assert read_table(outputs[f"CreateCohortMethylationFile.{output}"][0]) == expected


def ncr_size_bin(record):
    length = record["info"].get("allele_length")
    if length is None:
        return "missing"
    return "0-0" if abs(length) < 1 else "1-49" if abs(length) < 50 else "50-499" if abs(length) < 500 else "500+"


def trio_outcome(record, child, parents):
    """Transmission of a proband's alternate call: inherited, uninherited, skipped (a parent is a no-call) or None."""
    if not is_carrier(record, child):
        return None
    if any(is_carrier(record, p) for p in parents):
        return "inherited"
    return "skipped" if any(None in record["gts"][p] for p in parents) else "uninherited"


def test_evaluate_ncr_cutoffs_tabulates_trio_transmission_cumulatively_by_ncr_cutoff(run_workflow, run_dir, vcf_inputs, source):
    trios = [("HG00096", "HG00097", "HG00099"), ("HG00126", "HG00128", "HG00097")]
    ped = write_table(
        run_dir / "trios.ped", [[f"fam{i}", child, father, mother, 1, 0] for i, (child, father, mother) in enumerate(trios)]
    )

    outputs = run_workflow(
        f"{UTILS}/EvaluateNCRCutoffs.wdl",
        {**vcf_inputs, "contig": CONTIG, "prefix": "ncr", "ped": str(ped), "records_per_shard": 150},
    )

    cutoffs = [round(i * 0.01, 10) for i in range(101)]
    rows = read_table(outputs["EvaluateNCRCutoffs.ncr_cutoffs_tsv"])
    assert rows[0] == [
        "allele_type",
        "size_bin",
        "filter",
        "region",
        "cutoff",
        "prec_numerator",
        "prec_denominator",
        "prec_skipped",
        "rec_numerator",
        "rec_denominator",
        "sites_passing",
        "sites_total",
    ]
    by_stratum = collections.defaultdict(list)
    for row in rows[1:]:
        by_stratum[tuple(row[:4])].append(row[4:])

    strata = {
        ("all", "all", "all", "all"): lambda r: True,
        ("snv", "all", "all", "all"): lambda r: r["info"]["allele_type"] == "snv",
        ("del", "50-499", "all", "all"): lambda r: r["info"]["allele_type"] == "del" and ncr_size_bin(r) == "50-499",
        ("all", "all", "PASS", "SR"): lambda r: r["filter"] == ("PASS",) and r["info"]["REGION"] == "SR",
        ("trv", "missing", "all", "all"): lambda r: r["info"]["allele_type"] == "trv",
    }
    for stratum, member in strata.items():
        records = [r for r in source[1] if member(r)]
        assert records
        expected = []
        for cutoff in cutoffs:
            passing = [r for r in records if float(f"{r['info']['NCR']:.6g}") <= cutoff + 1e-9]
            outcomes = [trio_outcome(r, child, (father, mother)) for r in passing for child, father, mother in trios]
            inherited, uninherited = outcomes.count("inherited"), outcomes.count("uninherited")
            calls = sum(is_carrier(r, s) for r in passing for s in SAMPLES)
            expected.append(
                [f"{cutoff:g}", str(inherited), str(inherited + uninherited), str(outcomes.count("skipped")), str(calls)]
            )
            expected[-1] += [
                str(sum(is_carrier(r, s) for r in records for s in SAMPLES)),
                str(len(passing)),
                str(len(records)),
            ]
        assert by_stratum[stratum] == expected, stratum
    assert any(int(row[1]) > 0 for row in by_stratum[("all", "all", "all", "all")])
    assert any(int(row[3]) > 0 for row in by_stratum[("all", "all", "all", "all")])


@pytest.mark.parametrize("pass_only", [False, True])
def test_filter_ncr_cutoffs_flags_sites_above_their_stratum_cutoff_and_drops_sites_without_alt_calls(
    run_workflow, run_dir, subset_vcf, utils_docker, pass_only
):
    cutoff_rows = [
        ["snv", 0, 0, "ALL", 0.05],
        ["del", 1, 100000, "ALL", 0.1],
        ["ins", 1, 100000, "SR", 0.02],
        ["ins", 1, 100000, "US", 0.5],
    ]
    cutoffs = write_table(
        run_dir / "cutoffs.tsv", cutoff_rows, header=["allele_type", "min_length", "max_length", "region", "cutoff"]
    )

    # Every seventh site has no genotype calls at all
    count = itertools.count()

    def drop_calls(record):
        if next(count) % 7 == 0:
            for sample in record.samples.values():
                sample["GT"] = (None, None)

    vcf = rewrite_vcf(subset_vcf, run_dir / "input.vcf.gz", [], drop_calls)
    _, records = load_records(vcf)

    outputs = run_workflow(
        f"{UTILS}/FilterNCRCutoffs.wdl",
        {
            "vcf": str(vcf),
            "vcf_idx": f"{vcf}.tbi",
            "utils_docker": utils_docker,
            "prefix": "ncrfilter",
            "remove_no_calls": True,
            "ncr_cutoffs_tsv": str(cutoffs),
            "ncr_filter_pass_only": pass_only,
            "records_per_shard": 150,
        },
    )

    kept = [r for r in records if any(is_carrier(r, s) for s in SAMPLES)]
    assert 0 < len(kept) < len(records)
    got_header, got = load_records(outputs["FilterNCRCutoffs.filtered_vcf"])
    assert [site_key(r) for r in got] == [site_key(r) for r in kept]
    assert "HIGH_NCR" in pysam.VariantFile(outputs["FilterNCRCutoffs.filtered_vcf"]).header.filters

    flagged = 0
    for new, old in zip(got, kept):
        length = old["info"].get("allele_length")
        cutoff = next(
            (
                c
                for kind, lo, hi, region, c in cutoff_rows
                if length is not None
                and kind == old["info"]["allele_type"]
                and lo <= abs(length) <= hi
                and region in ("ALL", old["info"]["REGION"])
            ),
            None,
        )
        applies = cutoff is not None and (not pass_only or "PASS" in old["filter"])
        if applies and float(f"{old['info']['NCR']:.6g}") > cutoff + 1e-9:
            flagged += 1
            assert new["filter"] == tuple([f for f in old["filter"] if f != "PASS"] + ["HIGH_NCR"])
        else:
            assert new["filter"] == old["filter"]
        assert new["gts"] == old["gts"]
    assert flagged > 10


def test_integrate_trs_envelops_base_calls_in_tr_loci_and_tags_the_tr_records(run_workflow, run_dir, subset_vcf, utils_docker):
    def is_trv(record):
        return record.info["allele_type"] == "trv"

    # Every fifth TR record has no FILTER value
    count = itertools.count()

    def clear_some_filters(record):
        if next(count) % 5 == 0:
            record.filter.clear()

    tr_vcf = rewrite_vcf(
        write_vcf(subset_vcf, run_dir / "trs.vcf.gz", is_trv), run_dir / "trs.nofilter.vcf.gz", [], clear_some_filters
    )

    # Base calls: the non-TR records, no TR annotations, plus a copy of a base call placed inside every long TR locus
    # and a deletion hanging over the end of every second one
    with pysam.VariantFile(str(subset_vcf)) as vcf_in:
        header = vcf_in.header.copy()
        for key in ("TR_ENVELOPED", "TRID"):
            header.info.remove_header(key)
        records = list(vcf_in)
    snv = next(r for r in records if r.info["allele_type"] == "snv")
    deletion = next(r for r in records if r.info["allele_type"] == "del" and len(r.ref) > 5)
    base = [r for r in records if not is_trv(r)]
    long_trs = [r for r in records if is_trv(r) and len(r.ref) >= 10]
    for i, tr in enumerate(long_trs):
        template, pos = (snv, tr.pos + 3) if i % 2 == 0 else (deletion, tr.pos + len(tr.ref) - 2)
        clone = template.copy()
        clone.pos = pos
        base.append(clone)
    base.sort(key=lambda r: r.pos)
    base_plain = run_dir / "base.vcf"
    with pysam.VariantFile(str(base_plain), "w", header=header) as out:
        for record in base:
            for key in ("TR_ENVELOPED", "TRID"):
                if key in record.info:
                    del record.info[key]
            record.translate(header)
            out.write(record)
    pysam.tabix_index(str(base_plain), preset="vcf", force=True)
    base_vcf = Path(f"{base_plain}.gz")

    # Two catalogs, the later taking precedence for loci listed in both
    _, tr_loaded = load_records(tr_vcf)
    trids = [r["info"]["TRID"] for r in tr_loaded]
    catalogs = []
    for name, ids in (("catA", trids[::2]), ("catB", trids[::3])):
        path = run_dir / f"{name}.bed"
        write_table(path, [[CONTIG, 1, 2, f"ID={trid};MOTIFS=A"] for trid in ids])
        with gzip.open(f"{path}.gz", "wt") as handle:
            handle.write(path.read_text())
        catalogs.append(f"{path}.gz")
    sources = {trid: "catA" for trid in trids[::2]} | {trid: "catB" for trid in trids[::3]}

    outputs = run_workflow(
        f"{UTILS}/IntegrateTRs.wdl",
        {
            "vcf": str(base_vcf),
            "vcf_idx": f"{base_vcf}.tbi",
            "tr_vcf": str(tr_vcf),
            "tr_vcf_idx": f"{tr_vcf}.tbi",
            "contigs": [CONTIG],
            "prefix": "integrated",
            "sample_ids": SAMPLES,
            "tr_catalogs": catalogs,
            "tr_catalog_ids": ["catA", "catB"],
            "utils_docker": utils_docker,
        },
    )

    # TR records are renamed after their position and reference length
    tr_ids = [f"{r['contig']}-{r['pos']}-TRV-{len(r['ref']) - 1}" for r in tr_loaded]
    assert len(set(tr_ids)) == len(tr_ids)

    _, base_loaded = load_records(base_vcf)
    got_samples, got = load_records(outputs["IntegrateTRs.tr_annotated_vcf"])
    assert got_samples == SAMPLES
    got_by_key = {(*site_key(r), r["id"]): r for r in got}
    assert len(got_by_key) == len(got) == len(base_loaded) + len(tr_loaded)

    enveloped = 0
    for record in base_loaded:
        start, end = record["pos"], record["pos"] + len(record["ref"])
        tr_id = next((i for t, i in zip(tr_loaded, tr_ids) if t["pos"] <= start and end <= t["pos"] + len(t["ref"])), None)
        new = got_by_key[(*site_key(record), record["id"])]
        assert new["gts"] == record["gts"]
        if tr_id is None:
            assert "TR_ENVELOPED" not in new["info"] and "TRID" not in new["info"]
        else:
            enveloped += 1
            assert new["info"]["TR_ENVELOPED"] and new["info"]["TRID"] == tr_id
    assert enveloped >= len(long_trs) // 2

    for record, tr_id in zip(tr_loaded, tr_ids):
        new = got_by_key[(*site_key(record), tr_id)]
        assert new["filter"] == ("PASS",)
        assert new["info"]["allele_type"] == "trv"
        assert new["info"].get("SOURCE") == sources.get(record["info"]["TRID"], record["info"].get("SOURCE"))
        assert new["gts"] == record["gts"]


def test_summarize_trv_alleles_per_genome_counts_alt_alleles_overlapping_coding_sequence(
    run_workflow, run_dir, vcf_inputs, source, utils_docker
):
    trvs = [r for r in source[1] if r["info"]["allele_type"] == "trv"]
    # A CDS inside the span of every third TR locus, named after the locus; two loci share a gene
    cds = []
    for i, record in enumerate(trvs[::3]):
        gene = f"GENE{i // 2}"
        cds.append((record["pos"] + 2, record["pos"] + len(record["ref"]) // 2, gene))
    gtf = run_dir / "coding.gtf"
    gtf.write_text(
        "".join(f'{CONTIG}\tHAVANA\tCDS\t{s}\t{e}\t.\t+\t0\tgene_id "{g}"; gene_name "{g}";\n' for s, e, g in cds)
        + f'{CONTIG}\tHAVANA\texon\t1\t50818468\t.\t+\t.\tgene_id "X"; gene_name "NOTCODING";\n'
    )

    outputs = run_workflow(
        f"{UTILS}/SummarizeTRVAllelesPerGenome.wdl",
        {
            "vcfs": [vcf_inputs["vcf"]],
            "vcf_idxs": [vcf_inputs["vcf_idx"]],
            "prefix": "trvs",
            "coding_gtf": str(gtf),
            "records_per_shard": 30,
            "utils_docker": utils_docker,
        },
    )

    def trim(ref, alt):
        head = 0
        while head < min(len(ref), len(alt)) and ref[head] == alt[head]:
            head += 1
        ref, alt = ref[head:], alt[head:]
        tail = 0
        while tail < min(len(ref), len(alt)) and ref[-1 - tail] == alt[-1 - tail]:
            tail += 1
        return head, len(ref) - tail

    def genes_over(start, end, contained=False):
        return {g for s, e, g in cds if (s <= start and e >= end if contained else s <= end and e >= start)}

    def genes_of_allele(record, alt):
        head, ref_len = trim(record["ref"], alt)
        start = record["pos"] + head
        return genes_over(start, start + ref_len - 1) if ref_len else genes_over(start - 1, start, contained=True)

    def carried(record, sample):
        return {a for a in record["gts"][sample] if a not in (0, None)}

    n_alts = alt_sum = coding_alt_sum = locus_sum = 0
    coding_rows, sample_genes, coding_alleles = [], {s: set() for s in SAMPLES}, 0
    for record in trvs:
        n_alts += len(record["alts"])
        per_sample = {s: carried(record, s) for s in SAMPLES}
        alt_sum += sum(len(c) for c in per_sample.values())
        if genes_over(record["pos"], record["pos"] + len(record["ref"]) - 1):
            locus_sum += sum(len(c) for c in per_sample.values())
        for k, alt in enumerate(record["alts"], start=1):
            genes = genes_of_allele(record, alt)
            if not genes:
                continue
            coding_alleles += 1
            n_carriers = sum(k in c for c in per_sample.values())
            coding_alt_sum += n_carriers
            coding_rows.append([CONTIG, str(record["pos"]), record["id"], str(k), str(n_carriers), ",".join(sorted(genes))])
            for sample, c in per_sample.items():
                if k in c:
                    sample_genes[sample] |= genes
    gene_alleles = collections.Counter(g for row in coding_rows for g in row[5].split(","))
    assert coding_rows and len(gene_alleles) > 3

    summary = dict(row[:2] for row in read_table(outputs["SummarizeTRVAllelesPerGenome.summary_tsv"])[1:])
    n_samples = len(SAMPLES)
    carrier_count = sum(len(g) for g in sample_genes.values())
    assert summary == {
        "n_samples": str(n_samples),
        "n_trv_records": str(len(trvs)),
        "n_trv_alt_alleles": str(n_alts),
        "n_coding_trv_alt_alleles": str(coding_alleles),
        "n_coding_genes": str(len(gene_alleles)),
        "alt_alleles_per_genome": f"{alt_sum / n_samples:.4f}",
        "coding_alt_alleles_per_genome": f"{coding_alt_sum / n_samples:.4f}",
        "coding_locus_alt_alleles_per_genome": f"{locus_sum / n_samples:.4f}",
        "coding_genes_per_genome": f"{carrier_count / n_samples:.4f}",
    }
    assert outputs["SummarizeTRVAllelesPerGenome.summary_sentence"].strip() == (
        f"An average of {alt_sum / n_samples:,.1f} alternative alleles were observed per individual genome, "
        f"including {coding_alt_sum / n_samples:,.1f} that reside within coding sequences "
        f"of {carrier_count / n_samples:,.1f} genes"
    )

    rows = read_table(outputs["SummarizeTRVAllelesPerGenome.coding_alleles_tsv"])
    assert rows[0] == ["contig", "pos", "id", "alt_index", "n_carriers", "genes"]
    assert sorted(rows[1:], key=lambda r: (int(r[1]), int(r[3]))) == coding_rows
    carriers = collections.Counter(g for genes in sample_genes.values() for g in genes)
    assert read_table(outputs["SummarizeTRVAllelesPerGenome.coding_genes_tsv"]) == [
        ["gene_name", "n_coding_alt_alleles", "n_carrier_samples"],
        *[[g, str(gene_alleles[g]), str(carriers[g])] for g in sorted(gene_alleles)],
    ]


SUMMARY_SIZE_LABELS = ["0-0", "1-49", "50-499", "500+"]
SUMMARY_COLUMNS = ["SNV"] + [f"{t} {s}" for t in ("DEL", "INS") for s in SUMMARY_SIZE_LABELS] + ["TRV", "Other"]


def summary_column(record):
    kind, length = record["info"]["allele_type"], record["info"].get("allele_length")
    if kind in ("snv", "trv"):
        return kind.upper()
    if length is None or kind not in ("ins", "dup", "del"):
        return "Other"
    label = SUMMARY_SIZE_LABELS[sum(abs(length) >= edge for edge in (1, 50, 500))]
    return f"{'DEL' if kind == 'del' else 'INS'} {label}"


def summary_rows(record):
    """The (category, annotation) rows of the annotation summary that count a record."""
    info = record["info"]
    kind, subtype = info["allele_type"], info.get("allele_subtype") or ""
    rows = [("All", "All"), ("VEP", "All"), ("VEP", "Other"), ("Concordance", "All")]
    if subtype in ("alu_ins", "alu_del"):
        rows += [("ME", "All"), ("ME", "ALU")]
    if subtype == "tandem_dup" or subtype == "dup_interspersed":
        rows.append(("Duplication", "All"))
        rows.append(
            ("Duplication", {"dup": "Filtered Tandem", "ins": "SVAN Tandem"}[kind] if subtype == "tandem_dup" else "Interspersed")
        )
    if "TR_PARSED" in info:
        rows.append(("TR Parsed", "All"))
    if "PREDICTED_LOF" in info:
        rows += [("SVAnnotate", "All"), ("SVAnnotate", "LoF")]
    no_dbsnp, no_gnomad = "dbSNP_ID" not in info and "dbGaP_ID" not in info, "gnomAD_V4_match_ID" not in info
    rows += [("Concordance", name) for name, missing in (("dbSNP Missing", no_dbsnp), ("gnomAD Missing", no_gnomad)) if missing]
    if no_dbsnp and no_gnomad:
        rows.append(("Concordance", "dbSNP/gnomAD Missing"))
    return rows


def test_summarize_annotations_tabulates_site_sample_allele_list_and_gene_summaries_across_shards(
    run_workflow, vcf_inputs, source
):
    records = source[1]
    outputs = run_workflow(
        f"{UTILS}/SummarizeAnnotations.wdl",
        {
            "vcfs": [vcf_inputs["vcf"]],
            "vcf_idxs": [vcf_inputs["vcf_idx"]],
            "prefix": "summary",
            "utils_docker": vcf_inputs["utils_docker"],
            "create_per_sample": True,
            "create_per_allele": True,
            "create_list": True,
            "create_functional": True,
            "split_by_region": True,
            "records_per_shard": 150,
        },
    )
    prefix = "SummarizeAnnotations"
    alt_alleles = {id(r): [sum(a not in (0, None) for a in r["gts"][s]) for s in SAMPLES] for r in records}

    def expected_table(weight):
        table = collections.defaultdict(lambda: collections.defaultdict(float))
        for r in records:
            tr = "TR" if "TR_ENVELOPED" in r["info"] else "Not TR"
            for row in summary_rows(r):
                for status in ("All", tr):
                    for region in ("All", r["info"]["REGION"]):
                        table[(*row, status, region)][summary_column(r)] += weight(r)
        return table

    kept_types = {summary_column(r).split()[0] for r in records}
    columns = [c for c in SUMMARY_COLUMNS if c.split()[0] in kept_types]
    assert {"SNV", "DEL 50-499", "INS 1-49", "TRV"} <= set(columns)
    for name, weight, fmt in (
        ("summary_sites_tsv", lambda r: 1, lambda v: str(int(v))),
        ("summary_samples_tsv", lambda r: sum(a > 0 for a in alt_alleles[id(r)]) / len(SAMPLES), lambda v: f"{v:.2f}"),
        ("summary_alleles_tsv", lambda r: sum(alt_alleles[id(r)]) / len(SAMPLES), lambda v: f"{v:.2f}"),
    ):
        rows = read_table(outputs[f"{prefix}.{name}"])
        assert rows[0] == ["category", "sub_category", "tr_status", "region"] + columns
        got = {tuple(row[:4]): row[4:] for row in rows[1:]}
        expected = expected_table(weight)
        assert len(expected) > 100
        for key, values in expected.items():
            assert got[key] == [fmt(values.get(c, 0.0)) for c in columns], (name, key)
        assert got[("All", "All", "All", "All")] != ["0"] * len(columns)

    listing = {row[0]: row for row in read_table(outputs[f"{prefix}.summary_list_tsv"])[1:]}
    header = read_table(outputs[f"{prefix}.summary_list_tsv"])[0]
    assert header[:2] == ["variant_id", "classification"]
    column_of = {name: i for i, name in enumerate(header)}
    assert len(listing) == len({r["id"] for r in records})
    for r in records:
        row = listing[r["id"]]
        assert row[1] == summary_column(r)
        members = {annotation for _, annotation in summary_rows(r)}
        for label in ("ALU", "Filtered Tandem", "SVAN Tandem", "Interspersed", "dbSNP Missing", "gnomAD Missing"):
            assert row[column_of[label]] == ("1" if label in members else "0"), (r["id"], label)

    for name, weight in (
        ("summary_functional_tsv", lambda r: 1),
        ("summary_functional_samples_tsv", lambda r: sum(a > 0 for a in alt_alleles[id(r)])),
        ("summary_functional_alleles_tsv", lambda r: sum(alt_alleles[id(r)])),
    ):
        rows = read_table(outputs[f"{prefix}.{name}"])
        got = {row[0]: dict(zip(rows[0][1:], row[1:])) for row in rows[1:]}
        expected = collections.defaultdict(collections.Counter)
        for r in records:
            af, ac = r["info"]["AF"][0], abs(int(r["info"]["AC"][0]))
            buckets = ["all"] + ["singleton"] * (ac == 1) + ["rare"] * (af < 0.01) + ["common"] * (af > 0.05)
            for field in ("PREDICTED_INTRONIC", "PREDICTED_LOF"):
                value = r["info"].get(field)
                for gene in {g.strip() for item in ([value] if isinstance(value, str) else value or []) for g in item.split(",")}:
                    for bucket in buckets:
                        expected[gene][f"{field}.{bucket}"] += weight(r)
        assert len(expected) > 10
        for gene, counts in expected.items():
            for column, count in counts.items():
                assert got[gene][column] == str(count), (name, gene, column)


def trim_alleles(ref, alt):
    while len(ref) > 1 and len(alt) > 1 and ref[-1] == alt[-1]:
        ref, alt = ref[:-1], alt[:-1]
    while len(ref) > 1 and len(alt) > 1 and ref[0] == alt[0]:
        ref, alt = ref[1:], alt[1:]
    return ref, alt


def is_snv_ins_or_del(ref, alt):
    trimmed_ref, trimmed_alt = trim_alleles(ref, alt)
    return len(trimmed_ref) != len(trimmed_alt) or (len(trimmed_ref) == 1 and trimmed_ref != trimmed_alt)


def test_preprocess_vcfs_splits_multiallelics_types_alleles_filters_by_length_and_renames_variants(
    run_workflow, run_dir, window_vcf, utils_docker
):
    # Tandem repeat loci with same-length alleles never reach PreprocessVcfs, which fails on them
    snv_ins_del = write_vcf(
        window_vcf, run_dir / "snv_ins_del.vcf.gz", lambda r: all(is_snv_ins_or_del(r.ref, a) for a in r.alts)
    )
    stripped = rewrite_vcf(
        snv_ins_del, run_dir / "input.vcf.gz", ["allele_type", "allele_length", "allele_subtype", "SOURCE"], lambda r: None
    )
    swaps = write_table(run_dir / "swap.tsv", [["HG00096", "SWAPPED"]])
    samples = ["SWAPPED", *SAMPLES[1:]]

    outputs = run_workflow(
        f"{UTILS}/PreprocessVcfs.wdl",
        {
            "vcfs": [str(stripped)],
            "vcf_idxs": [f"{stripped}.tbi"],
            "prefix": "pre",
            "normalize_vcfs": [True],
            "convert_aou_phase2_symbolic_alleles": [False],
            "source_tags": ["HPRC"],
            "swap_sample_lists": [str(swaps)],
            "min_length_cutoffs": [5],
            "max_length_cutoffs": [100],
            "sample_ids": samples,
            "ref_fa": str(WINDOW_FASTA),
            "ref_fai": f"{WINDOW_FASTA}.fai",
            "utils_docker": utils_docker,
        },
    )

    _, original = load_records(stripped)
    got_samples, got = load_records(outputs["PreprocessVcfs.preprocessed_vcf"])
    assert sorted(got_samples) == sorted(samples)
    assert any(len(r["alts"]) > 1 for r in original)
    assert len(got) == sum(len(r["alts"]) for r in original)

    def kind_and_length(ref, alt):
        return ("snv" if len(alt) == len(ref) else "ins" if len(alt) > len(ref) else "del"), len(alt) - len(ref)

    expected = collections.Counter(kind_and_length(r["ref"], a) for r in original for a in r["alts"])
    assert collections.Counter((r["info"]["allele_type"], r["info"]["allele_length"]) for r in got) == expected

    def base_name(r):
        kind, length = r["info"]["allele_type"], abs(r["info"]["allele_length"])
        if kind == "snv":
            return f"{CONTIG}-{r['pos']}-{r['ref']}-{r['alts'][0]}"
        return f"{CONTIG}-{r['pos']}-{kind.upper()}-{length}"

    names = collections.Counter(base_name(r) for r in got)
    expected_ids = sorted(f"{n}_{i}" if c > 1 else n for n, c in names.items() for i in range(1, c + 1))
    assert sorted(r["id"] for r in got) == expected_ids
    for r in got:
        assert r["info"]["SOURCE"] == "HPRC"
        length = abs(r["info"]["allele_length"])
        assert ("SMALL_HPRC" in r["filter"]) == (length < 5) and ("LARGE_HPRC" in r["filter"]) == (length > 100)

    # Biallelic calls keep their genotypes under the swapped sample name, unless left-alignment moved the record
    by_site = {site_key(r): r for r in got}
    kept = [(r, by_site[site_key(r)]) for r in original if len(r["alts"]) == 1 and site_key(r) in by_site]
    assert len(kept) > 20
    for old, new in kept:
        assert new["gts"]["SWAPPED"] == old["gts"]["HG00096"] and new["gts"]["HG00128"] == old["gts"]["HG00128"]


def test_preprocess_vcfs_drops_out_of_range_lengths_when_asked(run_workflow, run_dir, window_vcf, utils_docker):
    snv_ins_del = write_vcf(
        window_vcf, run_dir / "snv_ins_del.vcf.gz", lambda r: all(is_snv_ins_or_del(r.ref, a) for a in r.alts)
    )
    stripped = rewrite_vcf(
        snv_ins_del, run_dir / "input.vcf.gz", ["allele_type", "allele_length", "allele_subtype", "SOURCE"], lambda r: None
    )

    outputs = run_workflow(
        f"{UTILS}/PreprocessVcfs.wdl",
        {
            "vcfs": [str(stripped)],
            "vcf_idxs": [f"{stripped}.tbi"],
            "prefix": "pre",
            "normalize_vcfs": [True],
            "convert_aou_phase2_symbolic_alleles": [False],
            "source_tags": ["HPRC"],
            "swap_sample_lists": [],
            "min_length_cutoffs": [5],
            "max_length_cutoffs": [100],
            "drop_out_of_range_lengths": True,
            "sample_ids": [],
            "ref_fa": str(WINDOW_FASTA),
            "ref_fai": f"{WINDOW_FASTA}.fai",
            "utils_docker": utils_docker,
        },
    )

    _, original = load_records(stripped)
    _, got = load_records(outputs["PreprocessVcfs.preprocessed_vcf"])
    in_range = [a for r in original for a in r["alts"] if 5 <= abs(len(a) - len(r["ref"])) <= 100]
    assert in_range and len(got) == len(in_range)
    for r in got:
        assert 5 <= abs(r["info"]["allele_length"]) <= 100
        assert not {"SMALL_HPRC", "LARGE_HPRC"} & set(r["filter"])


def write_window_type_vcf(path, rows, sv_info=True):
    """Write one-sample (pos, ref, alt, info text, GT) rows on the window FASTA contig, bgzipped and tabix-indexed."""
    with pysam.FastaFile(str(WINDOW_FASTA)) as fasta:
        length = fasta.get_reference_length(CONTIG)
    lines = ["##fileformat=VCFv4.2", f"##contig=<ID={CONTIG},length={length}>"]
    if sv_info:
        lines += [
            '##INFO=<ID=SVTYPE,Number=1,Type=String,Description="SVTYPE">',
            '##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="SVLEN">',
        ]
    lines += [
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1",
    ]
    lines += [f"{CONTIG}\t{pos}\t.\t{ref}\t{alt}\t.\tPASS\t{info or '.'}\tGT\t{gt}" for pos, ref, alt, info, gt in rows]
    path.write_text("\n".join(lines) + "\n")
    return index_table(path, "vcf")


def run_preprocess_on_window_types(run_workflow, path, utils_docker, convert_symbolic=True):
    return run_workflow(
        f"{UTILS}/PreprocessVcfs.wdl",
        {
            "vcfs": [str(path)],
            "vcf_idxs": [f"{path}.tbi"],
            "prefix": "types",
            "normalize_vcfs": [True],
            "convert_aou_phase2_symbolic_alleles": [convert_symbolic],
            "source_tags": [],
            "swap_sample_lists": [],
            "min_length_cutoffs": [],
            "max_length_cutoffs": [],
            "sample_ids": [],
            "ref_fa": str(WINDOW_FASTA),
            "ref_fai": f"{WINDOW_FASTA}.fai",
            "utils_docker": utils_docker,
        },
    )


def other_base(base):
    return "C" if base != "C" else "G"


def reverse_complement(sequence):
    return sequence[::-1].translate(str.maketrans("ACGTN", "TGCAN"))


def test_preprocess_vcfs_keeps_only_snv_ins_del_records(run_workflow, run_dir, utils_docker):
    with pysam.FastaFile(str(WINDOW_FASTA)) as fasta:
        seq = fasta.fetch(CONTIG).upper()

    def ref(pos, length=1):
        return seq[pos - 1 : pos - 1 + length]

    padded = ref(10000, 3)
    path = write_window_type_vcf(
        run_dir / "types.vcf",
        [
            (1000, ref(1000), other_base(ref(1000)), None, "0/1"),
            (2000, ref(2000, 6), ref(2000), None, "0/1"),
            (3000, ref(3000), ref(3000) + "ACGTACGT", None, "0/1"),
            # Long-SV representation from MergeCohortLongSVCallsets: REF=N with a symbolic ALT
            (4000, "N", "<DEL>", "SVTYPE=DEL;SVLEN=-300", "0/1"),
            (5000, "N", "<DUP>", "SVTYPE=DUP;SVLEN=200", "0/1"),
            (6000, "N", "<INV>", "SVTYPE=INV;SVLEN=500", "0/1"),
            (7000, ref(7000), "<INS>", "SVTYPE=INS;SVLEN=400", "0/1"),
            (8000, ref(8000), ref(8000) + f"[{CONTIG}:9500[", "SVTYPE=BND", "0/1"),
            # Sequence-resolved inversion, as FixSymbolicRecords writes it
            (9000, ref(9000, 101), ref(9000) + reverse_complement(ref(9001, 100)), "SVTYPE=INV", "0/1"),
            # Padded SNV beside a deletion at a multiallelic site, and an SNV beside a spanning deletion
            (10000, padded, f"{padded[:2]}{other_base(padded[2])},{padded[0]}", None, "1/2"),
            (11000, ref(11000), f"{other_base(ref(11000))},*", None, "1/2"),
        ],
    )

    outputs = run_preprocess_on_window_types(run_workflow, path, utils_docker)

    _, got = load_records(outputs["PreprocessVcfs.preprocessed_vcf"])
    assert collections.Counter((r["info"]["allele_type"], r["info"]["allele_length"]) for r in got) == collections.Counter(
        [("snv", 0), ("del", -5), ("ins", 8), ("del", -300), ("ins", 200), ("snv", 0), ("del", -2), ("snv", 0)]
    )
    for r in got:
        assert r["ref"] == ref(r["pos"], len(r["ref"]))
        assert not r["alts"][0].startswith("<") and r["alts"][0] != "*"
    assert {r["pos"] for r in got if r["info"]["allele_type"] == "snv"} == {1000, 10002, 11000}


def test_preprocess_vcfs_keeps_uppercased_records_without_an_svtype_header(run_workflow, run_dir, utils_docker):
    # DeepVariant VCFs declare no SVTYPE
    with pysam.FastaFile(str(WINDOW_FASTA)) as fasta:
        seq = fasta.fetch(CONTIG).upper()
    path = write_window_type_vcf(
        run_dir / "dv.vcf",
        # Lowercase bases, as normalizing against a soft-masked reference writes them
        [
            (1000, seq[999].lower(), other_base(seq[999]).lower(), None, "0/1"),
            (3000, seq[2999], seq[2999] + "acgtacgt", None, "0/1"),
        ],
        sv_info=False,
    )

    outputs = run_preprocess_on_window_types(run_workflow, path, utils_docker, convert_symbolic=False)

    _, got = load_records(outputs["PreprocessVcfs.preprocessed_vcf"])
    assert [(r["info"]["allele_type"], r["info"]["allele_length"]) for r in got] == [("snv", 0), ("ins", 8)]
    assert all(r["ref"].isupper() and all(a.isupper() for a in r["alts"]) for r in got)


def test_preprocess_vcfs_rejects_same_length_multibase_records(run_workflow, run_dir, utils_docker):
    with pysam.FastaFile(str(WINDOW_FASTA)) as fasta:
        mnp_ref = fasta.fetch(CONTIG, 11999, 12002).upper()
    mnp_alt = other_base(mnp_ref[0]) + mnp_ref[1] + other_base(mnp_ref[2])
    path = write_window_type_vcf(run_dir / "mnp.vcf", [(12000, mnp_ref, mnp_alt, None, "0/1")])

    with pytest.raises(WorkflowError):
        run_preprocess_on_window_types(run_workflow, path, utils_docker)


def test_preprocess_vcfs_concatenates_record_shards_of_several_vcfs_with_source_tags(
    run_workflow, run_dir, window_vcf, utils_docker
):
    stripped = rewrite_vcf(
        window_vcf, run_dir / "input.vcf.gz", ["allele_type", "allele_length", "allele_subtype", "SOURCE"], lambda r: None
    )
    counter = itertools.count()
    first = write_vcf(stripped, run_dir / "first.vcf.gz", lambda r: next(counter) % 2 == 0)
    counter = itertools.count()
    second = write_vcf(stripped, run_dir / "second.vcf.gz", lambda r: next(counter) % 2 != 0)

    outputs = run_workflow(
        f"{UTILS}/PreprocessVcfs.wdl",
        {
            "vcfs": [str(first), str(second)],
            "vcf_idxs": [f"{first}.tbi", f"{second}.tbi"],
            "prefix": "pre",
            "normalize_vcfs": [True, True],
            "convert_aou_phase2_symbolic_alleles": [False, False],
            "source_tags": ["A", "B"],
            "swap_sample_lists": [],
            "min_length_cutoffs": [],
            "max_length_cutoffs": [],
            "records_per_shard": 20,
            "sample_ids": [],
            "ref_fa": str(WINDOW_FASTA),
            "ref_fai": f"{WINDOW_FASTA}.fai",
            "utils_docker": utils_docker,
        },
    )

    halves = {"A": load_records(first), "B": load_records(second)}
    got_samples, got = load_records(outputs["PreprocessVcfs.preprocessed_vcf"])
    expected = {tag: sum(len(r["alts"]) for r in records) for tag, (_, records) in halves.items()}
    assert all(expected.values())
    assert len(got) == sum(expected.values())
    assert collections.Counter(r["info"]["SOURCE"] for r in got) == collections.Counter(expected)
    assert sorted(got_samples) == sorted(halves["A"][0]) == sorted(halves["B"][0])


def test_merge_vcfs_merges_overlapping_callsets_over_disjoint_samples_and_records_provenance(
    run_workflow, run_dir, window_vcf, utils_docker
):
    # Biallelic non-TR records, which Truvari can compare, and every TR record
    _, loaded = load_records(window_vcf)
    usable = [r for r in loaded if r["info"]["allele_type"] == "trv" or len(r["alts"]) == 1]
    assert len(usable) > 40
    first, second = usable[: len(usable) * 2 // 3], usable[len(usable) // 3 :]
    shared = {site_key(r) for r in first} & {site_key(r) for r in second}
    assert len(shared) > 10
    callsets = {"A": (first, SAMPLES[:3]), "B": (second, SAMPLES[3:])}

    paths = {}
    for name, (records, samples) in callsets.items():
        keep = {site_key(r) for r in records}
        full = write_vcf(
            window_vcf, run_dir / f"{name}.all.vcf.gz", lambda r, keep=keep: (r.contig, r.pos, r.ref, r.alts) in keep
        )
        plain = run_dir / f"{name}.vcf.gz"
        bcftools.view("-s", ",".join(samples), "-Oz", "-o", str(plain), str(full), catch_stdout=False)
        pysam.tabix_index(str(plain), preset="vcf", force=True)
        paths[name] = plain

    outputs = run_workflow(
        f"{UTILS}/MergeVcfs.wdl",
        {
            "vcfs": [str(paths["A"]), str(paths["B"])],
            "vcf_idxs": [f"{paths['A']}.tbi", f"{paths['B']}.tbi"],
            "vcf_names": ["A", "B"],
            "contig": CONTIG,
            "prefix": "merged",
            "min_truvari_match": 50,
            "truvari_breakpoint_window": 500,
            "truvari_reciprocal_overlap": 0.9,
            "truvari_sample_similarity": 0.0,
            "truvari_sequence_similarity": 0.9,
            "truvari_size_similarity": 0.9,
            "truvari_size_max": 100000,
            "truvari_size_min": 50,
            "ref_fa": str(WINDOW_FASTA),
            "ref_fai": f"{WINDOW_FASTA}.fai",
            "shard_bin_size": 30000,
            "utils_docker": utils_docker,
        },
    )

    got_samples, got = load_records(outputs["MergeVcfs.merged_vcf"])
    assert sorted(got_samples) == sorted(SAMPLES)
    by_key = {site_key(r): r for r in got}
    union = {site_key(r) for r in first} | {site_key(r) for r in second}
    assert set(by_key) == union and len(got) == len(union)
    by_original = {site_key(r): r for r in usable}
    for key, record in by_key.items():
        in_a, in_b = key in {site_key(r) for r in first}, key in {site_key(r) for r in second}
        assert record["info"]["MERGE_COUNT"] == in_a + in_b
        assert list(record["info"]["SOURCE_NAMES"]) == [n for n, present in (("A", in_a), ("B", in_b)) if present]
        for sample in SAMPLES[:3]:
            assert record["gts"][sample] == (by_original[key]["gts"][sample] if in_a else (None, None))
        for sample in SAMPLES[3:]:
            assert record["gts"][sample] == (by_original[key]["gts"][sample] if in_b else (None, None))
    assert sum(r["info"]["MERGE_COUNT"] == 2 for r in got) == len(shared)

    rows = read_table(outputs["MergeVcfs.merge_summary_tsv"])
    assert len(rows) > 2 and any("A" in cell for cell in rows[0]) and any("B" in cell for cell in rows[0])


def test_annotate_allele_type_overwrites_types_and_sets_origin_and_sub_family(run_workflow, run_dir, vcf_inputs, source):
    records = source[1]
    key = lambda r: [r["contig"], r["pos"], r["ref"], r["alts"][0], r["id"]]  # noqa: E731
    dups = [r for r in records if r["info"]["allele_type"] == "dup" and len(r["alts"]) == 1][:8]
    inserts = [r for r in records if r["info"]["allele_type"] == "ins" and len(r["alts"]) == 1][:8]
    deletions = [r for r in records if r["info"]["allele_type"] == "del" and len(r["alts"]) == 1][:8]
    # A dup with no source locus falls back to the final column
    dup_rows = [key(r) + ["tandem_dup", "" if i == 0 else f"chr22:{i}-{i + 99}_+", "fallback_origin"] for i, r in enumerate(dups)]
    mei_rows = [key(r) + ["Alu_ins", f"AluY{i}"] for i, r in enumerate(inserts)]
    med_rows = [key(r) + ["Alu_del", f"AluS{i}"] for i, r in enumerate(deletions)]

    outputs = run_workflow(
        f"{UTILS}/AnnotateAlleleType.wdl",
        {
            **vcf_inputs,
            "contigs": [CONTIG],
            "prefix": "types",
            "dup_tsv": str(write_table(run_dir / "dup.tsv", dup_rows)),
            "mei_tsv": str(write_table(run_dir / "mei.tsv", mei_rows)),
            "med_tsv": str(write_table(run_dir / "med.tsv", med_rows)),
            "dup_prefix": "p_",
            "dup_suffix": "_s",
            "dup_lowercase": True,
            "mei_prefix": "mei_",
            "med_lowercase": True,
        },
    )

    _, got = load_records(outputs["AnnotateAlleleType.allele_type_annotated_vcf"])
    assert [site_key(r) for r in got] == [site_key(r) for r in records]
    expected_type, expected_origin, expected_family = {}, {}, {}
    for i, r in enumerate(dups):
        expected_type[r["id"]] = "p_tandem_dup_s"
        expected_origin[r["id"]] = f"chr22:{i}-{i + 99}_+" if i else "fallback_origin"
    for i, r in enumerate(inserts):
        expected_type[r["id"]], expected_family[r["id"]] = "mei_Alu_ins", f"AluY{i}"
    for i, r in enumerate(deletions):
        expected_type[r["id"]], expected_family[r["id"]] = "alu_del", f"AluS{i}"
    for new, old in zip(got, records):
        assert new["info"]["allele_type"] == expected_type.get(old["id"], old["info"]["allele_type"])
        if old["id"] in expected_origin:
            assert new["info"]["ORIGIN"] == (expected_origin[old["id"]],)
        assert new["info"].get("SUB_FAMILY") == expected_family.get(old["id"], old["info"].get("SUB_FAMILY"))
        assert new["gts"] == old["gts"]
    assert len(dups) > 2 and len(inserts) > 2 and len(deletions) > 2


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
            "sv_source_tag": "HPRC_SV_Integration",
            "long_sv_source_tag": "HPRC_LongSV_Integration",
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
        if abs(r["info"]["allele_length"]) >= min_sv_length
        and r["info"]["SOURCE"] in ("HPRC_SV_Integration", "HPRC_LongSV_Integration")
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


def test_convert_vcf_to_bed_concatenates_shards_under_one_header(run_workflow, run_dir, vcf_inputs):
    info_columns = ["allele_type", "allele_length"]
    inputs = {
        "vcfs": [vcf_inputs["vcf"]],
        "vcf_idxs": [vcf_inputs["vcf_idx"]],
        "contigs": [CONTIG],
        "prefix": "converted",
        "records_per_shard": 100,
        "info_columns": info_columns,
        "utils_docker": vcf_inputs["utils_docker"],
    }
    outputs = run_workflow(f"{UTILS}/ConvertVcfToBed.wdl", inputs)

    expected = run_dir / "expected.bed"
    info_args = [arg for column in info_columns for arg in ("--info", column)]
    subprocess.run(
        [
            sys.executable,
            REPO_ROOT / "scripts/qc/vcf_qc/vcf2bed.py",
            *info_args,
            "--include-filters",
            vcf_inputs["vcf"],
            expected,
        ],
        check=True,
    )
    lines = open(outputs["ConvertVcfToBed.bed"]).read().splitlines()
    assert lines == expected.read_text().splitlines()
    assert lines[0] == "#chrom\tstart\tend\tname\tsvtype\tsamples\tallele_type\tallele_length\tFILTER"
    assert len(lines) == 386


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


def test_summarize_annotations_gatksv_splits_alleles_and_draws_the_cohort_plots(run_workflow, run_dir, window_vcf, utils_docker):
    ped = write_table(run_dir / "cohort.ped", [[s, s, 0, 0, 1 + i % 2, 0] for i, s in enumerate(SAMPLES)])
    outputs = run_workflow(
        f"{UTILS}/SummarizeAnnotationsGATKSV.wdl",
        {
            "vcfs": [str(window_vcf)],
            "vcf_idxs": [f"{window_vcf}.tbi"],
            "prefix": "qc",
            "ped": str(ped),
            "ref_fa": str(WINDOW_FASTA),
            "ref_fai": f"{WINDOW_FASTA}.fai",
            "records_per_shard": 30,
            "samples_per_shard": 2,
            "utils_docker": utils_docker,
        },
    )

    # Multi-allelic records are split into one row per ALT, named <ID>_<n>, so carriers are compared per original ID
    _, records = load_records(window_vcf)
    expected, got = collections.defaultdict(collections.Counter), collections.defaultdict(collections.Counter)
    for record in records:
        for allele in range(1, len(record["alts"]) + 1):
            expected[record["id"]][tuple(sorted(s for s, gt in record["gts"].items() if allele in gt))] += 1
    for row in read_dicts(outputs["SummarizeAnnotationsGATKSV.vcf2bed_bed"]):
        name = row["name"] if row["name"] in expected else row["name"].rsplit("_", 1)[0]
        got[name][tuple(sorted(filter(None, row["samples"].split(","))))] += 1
    assert got == expected

    with tarfile.open(outputs["SummarizeAnnotationsGATKSV.vcf_qc_tarball"]) as tar:
        names = {Path(n).name for n in tar.getnames()}
        member = next(m for m in tar.getmembers() if m.name.endswith("qc.samples_analyzed.list"))
        analyzed = tar.extractfile(member).read().decode().split()
    assert analyzed == SAMPLES
    assert {"qc.VCF_sites.stats.bed.gz", "qc.VCF_sites.stats.bed.gz.tbi"} <= names
    plots = {"size_distribution.pdf", "af_distribution.pdf", "gt_overall_distribution.png", "counts_per_genome_distribution.png"}
    assert plots <= names
