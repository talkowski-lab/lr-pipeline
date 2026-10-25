import collections
import datetime
import gzip
import itertools
import json

import pysam
import pytest
from conftest import (
    ANNOTATION,
    CONTIG,
    SAMPLES,
    load_records,
    read_table,
    rewrite_vcf,
    site_key,
    write_sites_vcf,
    write_table,
    write_vcf,
)

pytestmark = pytest.mark.wdl


def build_overlap_case(run_dir, subset_vcf):
    counter = itertools.count()
    # Truth lacks every fourth record and holds every seventh remaining one three bases away, so neither matches exactly
    shifted = itertools.count()
    shifted_sites = set()

    def shift(record):
        if next(shifted) % 7 == 0:
            shifted_sites.add((record.contig, record.pos, record.ref, record.alts))
            record.pos += 3

    kept = write_vcf(subset_vcf, run_dir / "kept.vcf.gz", lambda r: next(counter) % 4 != 0)
    truth = rewrite_vcf(kept, run_dir / "truth.vcf.gz", [], shift)
    return truth, load_records(truth)[1], shifted_sites


def overlap_inputs(subset_vcf, truth, utils_docker, **overrides):
    return {
        "vcf": str(subset_vcf),
        "vcf_idx": f"{subset_vcf}.tbi",
        "truth_snv_indel_vcf": str(truth),
        "truth_snv_indel_vcf_idx": f"{truth}.tbi",
        "truth_sv_vcf": str(truth),
        "truth_sv_vcf_idx": f"{truth}.tbi",
        "contig": CONTIG,
        "prefix": "overlap",
        "run_exact_matching": True,
        "run_truvari_matching": False,
        "run_bedtools_closest_matching": False,
        "min_sv_length_truvari_vcf": 50,
        "min_sv_length_truvari_truth_snv_indel_vcf": 50,
        "min_sv_length_bedtools_closest_vcf": 50,
        "min_sv_length_bedtools_closest_truth_vcf": 50,
        "convert_symbolic_truth_sv_vcf": False,
        "move_dup_to_origin": False,
        "subset_contig_truth_sv_vcf": False,
        "shard_bin_size_exact_match": 15000000,
        "args_string_vcf": "N_ALT==1",
        "args_string_truth_snv_indel_vcf": "N_ALT==1",
        "utils_docker": utils_docker,
        **overrides,
    }


def overlap_rows(outputs):
    header = [line[0] for line in read_table(outputs["AnnotateCallsetOverlap.annotations_header_benchmark"])]
    assert header[:4] == ["match_type", "truth_ID", "source_tag", "filter"]
    assert header[-3:] == ["N_HOMREF", "N_HET", "N_HOMALT"]
    return [
        dict(zip(["CHROM", "POS", "REF", "ALT", "ID"] + header, row))
        for row in read_table(outputs["AnnotateCallsetOverlap.annotations_tsv_benchmark"])
    ]


def row_site(row):
    return (row["CHROM"], int(row["POS"]), row["REF"], (row["ALT"],))


def assert_truth_fields(row, record, match_type):
    assert row["match_type"] == match_type
    assert row["source_tag"] == "SNV_indel"
    assert row["truth_ID"] == record["id"]
    # Genotype counts are derived from the truth's INFO allele counts: het = AC - 2 * nhomalt, homref = AN / 2 - het - nhomalt
    info = record["info"]
    het = info["AC"][0] - 2 * info["nhomalt"]
    assert (int(row["N_HOMREF"]), int(row["N_HET"]), int(row["N_HOMALT"])) == (
        int(info["AN"] / 2 - info["nhomalt"] - het),
        het,
        info["nhomalt"],
    )
    assert int(row["AC"]) == info["AC"][0]


def test_annotate_callset_overlap_exact_matches_biallelic_records_and_carries_truth_fields(
    run_workflow, run_dir, subset_vcf, utils_docker, source
):
    truth, truth_records, _ = build_overlap_case(run_dir, subset_vcf)
    outputs = run_workflow(f"{ANNOTATION}/AnnotateCallsetOverlap.wdl", overlap_inputs(subset_vcf, truth, utils_docker))

    truth_by_site = {site_key(r): r for r in truth_records}
    expected = {site_key(r): r for r in source[1] if len(r["alts"]) == 1 and site_key(r) in truth_by_site}
    assert 0 < len(expected) < sum(len(r["alts"]) == 1 for r in source[1])

    rows = overlap_rows(outputs)
    assert {row_site(r) for r in rows} == set(expected)
    for row in rows:
        assert_truth_fields(row, truth_by_site[row_site(row)], "EXACT")


def test_annotate_callset_overlap_truvari_matches_shifted_svs_across_record_shards(
    run_workflow, run_dir, subset_vcf, utils_docker, source
):
    truth, truth_records, shifted_sites = build_overlap_case(run_dir, subset_vcf)
    outputs = run_workflow(
        f"{ANNOTATION}/AnnotateCallsetOverlap.wdl",
        overlap_inputs(subset_vcf, truth, utils_docker, run_truvari_matching=True, shard_records_truvari_match=10),
    )

    truth_by_site = {site_key(r): r for r in truth_records}
    exact = {site_key(r) for r in source[1] if len(r["alts"]) == 1 and site_key(r) in truth_by_site}
    shifted_svs = {
        site_key(r)
        for r in source[1]
        if len(r["alts"]) == 1 and abs(r["info"].get("allele_length") or 0) >= 50 and site_key(r) in shifted_sites
    }
    assert shifted_svs
    rows = overlap_rows(outputs)
    by_site = {row_site(r): r for r in rows}
    assert exact <= set(by_site)
    for site in exact:
        assert by_site[site]["match_type"] == "EXACT"
        assert_truth_fields(by_site[site], truth_by_site[site], "EXACT")

    # Each shifted SV matches its own truth record through Truvari, whichever record shard held it
    source_by_site = {site_key(r): r for r in source[1]}
    id_to_truth = {r["id"]: r for r in truth_records}
    for site in shifted_svs:
        record = source_by_site[site]
        assert site in by_site
        assert_truth_fields(by_site[site], id_to_truth[record["id"]], "TRUVARI_0.9")
    for site, row in by_site.items():
        if site not in exact:
            assert row["match_type"].startswith("TRUVARI_")


POPS = {"HG00096": "nfe", "HG00097": "nfe", "HG00099": "afr", "HG00126": "eas", "HG00128": "afr"}
SEXES = {"HG00096": 1, "HG00097": 2, "HG00099": 1, "HG00126": 2, "HG00128": 1}


def write_af_inputs(run_dir):
    pops = write_table(run_dir / "pops.tsv", [[s, p] for s, p in POPS.items()])
    ped = write_table(run_dir / "samples.ped", [[f"fam{i}", s, 0, 0, SEXES[s], 0] for i, s in enumerate(SAMPLES)])
    return pops, ped


def allele_stats(genotypes, hemi=False):
    """AN, AC, AF, NCR and genotype counts of one biallelic record over genotypes, as compute_AFs defines them."""
    called = [[a for a in gt if a is not None] for gt in genotypes]
    an = sum(len(gt) for gt in called)
    ac = sum(a > 0 for gt in called for a in gt)
    homref = sum(gt == (0, 0) for gt in genotypes)
    het = sum(sorted(g for g in gt if g is not None) == [0, 1] for gt in genotypes)
    homalt = sum(len([a for a in gt if a is not None and a > 0]) == 2 for gt in genotypes)
    if hemi:
        an //= 2
        ac = het + homalt
    ncr = round(1 - an / (len(genotypes) * (1 if hemi else 2)), 6) if genotypes else 1.0
    return {
        "AN": an,
        "AC": ac,
        "AF": round(ac / an, 6) if an else 0,
        "NCR": ncr,
        "nhomref": homref,
        "nhet": het,
        "nhomalt": homalt,
    }


def test_annotate_af_computes_global_population_and_sex_frequencies(run_workflow, run_dir, vcf_inputs, source):
    pops, ped = write_af_inputs(run_dir)
    par_bed = write_table(run_dir / "par.bed", [["chrX", 10001, 2781479]])
    outputs = run_workflow(
        f"{ANNOTATION}/AnnotateAF.wdl",
        {
            **vcf_inputs,
            "contigs": [CONTIG],
            "prefix": "af",
            "records_per_shard": 150,
            "sample_pop_assignments": str(pops),
            "ped": str(ped),
            "par_bed": str(par_bed),
        },
    )

    _, got = load_records(outputs["AnnotateAF.af_annotated_vcf"])
    assert [site_key(r) for r in got] == [site_key(r) for r in source[1]]
    groups = {"": SAMPLES, "_XY": [s for s in SAMPLES if SEXES[s] == 1], "_XX": [s for s in SAMPLES if SEXES[s] == 2]}
    for pop in set(POPS.values()):
        in_pop = [s for s in SAMPLES if POPS[s] == pop]
        groups[f"_{pop}"] = in_pop
        groups[f"_{pop}_XY"] = [s for s in in_pop if SEXES[s] == 1]
        groups[f"_{pop}_XX"] = [s for s in in_pop if SEXES[s] == 2]
    biallelic = 0
    for new, old in zip(got, source[1]):
        assert "par" not in new["info"]
        if len(old["alts"]) != 1:
            # Multiallelic sites carry per-alt counts and the genotyped allele total
            alleles = [a for s in SAMPLES for a in old["gts"][s] if a is not None]
            assert new["info"]["AN"] == len(alleles)
            assert new["info"]["AC"] == tuple(alleles.count(i) for i in range(1, len(old["alts"]) + 1))
            continue
        biallelic += 1
        for suffix, members in groups.items():
            expected = allele_stats([old["gts"][s] for s in members])
            for key, value in expected.items():
                got_value = new["info"][f"{key}{suffix}"]
                got_value = got_value[0] if isinstance(got_value, tuple) else got_value
                assert got_value == pytest.approx(value, abs=1e-6), (old["id"], key, suffix)
        afs = {pop: new["info"][f"AF_{pop}"][0] for pop in set(POPS.values())}
        assert new["info"]["AF_grpmax"] == pytest.approx(max(afs.values()))
        assert afs[new["info"]["grpmax"]] == max(afs.values())
    assert biallelic > 200


def test_annotate_af_uses_hemizygous_counts_for_males_outside_par(run_workflow, run_dir, subset_vcf, utils_docker, source):
    # Relabel the fixture as chrX, with the first half of its records inside a pseudoautosomal region
    plain = run_dir / "chrx.vcf"
    with gzip.open(subset_vcf, "rt") as handle, open(plain, "w") as out:
        out.writelines(line.replace("chr22", "chrX") for line in handle)
    pysam.tabix_index(str(plain), preset="vcf", force=True)
    chrx = run_dir / "chrx.vcf.gz"
    _, records = load_records(chrx)
    par_end = records[len(records) // 2]["pos"]
    par_bed = write_table(run_dir / "par.bed", [["chrX", 0, par_end]])
    pops, ped = write_af_inputs(run_dir)

    outputs = run_workflow(
        f"{ANNOTATION}/AnnotateAF.wdl",
        {
            "vcf": str(chrx),
            "vcf_idx": f"{chrx}.tbi",
            "contigs": ["chrX"],
            "prefix": "afx",
            "records_per_shard": 1000,
            "ped": str(ped),
            "par_bed": str(par_bed),
            "utils_docker": utils_docker,
        },
    )

    _, got = load_records(outputs["AnnotateAF.af_annotated_vcf"])
    males = [s for s in SAMPLES if SEXES[s] == 1]
    females = [s for s in SAMPLES if SEXES[s] == 2]
    in_par = outside = 0
    for new, old in zip(got, records):
        if len(old["alts"]) != 1:
            continue
        is_par = old["pos"] <= par_end
        assert ("par" in new["info"]) == is_par
        male = allele_stats([old["gts"][s] for s in males], hemi=not is_par)
        female = allele_stats([old["gts"][s] for s in females])
        assert new["info"]["AN_XY"] == male["AN"] and new["info"]["AC_XY"] == (male["AC"],)
        assert new["info"]["AN_XX"] == female["AN"] and new["info"]["AC_XX"] == (female["AC"],)
        if is_par:
            whole = allele_stats([old["gts"][s] for s in SAMPLES])
            assert new["info"]["AN"] == whole["AN"] and new["info"]["AC"] == (whole["AC"],)
            in_par += 1
        else:
            # Outside the PAR the total is rebuilt from the hemizygous male and diploid female counts
            assert new["info"]["AN"] == male["AN"] + female["AN"]
            assert new["info"]["AC"] == (male["AC"] + female["AC"],)
            assert new["info"]["nhemialt_XY"] == male["nhet"] + male["nhomalt"]
            outside += 1
    assert in_par > 50 and outside > 50


def test_annotate_dbsnp_reports_rsids_for_variants_matching_on_position_and_alleles(run_workflow, run_dir, vcf_inputs, source):
    snvs = [r for r in source[1] if r["info"]["allele_type"] == "snv"]
    flipped = {"A": "C", "C": "G", "G": "T", "T": "A"}
    rows, expected = [], {}
    for i, r in enumerate(snvs):
        if i % 3 == 0:
            rows.append((r["pos"], f"rs{r['pos']}", r["ref"], r["alts"][0], ""))
            expected[r["id"]] = f"rs{r['pos']}"
        elif i % 3 == 1:
            # Same position with a different alternate allele must not match
            wrong = flipped[r["alts"][0]] if flipped[r["alts"][0]] != r["ref"] else flipped[flipped[r["alts"][0]]]
            rows.append((r["pos"], f"rs_wrong_{r['pos']}", r["ref"], wrong, ""))
    dbsnp = write_sites_vcf(run_dir / "dbsnp.vcf", rows)

    outputs = run_workflow(
        f"{ANNOTATION}/AnnotateDbSNP.wdl",
        {
            **vcf_inputs,
            "dbsnp_vcf": str(dbsnp),
            "dbsnp_vcf_idx": f"{dbsnp}.tbi",
            "contigs": [CONTIG],
            "prefix": "dbsnp",
            "shard_bin_size": 15000000,
        },
    )

    by_id = {r["id"]: r for r in snvs}
    got = {row[4]: row for row in read_table(outputs["AnnotateDbSNP.annotations_tsv_dbsnp"])}
    assert {k: v[5] for k, v in got.items()} == expected
    for rid, row in got.items():
        assert row[:4] == [CONTIG, str(by_id[rid]["pos"]), by_id[rid]["ref"], by_id[rid]["alts"][0]]


def test_annotate_dbvar_matches_deletions_and_insertions_within_breakpoint_windows(run_workflow, run_dir, vcf_inputs, source):
    min_length = 50
    eligible = [
        r for r in source[1] if r["info"]["allele_type"] != "trv" and abs(r["info"].get("allele_length", 0)) >= min_length
    ]
    # Keep only calls with no other eligible call nearby, so that each dbVar record can match exactly one of them
    candidates = [
        r
        for r in eligible
        if r["info"]["allele_type"] in ("del", "ins")
        and len(r["alts"]) == 1
        and not any(o is not r and abs(o["pos"] - r["pos"]) < 1500 for o in eligible)
    ]
    rows, expected = [], {}
    for i, r in enumerate(candidates):
        length = abs(r["info"]["allele_length"])
        kind = r["info"]["allele_type"].upper()
        end = r["pos"] + length if kind == "DEL" else r["pos"] + 1
        if i % 2 == 0:
            rows.append((r["pos"] + 3, f"dbvar_{i}", "N", f"<{kind}>", f"SVTYPE={kind};SVLEN={length};END={end + 3}"))
            expected[r["id"]] = f"dbvar_{i}"
        elif i % 4 == 1:
            # Far outside every breakpoint window
            rows.append((r["pos"] + 4000, f"dbvar_far_{i}", "N", f"<{kind}>", f"SVTYPE={kind};SVLEN={length};END={end + 4000}"))
    assert len(expected) > 10
    dbvar = write_sites_vcf(
        run_dir / "dbvar.vcf", sorted(rows), info=[("SVTYPE", "1", "String"), ("SVLEN", "1", "Integer"), ("END", "1", "Integer")]
    )

    outputs = run_workflow(
        f"{ANNOTATION}/AnnotateDbVaR.wdl",
        {
            **vcf_inputs,
            "dbvar_vcf": str(dbvar),
            "dbvar_vcf_idx": f"{dbvar}.tbi",
            "contigs": [CONTIG],
            "prefix": "dbvar",
            "min_length": min_length,
            "records_per_shard": 100,
        },
    )

    by_id = {r["id"]: r for r in candidates}
    got = {row[4]: row for row in read_table(outputs["AnnotateDbVaR.annotations_tsv_dbvar"])}
    assert {k: v[5] for k, v in got.items()} == expected
    for rid, row in got.items():
        record = by_id[rid]
        assert row[:4] == [CONTIG, str(record["pos"]), record["ref"], record["alts"][0]]


def test_annotate_gnomad_str_matches_enveloped_loci_and_motif_matched_overlaps(run_workflow, run_dir, vcf_inputs, source):
    trvs = [r for r in source[1] if r["info"]["allele_type"] == "trv"]
    catalog, expected = [], {}
    for i, r in enumerate(trvs):
        start, end = r["pos"] - 1, r["stop"]
        shift = max(1, (end - start) // 10)
        locus = f"locus_{i}"
        motif = r["info"]["MOTIFS"][0]
        if i % 4 == 0:
            # The locus envelops the variant, whatever its motif
            catalog.append({"LocusId": locus, "MainReferenceRegion": f"{CONTIG}:{start - 5}-{end + 5}", "RepeatUnit": "ZZZ"})
            expected[r["id"]] = locus
        elif i % 4 == 1:
            # A reciprocal overlap of at least 0.7 needs a matching motif
            catalog.append(
                {"LocusId": locus, "MainReferenceRegion": f"{CONTIG}:{start + shift}-{end + shift}", "RepeatUnit": motif}
            )
            expected[r["id"]] = locus
        elif i % 4 == 2:
            catalog.append(
                {"LocusId": locus, "MainReferenceRegion": f"{CONTIG}:{start + shift}-{end + shift}", "RepeatUnit": "ZZZ"}
            )
    catalog_json = run_dir / "catalog.json"
    catalog_json.write_text(json.dumps(catalog))

    outputs = run_workflow(
        f"{ANNOTATION}/AnnotateGnomADSTR.wdl",
        {**vcf_inputs, "contigs": [CONTIG], "prefix": "str", "records_per_shard": 40, "gnomad_tr_json": str(catalog_json)},
    )

    got = {row[4]: row for row in read_table(outputs["AnnotateGnomADSTR.annotations_tsv_gnomad_str"])}
    assert {k: v[5] for k, v in got.items()} == expected
    by_id = {r["id"]: r for r in trvs}
    for rid, row in got.items():
        assert row[:4] == [CONTIG, str(by_id[rid]["pos"]), by_id[rid]["ref"], ",".join(by_id[rid]["alts"])]


def test_annotate_region_labels_breakpoints_by_repeat_priority_and_large_calls_by_body_coverage(
    run_workflow, run_dir, vcf_inputs, source
):
    def normalized(record):
        kind = record["info"]["allele_type"].upper()
        return "SNV" if kind in ("SNV", "TRV") else "DEL" if "DEL" in kind else "INS" if "INS" in kind else "DUP"

    sr, sd, rm = [], [], []
    for i, r in enumerate(source[1]):
        left, right = r["pos"] - 1, r["stop"]
        if i % 5 == 0:
            sr.append((CONTIG, left - 5, left + 5))
        if i % 5 in (0, 2):
            rm.append((CONTIG, left - 5, left + 5))
        if i % 5 == 1:
            sd.append((CONTIG, right - 5, right + 5))
        if i % 5 == 3 and normalized(r) in ("DEL", "DUP") and right - left > 5000:
            sd.append((CONTIG, left, right))
    beds = [write_table(run_dir / f"{name}.bed", sorted(rows)) for name, rows in (("sr", sr), ("sd", sd), ("rm", rm))]

    outputs = run_workflow(
        f"{ANNOTATION}/AnnotateRegion.wdl",
        {
            **vcf_inputs,
            "contigs": [CONTIG],
            "prefix": "region",
            "records_per_shard": 150,
            "simple_repeats_bed": str(beds[0]),
            "seg_dup_bed": str(beds[1]),
            "repeat_masked_bed": str(beds[2]),
        },
    )

    def covered(intervals, point):
        return any(s <= point < e or s == point for _, s, e in intervals)

    got = {row[4]: row[5] for row in read_table(outputs["AnnotateRegion.annotations_tsv_region"])}
    assert set(got) == {r["id"] for r in source[1]}
    labels = collections.Counter()
    for r in source[1]:
        left, right = r["pos"] - 1, r["stop"]
        label = "US"
        for name, intervals in (("RM", rm), ("SD", sd), ("SR", sr)):
            if covered(intervals, left) or covered(intervals, right):
                label = name
        if normalized(r) in ("DEL", "DUP") and right - left > 5000:
            label = "SD" if any(s <= left and right <= e for _, s, e in sd) else "US"
        assert got[r["id"]] == label, (r["id"], label)
        labels[label] += 1
    assert {"US", "SR", "SD", "RM"} <= set(labels)


def test_annotate_age_metrics_bins_heterozygous_and_homozygous_carriers_by_age(run_workflow, run_dir, vcf_inputs, source):
    reference_date, bins = datetime.date(2026, 1, 1), [0, 30, 50, 70]
    birth_dates = {
        "HG00096": "1990-06-15T00:00:00+00:00",
        "HG00097": "2015-01-01T00:00:00+00:00",
        "HG00099": "1940-03-01T00:00:00+00:00",
        "HG00126": "1980-01-01T00:00:00+00:00",
        "HG00128": "2027-01-01T00:00:00+00:00",
    }
    ages = {s: (reference_date - datetime.date.fromisoformat(d[:10])).days / 365.25 for s, d in birth_dates.items()}
    age_data = run_dir / "ages.csv"
    age_data.write_text("person_id,date_of_birth\n" + "".join(f"{s},{d}\n" for s, d in birth_dates.items()))

    def bucket(age):
        if age < bins[0]:
            return "smaller"
        if age >= bins[-1]:
            return "larger"
        return max(i for i in range(len(bins) - 1) if bins[i] <= age)

    outputs = run_workflow(
        f"{ANNOTATION}/AnnotateAgeMetrics.wdl",
        {
            **vcf_inputs,
            "contigs": [CONTIG],
            "prefix": "age",
            "records_per_shard": 150,
            "age_data": str(age_data),
            "age_bins": bins,
            "reference_date": str(reference_date),
        },
    )

    rows = read_table(outputs["AnnotateAgeMetrics.annotations_tsv_age"])
    assert len(rows) == len(source[1])
    n_bins = len(bins) - 1
    seen = collections.Counter()
    for row, r in zip(rows, source[1]):
        assert row[:5] == [r["contig"], str(r["pos"]), r["ref"], ",".join(r["alts"]), r["id"]]
        for alt_index in range(len(r["alts"])):
            counts = {"het": [[0] * n_bins, 0, 0], "hom": [[0] * n_bins, 0, 0]}
            for sample in SAMPLES:
                gt = r["gts"][sample]
                if None in gt:
                    continue
                allele = alt_index + 1
                kind = (
                    "het"
                    if gt.count(allele) == 1 and gt.count(0) == len(gt) - 1
                    else "hom"
                    if gt.count(allele) == len(gt)
                    else None
                )
                if kind is None:
                    continue
                where = bucket(ages[sample])
                if where == "smaller":
                    counts[kind][1] += 1
                elif where == "larger":
                    counts[kind][2] += 1
                else:
                    counts[kind][0][where] += 1
                seen[(kind, where)] += 1
            for offset, kind in ((5, "het"), (8, "hom")):
                assert row[offset].split(",")[alt_index] == "|".join(map(str, counts[kind][0]))
                assert row[offset + 1].split(",")[alt_index] == str(counts[kind][1])
                assert row[offset + 2].split(",")[alt_index] == str(counts[kind][2])
    assert {"smaller", "larger", 0, 1} <= {where for _, where in seen}


# The default overlap of 0 skips the bedtools reciprocal-overlap filter, which rejects -f 0
@pytest.mark.parametrize("overlap", [{}, {"del_reciprocal_overlap": 0.5}], ids=["default", "half"])
def test_annotate_meds_labels_deletions_matching_catalog_elements_by_family(run_workflow, run_dir, vcf_inputs, source, overlap):
    deletions = [r for r in source[1] if r["info"]["allele_type"] == "del" and len(r["ref"]) > 50]
    # Only deletions with no other deletion nearby get a catalog element, so that none can match a neighbour's element
    isolated = [r for r in deletions if all(abs(o["pos"] - r["pos"]) > 1000 for o in deletions if o is not r)]
    designations = ["Alu", "L1", "Retroposon"]
    me_types = {"Alu": "ALU", "L1": "LINE", "Retroposon": "SVA"}
    catalog, expected = [], {}
    for i, r in enumerate(isolated):
        if i % 4 == 3:
            continue
        # The catalog element covers the deletion exactly; elements with an unrecognised designation are never reported
        designation = "Other" if i % 4 == 2 else designations[i % 3]
        catalog.append((CONTIG, r["pos"], r["stop"], f"element_{i}", r["ref"][1:], designation, f"sub{i}"))
        if designation != "Other":
            expected[r["id"]] = (me_types[designation], f"sub{i}")
    assert len(expected) > 5
    mei_catalog = write_table(run_dir / "meis.tsv", catalog)

    outputs = run_workflow(
        f"{ANNOTATION}/AnnotateMEDs.wdl",
        {
            **vcf_inputs,
            "contigs": [CONTIG],
            "prefix": "meds",
            "records_per_shard": 20,
            **overlap,
            "mei_catalog": str(mei_catalog),
        },
    )

    got = {row[4]: row for row in read_table(outputs["AnnotateMEDs.annotations_tsv_meds"])}
    assert {k: (v[5], v[6]) for k, v in got.items()} == expected
    by_id = {r["id"]: r for r in isolated}
    for rid, row in got.items():
        assert row[:4] == [CONTIG, str(by_id[rid]["pos"]), by_id[rid]["ref"], by_id[rid]["alts"][0]]


def test_annotate_meis_keeps_line_calls_and_only_corroborated_alu_and_sva_calls(run_workflow, run_dir, utils_docker, source):
    inserts = [r for r in source[1] if r["info"]["allele_type"] == "ins" and len(r["alts"]) == 1][:30]
    key = lambda r: [r["contig"], r["pos"], r["ref"], r["alts"][0], r["id"]]  # noqa: E731
    types = ["LINE", "ALU", "SVA"]
    l1meaid = [key(r) + [types[i % 3], "L1MEAID"] for i, r in enumerate(inserts)]
    palmer = [key(r) + ["ALU", "PALMER"] for i, r in enumerate(inserts) if i % 3 == 1 and i % 2 == 0]
    # SVAN tables hold named columns after the five key columns; FAM_N is the second of these
    svan = [key(r) + ["TSD", "sva" if i % 3 == 2 and i % 2 == 0 else "alu", "+"] for i, r in enumerate(inserts)]
    header = run_dir / "svan_header.txt"
    header.write_text("TSD\nFAM_N\nSTRAND\n")

    outputs = run_workflow(
        f"{ANNOTATION}/AnnotateMEIs.wdl",
        {
            "annotations_tsv_l1meaid": str(write_table(run_dir / "l1meaid.tsv", l1meaid)),
            "annotations_tsv_palmer": str(write_table(run_dir / "palmer.tsv", palmer)),
            "annotations_tsv_svan": str(write_table(run_dir / "svan.tsv", svan)),
            "annotations_header_svan": str(header),
            "contigs": [CONTIG],
            "prefix": "meis",
            "utils_docker": utils_docker,
        },
    )

    palmer_keys = {tuple(map(str, row[:6])) for row in palmer}
    svan_keys = {tuple(map(str, row[:5])) + (row[6].upper(),) for row in svan}
    expected = [
        [str(v) for v in row] for row in l1meaid if row[5] == "LINE" or tuple(map(str, row[:6])) in palmer_keys | svan_keys
    ]
    assert len(expected) > 10 and len(expected) < len(l1meaid)
    assert read_table(outputs["AnnotateMEIs.annotations_tsv_meis"]) == expected
