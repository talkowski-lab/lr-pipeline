#!/usr/bin/env python3
"""Prepare and interpret Paraphase 4.0.0 outputs for Mendelian candidates.

Standard-library only. WDL supplies offline VEP JSONL and an allele-normalized
ClinVar-annotated VCF between the prepare and report stages. No live APIs.
"""

import argparse
import csv
import gzip
import hashlib
import json
import re
import shutil
import subprocess
import sys
import tarfile
import tempfile
from pathlib import Path, PurePosixPath
from urllib.parse import unquote

VERSION = "0.1.1"
INTERPRETERS = {
    "smn1": ("SMN1", "smn1", {"autosomal_recessive"}),
    "pms2": ("PMS2", "pms2", {"autosomal_recessive", "autosomal_dominant"}),
    "rccx": ("CYP21A2", "rccx", {"autosomal_recessive"}),
    "f8": ("F8", "f8", {"x_linked"}),
}


def read_json(path):
    with open(path) as handle:
        return json.load(handle)


def write_json(path, data):
    Path(path).write_text(json.dumps(data, indent=2, sort_keys=True) + "\n")


def digest(path):
    sha = hashlib.sha256()
    with open(path, "rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            sha.update(chunk)
    return sha.hexdigest()


def load_rules(path):
    rules = read_json(path)
    if (rules["schema_version"], rules["paraphase_version"], rules["genome_build"]) != (
        1,
        "4.0.0",
        "GRCh38",
    ):
        raise ValueError("Only rules schema 1 / Paraphase 4.0.0 / GRCh38 supported")
    if rules["report_carriers"] is not False:
        raise ValueError("Carrier reporting is not supported")
    if rules["min_depth"] < 1 or not 0 < rules["min_alt_fraction"] <= 1:
        raise ValueError("Invalid quality thresholds")
    ids = set()
    for rule in rules["associations"]:
        if rule["id"] in ids:
            raise ValueError("Duplicate association id")
        ids.add(rule["id"])
        gene, region, inheritance = INTERPRETERS[rule["interpreter"]]
        if (rule["gene"], rule["region"]) != (gene, region):
            raise ValueError("Interpreter does not support this gene/region")
        if rule["inheritance"] not in inheritance:
            raise ValueError("Unsupported inheritance for interpreter")
        if rule["mechanism"] != "loss_of_function" or not rule["sources"]:
            raise ValueError("Only sourced loss-of-function associations supported")
        if rule["interpreter"] != "f8":
            if not rule["transcripts"] or len(rule["transcript_interval"]) != 3:
                raise ValueError("Transcript and GRCh38 interval required")
    return rules


def vcf_rows(path):
    opener = gzip.open if str(path).endswith((".gz", ".bgz")) else open
    columns = None
    with opener(path, "rt") as handle:
        for line in handle:
            if line.startswith("##"):
                continue
            if line.startswith("#CHROM"):
                columns = line.rstrip().split("\t")
                if len(columns) != len(set(columns)):
                    raise ValueError("Duplicate VCF columns")
                yield columns, None
            elif line.strip():
                fields = line.rstrip("\n").split("\t")
                if columns is None or len(fields) != len(columns):
                    raise ValueError(f"Malformed VCF: {path}")
                yield columns, fields
    if columns is None:
        raise ValueError(f"VCF header missing: {path}")


def info_fields(value):
    return dict(part.split("=", 1) for part in value.split(";") if "=" in part)


def base_hap(name):
    return name[:-4] if name.endswith("_cp2") else name


def variant_id(chrom, pos, ref, alt):
    key = f"{chrom}:{pos}:{ref}:{alt}"
    return "pp_" + hashlib.sha256(key.encode()).hexdigest()[:24]


def unpack_vcfs(archive_path, sample_id, regions):
    """Read the runner's <sample>_paraphase_vcfs/<sample>_<region>.vcf layout."""
    if (
        not sample_id
        or sample_id in {".", ".."}
        or "/" in sample_id
        or "\\" in sample_id
    ):
        raise ValueError("Archive sample_id must be a single filename component")
    directory = sample_id + "_paraphase_vcfs"
    expected = {
        f"{directory}/{sample_id}_{region}{suffix}": region
        for region in regions
        for suffix in (".vcf", ".vcf.gz", ".vcf.bgz")
    }
    destination = Path(tempfile.mkdtemp(prefix="paraphase_vcfs_", dir="."))
    vcfs, members = {}, {}
    with tarfile.open(archive_path, "r:gz") as archive:
        for member in archive:
            path = PurePosixPath(member.name)
            if path.is_absolute() or ".." in path.parts:
                raise ValueError(f"Unsafe VCF archive path: {member.name}")
            if member.isdir() and str(path) in {".", directory}:
                continue
            if not member.isfile():
                raise ValueError(
                    f"VCF archive must contain regular files: {member.name}"
                )
            region = expected.get(str(path))
            if region is None:
                raise ValueError(
                    "VCF archive member does not match sample/JSON regions: "
                    f"{member.name}"
                )
            if region in vcfs:
                raise ValueError(f"Duplicate region in VCF archive: {region}")
            target = destination / path.name
            with archive.extractfile(member) as source, target.open("xb") as out:
                shutil.copyfileobj(source, out)
            vcfs[region] = str(target)
            members[region] = member.name
    if not vcfs:
        raise ValueError("Supplied VCF archive contains no region VCFs")
    return vcfs, members


def prepare(request):
    request = dict(request)
    rules = load_rules(request["rules"])
    request.setdefault("region_vcfs", [])
    if isinstance(request["region_vcfs"], list):
        names = request.get("region_names", [])
        if len(names) != len(request["region_vcfs"]) or len(set(names)) != len(names):
            raise ValueError("region_names must be unique and match region_vcfs length")
        request["region_vcfs"] = dict(zip(names, request["region_vcfs"]))
    for field in ("paraphase_version", "genome_build"):
        if request[field] != rules[field]:
            raise ValueError(f"Unsupported {field}: {request[field]}")
    if not re.fullmatch(r"[A-Za-z0-9_.-]+", request["prefix"]):
        raise ValueError("Prefix must contain only letters, digits, _, . and -")
    raw = read_json(request["paraphase_json"])
    if not isinstance(raw, dict) or not raw:
        raise ValueError("Expected a nonempty region-keyed Paraphase JSON object")
    archive_members = {}
    if request.get("paraphase_vcfs"):
        if request["region_vcfs"]:
            raise ValueError("Supply a VCF tarball or explicit region VCFs, not both")
        request["region_vcfs"], archive_members = unpack_vcfs(
            request["paraphase_vcfs"], request["sample_id"], raw
        )
    payload = {"sample_id": request["sample_id"], "regions": {}, "sites": {}}
    for region, call in raw.items():
        if not isinstance(call, dict):
            raise ValueError(f"Invalid region object: {region}")
        required = {
            "region_name",
            "failed_for_coverage",
            "final_haplotypes",
            "two_copy_haplotypes",
            "region_specific_info",
            "total_cn",
        }
        if not required <= call.keys() or call["region_name"] != region:
            raise ValueError(f"Not a v4.0.0 region object: {region}")
        if not isinstance(call["failed_for_coverage"], bool):
            raise ValueError("failed_for_coverage must be boolean")
        info = call["region_specific_info"]
        if not call["failed_for_coverage"]:
            region_required = {
                "smn1": {"smn1_cn", "smn1_haplotypes"},
                "pms2": {"gene_cn"},
                "rccx": {"annotated_alleles", "phasing_success"},
                "f8": {"sv_called"},
            }.get(region, set())
            if not region_required <= info.keys():
                raise ValueError(f"Missing v4.0.0 region-specific fields: {region}")
        for cn in [call["total_cn"]] + [
            info.get(k) for k in ("gene_cn", "smn1_cn", "smn2_cn")
        ]:
            if cn is not None and (type(cn) is not int or cn < 0):
                raise ValueError("Copy number must be a nonnegative integer or null")
        names = list(call["final_haplotypes"].values())
        if len(set(names)) != len(names):
            raise ValueError(f"Duplicate haplotype names in {region}")
        duplicates = call["two_copy_haplotypes"]
        if not set(duplicates) <= set(names) or len(set(duplicates)) != len(duplicates):
            raise ValueError(f"Invalid two_copy_haplotypes in {region}")
        keep = {
            k: v
            for k, v in call.items()
            if k
            not in {
                "read_details",
                "unique_supporting_reads",
                "nonunique_supporting_reads",
                "assembled_haplotypes",
                "sites_for_phasing",
            }
        }
        keep.update(
            {
                "vcf_present": False,
                "vcf_haplotypes": [],
                "calls": [],
                "no_calls": [],
                "vcf_alleles": [],
            }
        )
        payload["regions"][region] = keep
    for region, path in request["region_vcfs"].items():
        if region not in raw:
            raise ValueError(f"VCF region absent from JSON: {region}")
        call = payload["regions"][region]
        names = set(call["final_haplotypes"].values())
        allowed = names | {h + "_cp2" for h in call["two_copy_haplotypes"]}
        for header, fields in vcf_rows(path):
            if fields is None:
                if set(header[9:]) - allowed:
                    raise ValueError(f"VCF haplotypes absent from JSON: {region}")
                call["vcf_present"] = True
                call["vcf_haplotypes"] = header[9:]
                continue
            info = info_fields(fields[7])
            if "ALLELE" in info:
                alleles = [x.split("+") for x in info["ALLELE"].split(",")]
                previous = call["vcf_alleles"]
                if previous and previous != alleles:
                    raise ValueError(f"Inconsistent VCF ALLELE values: {region}")
                call["vcf_alleles"] = alleles
            formats = fields[8].split(":") if len(fields) > 8 else []
            alts = fields[4].split(",")
            for hap, value in zip(header[9:], fields[9:]):
                data = dict(zip(formats, value.split(":")))
                gt = data.get("GT", ".")
                if gt == ".":
                    call["no_calls"].append(hap)
                    continue
                if not gt.isdigit() or int(gt) > len(alts):
                    raise ValueError(f"Expected haploid GT in {region}: {gt}")
                if int(gt) == 0:
                    continue
                alt = alts[int(gt) - 1]
                key = variant_id(fields[0], fields[1], fields[3], alt)
                site = {
                    "id": key,
                    "chrom": fields[0],
                    "pos": int(fields[1]),
                    "ref": fields[3],
                    "alt": alt,
                }
                site["small_variant"] = bool(
                    re.fullmatch("[ACGTN]+", alt)
                    and re.fullmatch("[ACGTN]+", fields[3])
                )
                payload["sites"][key] = site
                dp = int(data["DP"]) if data.get("DP", ".").isdigit() else None
                ad = data.get("AD", ".").split(",")
                support = (
                    int(ad[int(gt)])
                    if len(ad) > int(gt) and ad[int(gt)].isdigit()
                    else None
                )
                call["calls"].append(
                    {
                        "id": key,
                        "haplotype": hap,
                        "filter": fields[6],
                        "dp": dp,
                        "alt_depth": support,
                    }
                )
    contigs = {}
    with open(request["ref_fai"]) as handle:
        for line in handle:
            name, length, *_ = line.split()
            contigs[name] = int(length)
    prefix = request["prefix"]
    with open(prefix + ".sites.vcf", "w") as handle:
        handle.write("##fileformat=VCFv4.2\n")
        for chrom, length in contigs.items():
            handle.write(f"##contig=<ID={chrom},length={length}>\n")
        handle.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
        for key, site in sorted(payload["sites"].items()):
            if site["chrom"] not in contigs:
                raise ValueError(f"Contig absent from reference: {site['chrom']}")
            if site["small_variant"]:
                handle.write(
                    "\t".join(
                        map(
                            str,
                            [
                                site["chrom"],
                                site["pos"],
                                key,
                                site["ref"],
                                site["alt"],
                                ".",
                                "PASS",
                                ".",
                            ],
                        )
                    )
                    + "\n"
                )
    files = {
        "paraphase_json": request["paraphase_json"],
        "rules": request["rules"],
        "ref_fai": request["ref_fai"],
        **request["region_vcfs"],
    }
    if request.get("paraphase_vcfs"):
        files["paraphase_vcfs"] = request["paraphase_vcfs"]
    payload["provenance"] = {
        "parser_version": VERSION,
        "sample_id_source": "user_supplied_manifest",
        "parser_sha256": digest(__file__),
        "python_version": sys.version,
        "container": request.get("container"),
        "rules_version": rules["rules_version"],
        "rules_sha256": digest(request["rules"]),
        "paraphase_version": request["paraphase_version"],
        "genome_build": request["genome_build"],
        "gene1only": request["gene1only"],
        "targeted": request["targeted"],
        "input_sha256": {k: digest(v) for k, v in files.items()},
        "vcf_archive_members": archive_members,
    }
    write_json(prefix + ".prepared.json", payload)


def annotations(vcf, vep_json, sites):
    result = {}
    for _, row in vcf_rows(vcf):
        if row is not None:
            if row[2] not in sites or row[2] in result:
                raise ValueError("Unexpected or duplicated annotation variant ID")
            result[row[2]] = {
                "normalized_variant": ":".join([row[0], row[1], row[3], row[4]]),
                "clinvar": info_fields(row[7]),
            }
    seen = set()
    with open(vep_json) as handle:
        for line in handle:
            obj = json.loads(line)
            key = obj["input"].split("\t")[2]
            if key not in result or key in seen:
                raise ValueError("Unexpected or duplicated VEP variant ID")
            seen.add(key)
            result[key]["transcripts"] = obj.get("transcript_consequences", [])
    expected = {k for k, v in sites.items() if v["small_variant"]}
    if set(result) != expected or seen != expected:
        raise ValueError("Annotation incomplete: expected every small-variant ID")
    return result


def annotate(request):
    """Run local annotation tools; a pre-normalized GRCh38 ClinVar is required."""

    def run(args):
        subprocess.run(args, check=True)

    prefix = request["prefix"]
    Path("reference.fa").symlink_to(Path(request["ref_fa"]).resolve())
    Path("reference.fa.fai").symlink_to(Path(request["ref_fai"]).resolve())
    # Keep indexes beside the localized VCF without modifying input files.
    Path("clinvar.vcf.gz").symlink_to(Path(request["clinvar_vcf"]).resolve())
    Path("clinvar.vcf.gz.tbi").symlink_to(Path(request["clinvar_vcf_idx"]).resolve())
    cv_header = subprocess.check_output(
        ["bcftools", "view", "-h", "clinvar.vcf.gz"], text=True
    )
    cv_contigs = set(re.findall(r"##contig=<ID=([^,>]+)", cv_header))
    input_contigs = {r[0] for _, r in vcf_rows(request["sites_vcf"]) if r}
    if not input_contigs <= cv_contigs:
        raise ValueError("ClinVar contigs do not cover the input contig names")
    run(
        [
            "bcftools",
            "norm",
            "-f",
            "reference.fa",
            "--check-ref",
            "e",
            "-Ov",
            "-o",
            "normalized.vcf",
            request["sites_vcf"],
        ]
    )
    run(["bcftools", "sort", "-Oz", "-o", "sorted.vcf.gz", "normalized.vcf"])
    run(["bcftools", "index", "-t", "sorted.vcf.gz"])
    # ID is deliberately retained: it joins normalized alleles to original calls.
    annotated = prefix + ".annotated.vcf.gz"
    run(
        [
            "bcftools",
            "annotate",
            "-a",
            "clinvar.vcf.gz",
            "-c",
            "INFO/CLNSIG,INFO/CLNREVSTAT,INFO/CLNDN,INFO/CLNDISDB,INFO/ALLELEID",
            "--pair-logic",
            "exact",
            "-Oz",
            "-o",
            annotated,
            "sorted.vcf.gz",
        ]
    )
    run(["bcftools", "index", "-t", annotated])
    vep_output = prefix + ".vep.jsonl"
    has_sites = any(row is not None for _, row in vcf_rows(annotated))
    if has_sites:
        Path("cache").mkdir()
        with tarfile.open(request["ref_vep_cache"], "r:gz") as archive:
            # Cache archives need regular files/directories only. Refuse links
            # and traversal rather than trusting paths from an external archive.
            root = Path("cache").resolve()
            for member in archive.getmembers():
                dest = (root / member.name).resolve()
                if (root not in dest.parents and dest != root) or not (
                    member.isfile() or member.isdir()
                ):
                    raise ValueError("Unsafe path or link in VEP cache archive")
            if hasattr(tarfile, "data_filter"):
                archive.extractall(root, filter="data")
            else:
                archive.extractall(root)
        run(
            [
                request["vep_executable"],
                "--input_file",
                annotated,
                "--format",
                "vcf",
                "--output_file",
                vep_output,
                "--json",
                "--cache",
                "--offline",
                "--dir_cache",
                "cache",
                "--cache_version",
                str(request["vep_cache_version"]),
                "--assembly",
                "GRCh38",
                "--fasta",
                "reference.fa",
                "--symbol",
                "--hgvs",
                "--no_stats",
                "--force_overwrite",
            ]
        )
    else:
        Path(vep_output).write_text("")
    write_json(
        prefix + ".resources.json",
        {
            "clinvar_release": request["clinvar_release"],
            "clinvar_requires": (
                "biallelic, left-normalized GRCh38; matching contig names; tabix index"
            ),
            "vep_cache_version": request["vep_cache_version"],
            "container": request.get("container"),
            "vep_executable": request["vep_executable"],
            "bcftools_version": subprocess.check_output(
                ["bcftools", "--version"], text=True
            ),
            "resource_sha256": {
                k: digest(request[k])
                for k in ("ref_fa", "ref_fai", "clinvar_vcf", "ref_vep_cache")
            },
        },
    )


def gene_haplotypes(call, rule):
    names = list(call["final_haplotypes"].values())
    if rule["interpreter"] == "smn1":
        return set(call["region_specific_info"].get("smn1_haplotypes", {}).values())
    if rule["interpreter"] == "pms2":
        return {h for h in names if re.fullmatch(r"pms2_pms2hap\d+", h)}
    return set(names)


def boundary_state(call, hap, rule):
    detail = call.get("haplotype_details", {}).get(base_hap(hap), {})
    boundary = detail.get("boundary")
    if not boundary or len(boundary) != 2:
        return "unknown"
    if detail.get("is_truncated"):
        return "truncated"
    start, end = rule["transcript_interval"][1:]
    if boundary[1] < start or boundary[0] > end:
        return "no_overlap"
    return "full_span" if boundary[0] <= start and boundary[1] >= end else "partial"


def evidence(annotation, rule, rules):
    tx = [
        t
        for t in annotation.get("transcripts", [])
        if t.get("gene_symbol") == rule["gene"]
        and t.get("transcript_id", "").split(".")[0] in rule["transcripts"]
    ]
    consequences = sorted({c for t in tx for c in t.get("consequence_terms", [])})
    cv = annotation.get("clinvar", {})
    significance = unquote(cv.get("CLNSIG", "")).lower()
    terms = set(re.split(r"[|/,]", significance))
    plp = bool(terms & {"pathogenic", "likely_pathogenic"})
    conflict = "conflict" in significance or (
        plp and bool(terms - {"pathogenic", "likely_pathogenic"})
    )
    condition = unquote(cv.get("CLNDN", "")).replace("_", " ").lower()
    condition_match = any(term in condition for term in rule["condition_terms"])
    reviewed = cv.get("CLNREVSTAT") in rules["clinvar_review_statuses"]
    kind = "none"
    if tx and plp and not conflict:
        kind = "clinvar_plp" if condition_match and reviewed else "clinvar_plp_review"
    elif tx and set(consequences) & set(rules["plof_consequences"]):
        kind = "predicted_lof_review"
    elif tx and conflict:
        kind = "clinvar_conflict_review"
    return {
        "evidence": kind,
        "consequences": consequences,
        "transcripts": tx,
        "clinvar": cv,
        "condition_match": condition_match,
        "clinvar_conflict": conflict,
    }


def allele_groups(call):
    groups = call["region_specific_info"].get("alleles_final") or call["vcf_alleles"]
    if not groups or len(groups) != 2 or not all(groups):
        return []
    names = set(call["final_haplotypes"].values())
    if any(base_hap(h) not in names for group in groups for h in group):
        raise ValueError("Allele assignment refers to an unknown haplotype")
    bases = [base_hap(h) for group in groups for h in group]
    for hap in set(bases):
        maximum = 2 if hap in call["two_copy_haplotypes"] else 1
        if bases.count(hap) > maximum:
            raise ValueError("Allele assignments exceed copy multiplicity")
    vcf_groups = call["vcf_alleles"]
    if vcf_groups:

        def canonical(gs):
            return sorted(sorted(base_hap(h) for h in g) for g in gs)

        if canonical(groups) != canonical(vcf_groups):
            raise ValueError("JSON and VCF allele assignments disagree")
    return groups


def sequence_result(call, rule, rules, annotated):
    selected = gene_haplotypes(call, rule)
    copies = sorted(h for h in call["vcf_haplotypes"] if base_hap(h) in selected)
    expected = selected | {
        h + "_cp2" for h in call["two_copy_haplotypes"] if h in selected
    }
    cn_key = "smn1_cn" if rule["interpreter"] == "smn1" else "gene_cn"
    cn = call["region_specific_info"].get(cn_key)
    warnings = []
    if set(copies) != expected or (cn is not None and cn != len(expected)):
        warnings.append("copy_number_or_vcf_mismatch")
    if cn is None:
        warnings.append("copy_number_unresolved")
    rows = []
    for item in call["calls"]:
        ann = annotated.get(item["id"], {})
        ev = evidence(ann, rule, rules)
        if ev["evidence"] == "none":
            continue
        quality = (
            item["filter"] == "PASS"
            and item["dp"] is not None
            and item["dp"] >= rules["min_depth"]
            and item["alt_depth"] is not None
            and item["alt_depth"] / item["dp"] >= rules["min_alt_fraction"]
        )
        rows.append(
            {
                **item,
                **ev,
                "assigned_gene": base_hap(item["haplotype"]) in selected,
                "quality_pass": quality,
                "normalized_variant": ann.get("normalized_variant"),
            }
        )
    assigned = [r for r in rows if r["assigned_gene"]]
    affected = {r["haplotype"] for r in assigned}
    # Flag overlapping/nearby changes instead of declaring independent effects.
    for hap in affected:
        positions = sorted(
            int(annotated[r["id"]]["normalized_variant"].split(":")[1])
            for r in call["calls"]
            if r["haplotype"] == hap and r["id"] in annotated
        )
        if any(b - a <= 10 for a, b in zip(positions, positions[1:])):
            warnings.append("nearby_variants_require_combined_effect_review")
    if any(not r["quality_pass"] for r in assigned):
        warnings.append("low_quality_variant")
    if any(r["evidence"] != "clinvar_plp" for r in assigned):
        warnings.append("variant_evidence_requires_review")
    status, reason = "not_prioritized", "no_qualifying_configuration"
    if rule["interpreter"] == "smn1" and cn == 0 and not expected:
        status, reason = "candidate", "zero_smn1_copies"
    elif rule["interpreter"] == "pms2" and cn == 0 and not expected:
        status, reason = (
            "unresolved_candidate",
            "zero_pms2_copies_requires_deletion_review",
        )
    elif rule["inheritance"] == "autosomal_dominant" and assigned:
        status, reason = "candidate", "dominant_qualifying_variant"
    elif rule["inheritance"] == "autosomal_dominant" and cn == 1:
        status, reason = (
            "unresolved_candidate",
            "reduced_pms2_copies_requires_deletion_review",
        )
    elif assigned:
        groups = allele_groups(call)
        affected_bases = {base_hap(h) for h in affected}
        # A chromosome is affected only if ALL assigned gene copies on it qualify.
        # Repeated base names may represent identical copies on both chromosomes.
        group_genes = [{base_hap(h) for h in g} & selected for g in groups]
        fully_affected = (
            len(group_genes) == 2
            and all(group_genes)
            and all(g <= affected_bases for g in group_genes)
            and expected <= affected
        )
        if fully_affected:
            status, reason = "candidate", "qualifying_copies_on_both_chromosomes"
        elif rule["interpreter"] == "smn1" and cn == 1 and expected <= affected:
            status, reason = "candidate", "one_smn1_copy_with_qualifying_variant"
        elif rule["interpreter"] == "pms2" and cn == 1 and expected <= affected:
            status, reason = "unresolved_candidate", "possible_copy_loss_plus_variant"
        elif len(affected) >= 2:
            same_chromosome = (
                len(groups) == 2
                and any(affected_bases <= g for g in group_genes)
                and not any(
                    affected_bases & g for g in group_genes if not affected_bases <= g
                )
            )
            status = "not_prioritized" if same_chromosome else "unresolved_candidate"
            reason = (
                "qualifying_copies_in_cis" if same_chromosome else "phase_unresolved"
            )
        else:
            status, reason = "not_prioritized", "single_recessive_allele"
            other = set(copies) - affected
            if (
                warnings
                or not other
                or any(
                    boundary_state(call, h, rule) != "full_span"
                    or h in call["no_calls"]
                    for h in other
                )
            ):
                status, reason = (
                    "unresolved_candidate",
                    "other_allele_incompletely_assessed",
                )
    if any("unknownhap" in r["haplotype"] for r in rows):
        warnings.append("gene_assignment_unresolved")
        if status == "not_prioritized":
            status, reason = "unresolved_candidate", "gene_assignment_unresolved"
    if status == "candidate" and warnings:
        status = "unresolved_candidate"
    return status, reason, sorted(set(warnings)), rows, cn, copies


def structural_result(call, rule, x_copy_number):
    info = call["region_specific_info"]
    status, reason = "not_prioritized", "no_supported_structural_configuration"
    if rule["interpreter"] == "f8":
        events = info.get("sv_called") or {}
        qualifying = {h: e for h, e in events.items() if e in {"inversion", "deletion"}}
        if qualifying:
            if x_copy_number == 1:
                status, reason = "candidate", "f8_event_in_one_x_sample"
            elif x_copy_number == 2 and len(qualifying) == 1:
                status, reason = "not_prioritized", "single_x_linked_event"
            else:
                status, reason = (
                    "unresolved_candidate",
                    "f8_chromosome_dosage_or_phase_unresolved",
                )
        return status, reason, events
    labels = info.get("annotated_alleles") or []
    # These labels only summarize known differentiating sites. Never call a
    # whole chromosome damaged if an additional WT-like CYP21A2 copy is present.
    intact_labels = {
        "WT",
        "pseudogene_deletion",
        "pseudogene_duplication",
        "gene_duplication",
    }
    damaged = [
        label
        for label in labels
        if label
        and label not in intact_labels
        and not label.startswith("duplication_WT_plus_")
    ]
    if len(damaged) >= 2:
        status, reason = (
            "unresolved_candidate",
            "rccx_two_flagged_alleles_require_review",
        )
    elif damaged and (
        len(labels) != 2 or None in labels or not info.get("phasing_success")
    ):
        status, reason = "unresolved_candidate", "rccx_other_allele_unresolved"
    elif damaged:
        reason = "single_recessive_allele"
    return (
        status,
        reason,
        {
            "annotated_alleles": labels,
            "phasing_success": info.get("phasing_success"),
            "hap_variants": info.get("hap_variants"),
        },
    )


def table(path, columns, rows):
    with open(path, "w") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=columns,
            delimiter="\t",
            extrasaction="ignore",
            lineterminator="\n",
        )
        writer.writeheader()
        for row in rows:
            writer.writerow(
                {
                    k: (
                        json.dumps(v, sort_keys=True)
                        if isinstance(v, (dict, list))
                        else "" if v is None else v
                    )
                    for k, v in row.items()
                }
            )


def report(request):
    rules = load_rules(request["rules"])
    if request.get("x_copy_number") not in (None, 1, 2):
        raise ValueError("x_copy_number must be 1, 2, or omitted")
    payload = read_json(request["prepared"])
    if digest(request["rules"]) != payload["provenance"]["rules_sha256"]:
        raise ValueError("Rules changed between prepare and report")
    annotated = annotations(
        request["annotated_vcf"], request["vep_json"], payload["sites"]
    )
    summaries, variant_rows, hap_rows, qc = [], [], [], []
    regions = payload["regions"]
    for name, call in regions.items():
        if name not in {r["region"] for r in rules["associations"]}:
            qc.append({"region": name, "code": "no_gene_disease_rule"})
        if name != "f8" and any(
            not payload["sites"][item["id"]]["small_variant"] for item in call["calls"]
        ):
            qc.append({"region": name, "code": "symbolic_vcf_alleles_not_interpreted"})
        for hap in call["final_haplotypes"].values():
            hap_rows.append(
                {
                    "sample_id": payload["sample_id"],
                    "region": name,
                    "haplotype": hap,
                    "copy_count": 2 if hap in call["two_copy_haplotypes"] else 1,
                    "details": call.get("haplotype_details", {}).get(hap),
                    "region_specific_info": call["region_specific_info"],
                }
            )
    for rule in rules["associations"]:
        row = {
            k: rule[k]
            for k in ("gene", "disease", "inheritance", "mechanism", "region")
        }
        row.update(
            {
                "sample_id": payload["sample_id"],
                "association_id": rule["id"],
                "status": "not_assessed",
                "warnings": [],
                "evidence": [],
                "gene_copy_number": None,
                "affected_copies": [],
                "assessment": "not_a_negative_disease_test",
            }
        )
        call = regions.get(rule["region"])
        if call is None:
            row["reason"] = "region_missing"
        elif call["failed_for_coverage"]:
            row["reason"] = "failed_for_coverage"
        elif rule["interpreter"] in {"rccx", "f8"}:
            row["status"], row["reason"], events = structural_result(
                call, rule, request.get("x_copy_number")
            )
            row["evidence"] = events
            row["assessment"] = "selected_structural_or_marker_events_only"
            # Generic RCCX annotation cannot assign a functional copy. Preserve
            # evidence, but do not promote routine pseudogene differences.
            if rule["interpreter"] == "rccx":
                _, _, _, variants, _, _ = sequence_result(call, rule, rules, annotated)
                if variants:
                    row["warnings"].append(
                        "rccx_small_variants_not_genotype_interpreted"
                    )
                    for variant in variants:
                        variant["assigned_gene"] = None
                        variant_rows.append({**row, **variant})
        elif not call["vcf_present"]:
            row["reason"] = "vcf_missing"
        else:
            status, reason, warnings, variants, cn, copies = sequence_result(
                call, rule, rules, annotated
            )
            row.update(
                {
                    "status": status,
                    "reason": reason,
                    "warnings": warnings,
                    "gene_copy_number": cn,
                    "affected_copies": sorted(
                        {r["haplotype"] for r in variants if r["assigned_gene"]}
                    ),
                    "evidence": variants,
                    "assessment": {h: boundary_state(call, h, rule) for h in copies},
                }
            )
            for variant in variants:
                variant_rows.append({**row, **variant})
        if row["status"] == "not_assessed":
            qc.append({"region": rule["region"], "code": row["reason"]})
        for warning in row["warnings"]:
            qc.append({"region": rule["region"], "code": warning})
        summaries.append(row)
    prefix = request["prefix"]
    columns = [
        "sample_id",
        "association_id",
        "gene",
        "disease",
        "inheritance",
        "mechanism",
        "region",
        "status",
        "reason",
        "gene_copy_number",
        "affected_copies",
        "assessment",
        "warnings",
        "evidence",
    ]
    prioritized = [
        r for r in summaries if r["status"] in {"candidate", "unresolved_candidate"}
    ]
    table(prefix + ".findings.tsv", columns, prioritized)
    table(prefix + ".gene_summary.tsv", columns, summaries)
    table(
        prefix + ".variants.tsv",
        columns[:7]
        + [
            "id",
            "normalized_variant",
            "haplotype",
            "assigned_gene",
            "evidence",
            "quality_pass",
            "dp",
            "alt_depth",
            "filter",
            "consequences",
            "clinvar",
            "condition_match",
            "clinvar_conflict",
            "transcripts",
        ],
        variant_rows,
    )
    table(
        prefix + ".haplotypes.tsv",
        [
            "sample_id",
            "region",
            "haplotype",
            "copy_count",
            "details",
            "region_specific_info",
        ],
        hap_rows,
    )
    table(prefix + ".qc.tsv", ["region", "code"], qc)
    provenance = {
        **payload["provenance"],
        "annotation_resources": read_json(request["resources"]),
        "x_copy_number": request.get("x_copy_number"),
        "annotated_vcf_sha256": digest(request["annotated_vcf"]),
        "vep_json_sha256": digest(request["vep_json"]),
    }
    write_json(prefix + ".provenance.json", provenance)
    write_json(
        prefix + ".report.json",
        {"associations": summaries, "qc": qc, "provenance": provenance},
    )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "stage", choices=["prepare", "annotate", "report", "validate-rules"]
    )
    parser.add_argument(
        "--request",
        required=True,
        help="JSON request file (rules file for validate-rules)",
    )
    args = parser.parse_args()
    if args.stage == "validate-rules":
        load_rules(args.request)
    else:
        {"prepare": prepare, "annotate": annotate, "report": report}[args.stage](
            read_json(args.request)
        )


if __name__ == "__main__":
    main()
