#!/usr/bin/env python3
"""Generate the example input JSON of every directly-run workflow into inputs/<dir>/<Workflow>.json.

Each file holds the required inputs, every `*_docker` input and a few common optional ones, filled with
standard placeholder values chosen by input name first and then by type. Reference and image values come
from docs/references.md and docs/dockers.md. Run after changing any workflow's `input` block; `--check`
exits 1 when a committed file differs from the regenerated one.
"""

import argparse
import json
import re
import sys

from wdl_meta import DOCKER_INPUT_RE, REPO_ROOT, parse_workflow

INPUTS_DIR = REPO_ROOT / "inputs"
COHORT = "gs://<cohort-bucket>/"
REFERENCES = "gs://<references-bucket>/"
REGISTRY = "us-central1-docker.pkg.dev/talkowski-sv-gnomad/kj-dockers/"

# Optional inputs worth showing even though they have a default
COMMON_OPTIONAL = ["records_per_shard", "contigs", "ped"]

# Placeholders chosen by exact input name; the first match wins before any type rule
NAMED = {
    "prefix": "cohort",
    "contig": "chr22",
    "contigs": [*(f"chr{i}" for i in range(1, 23)), "chrX", "chrY"],
    "records_per_shard": 10000,
    "samples_per_shard": 100,
    "bin_size": 1000,
    "sample_id": "sample1",
    "sample": "sample1",
    "sample_ids": ["sample1", "sample2"],
    "sex": "female",
    "sexes": ["female", "male"],
    "sample_sexes": ["female", "male"],
    "remote_indir": COHORT + "indir",
    "remote_outdir": COHORT + "outdir",
    "vcf": COHORT + "cohort.vcf.gz",
    "cohort_vcf": COHORT + "cohort.vcf.gz",
    "ped": COHORT + "cohort.ped",
    "ped_file": COHORT + "cohort.ped",
    "pedigree": COHORT + "cohort.ped",
    "bam": COHORT + "sample.bam",
    "bai": COHORT + "sample.bam.bai",
}

# Reference-like inputs absent from docs/references.md; they live in the references bucket
EXTRA_REFERENCES = {
    "autosomes_bed", "standard_chromosomes_bed", "exclude_bed", "expected_cn_female", "expected_cn_male",
    "ref_agp", "ref_map_file", "ref_bundle_json_file", "training_resource_bed", "training_resource_vcf",
    "mei_bwa_idx", "ref_bwa_idx", "intervals",
}  # fmt: skip

# Extension of a File input is taken from its name suffix, e.g. `annotations_tsv` -> annotations.tsv
SUFFIX_EXT = {"_tsv": ".tsv", "_bed": ".bed", "_json": ".json", "_csv": ".csv", "_gtf": ".gtf", "_script": ".py", "_agp": ".agp"}

# Array inputs whose placeholder list has two elements; the rest get one
TWO_ELEMENTS = re.compile(r"(vcfs|vcf_idxs|bams|bais)$")


def read_table(path, section_start, section_end=None):
    """Rows of the first-column-backticked markdown tables between two headings."""
    text = (REPO_ROOT / path).read_text()
    start = text.index(section_start)
    end = text.index(section_end, start) if section_end else len(text)
    return [line for line in text[start:end].splitlines() if re.match(r"^\| `\w+` \|", line)]


def load_references():
    """Input name -> file name of the copy documented in docs/references.md."""
    refs = {}
    for row in read_table("docs/references.md", "## Shared References", "## Released Callset Inputs"):
        url = re.search(r"\]\((gs://[^)]+)\)", row)
        if url:
            refs[re.match(r"^\| `(\w+)`", row).group(1)] = url.group(1).rstrip("/").rsplit("/", 1)[1]
    return refs


def load_dockers():
    """`*_docker` input name -> image URI documented in docs/dockers.md."""
    dockers = {}
    for row in read_table("docs/dockers.md", "## Repository"):
        cells = [c.strip().strip("`") for c in row.strip("|").split("|")]
        name, image = cells[0], cells[1]
        dockers[name] = image if "/" in image or ":" in image else f"{REGISTRY}{image}:latest"
    return dockers


REFS = load_references()
DOCKERS = load_dockers()


def singular(name):
    return name[:-2] if name.endswith("sexes") else name[:-1] if name.endswith("s") else name


def file_value(name):
    if name in NAMED:
        return NAMED[name]
    if name in REFS:
        return REFERENCES + REFS[name]
    if name.endswith("_idx") and name not in EXTRA_REFERENCES:
        return file_value(name[: -len("_idx")]) + ".tbi"
    if name.endswith("_fai") or name.endswith("_fai_gz"):
        return file_value(name.replace("_fai", "_fa")) + ".fai"
    if re.search(r"^vcf_\w$|_vcf$", name):
        return COHORT + re.sub(r"_vcf$", "", name) + ".vcf.gz"
    if name.endswith(("_bam", "_bai")):
        return COHORT + name[:-4] + (".bam" if name.endswith("_bam") else ".bam.bai")
    bucket = REFERENCES if name in EXTRA_REFERENCES else COHORT
    for suffix, ext in SUFFIX_EXT.items():
        if name.endswith(suffix):
            return bucket + name[: -len(suffix)] + ext
    return bucket + name


def scalar(name, type_name):
    if type_name == "Boolean":
        return False
    if type_name == "Float":
        return 0.0
    if type_name == "Int":
        return NAMED.get(name, 0)
    if type_name == "File" or re.search(r"(bam|bai|vcf|vcf_idx)$", name):
        return file_value(name)
    if DOCKER_INPUT_RE.match(name):
        if name not in DOCKERS:
            raise KeyError(f"{name} is not listed in docs/dockers.md")
        return DOCKERS[name]
    if name in NAMED:
        return NAMED[name]
    return REFERENCES + REFS[name] if name in REFS else name


def value(name, type_name, plural=True):
    type_name = type_name.rstrip("?+")
    if name in NAMED and type_name.startswith("Array"):
        return NAMED[name]
    if type_name.startswith("Array["):
        inner = type_name[len("Array[") : -1]
        element = value(singular(name) if plural else name, inner, plural=False)
        if TWO_ELEMENTS.search(name) and isinstance(element, str):
            stem, dot, rest = element.partition(".")
            return [f"{stem}_{i}{dot}{rest}" for i in (1, 2)]
        return [element]
    return scalar(name, type_name)


def selected(decl):
    return (not decl.optional and decl.default is None) or DOCKER_INPUT_RE.match(decl.name) or decl.name in COMMON_OPTIONAL


def dockstore_workflows():
    """Descriptor paths of every workflow registered in .dockstore.yml."""
    text = (REPO_ROOT / ".dockstore.yml").read_text()
    return [REPO_ROOT / p.lstrip("/") for p in re.findall(r"^\s*primaryDescriptorPath:\s*(\S+)", text, re.MULTILINE)]


def example_inputs():
    """Map of output path -> rendered JSON for every registered workflow."""
    files = {}
    for path in dockstore_workflows():
        workflow = parse_workflow(path)
        data = {f"{workflow.name}.{d.name}": value(d.name, d.type) for d in workflow.inputs if selected(d)}
        out = INPUTS_DIR / path.parent.relative_to(REPO_ROOT / "wdl") / f"{workflow.name}.json"
        files[out] = json.dumps(data, indent=4) + "\n"
    return files


def main():
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--check", action="store_true", help="exit 1 if committed files differ from regenerated output")
    args = parser.parse_args()

    expected = example_inputs()
    existing = {p: p.read_text() for p in INPUTS_DIR.rglob("*.json")}
    diffs = sorted(p for p in expected.keys() | existing.keys() if expected.get(p) != existing.get(p))
    if args.check:
        for path in diffs:
            print(f"out of date: {path.relative_to(REPO_ROOT)}")
        if diffs:
            print("Run: python .github/scripts/generate_example_inputs.py")
        return 1 if diffs else 0

    for path in diffs:
        if path in expected:
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_text(expected[path])
        else:
            path.unlink()
    print(f"Wrote {len(expected)} example inputs ({len(diffs)} changed)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
