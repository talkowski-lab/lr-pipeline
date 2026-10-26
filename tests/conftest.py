import gzip
import importlib.util
import json
import os
import re
import shutil
import subprocess
import sys
from pathlib import Path

import pysam
import pytest

REPO_ROOT = Path(__file__).resolve().parents[1]
RUNS_DIR = REPO_ROOT / "tests" / ".runs"
SUBSET_VCF = REPO_ROOT / "tests" / "data" / "hprc_hgsvc_v10.chr22.subset.vcf.gz"
FAI = REPO_ROOT / "tests" / "data" / "GRCh38.chr22.fa.fai"
CONTIG = "chr22"
CONTIG_LENGTH = 50818468
SAMPLES = ["HG00096", "HG00097", "HG00099", "HG00126", "HG00128"]
UTILS = "wdl/annotation_utils"
ANNOTATION = "wdl/annotation"
TOOLS = "wdl/tools"

# The checker scripts import one another by bare module name
sys.path.insert(0, str(REPO_ROOT / ".github" / "scripts"))


def load_script(rel_path):
    """Import a script under the repo root by path, without executing its __main__ block."""
    path = REPO_ROOT / rel_path
    spec = importlib.util.spec_from_file_location(path.stem, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


@pytest.fixture(scope="session")
def repo_root():
    return REPO_ROOT


@pytest.fixture(scope="session")
def subset_vcf():
    return SUBSET_VCF


@pytest.fixture(scope="session")
def utils_docker():
    image = os.environ.get("UTILS_DOCKER")
    if not image:
        pytest.skip("UTILS_DOCKER is not set")
    return image


class WorkflowError(Exception):
    pass


@pytest.fixture
def run_dir(request):
    """Return a clean per-test directory under tests/.runs/, which is git-ignored."""
    path = RUNS_DIR / re.sub(r"[^\w.-]", "_", request.node.name)
    shutil.rmtree(path, ignore_errors=True)
    path.mkdir(parents=True)
    return path


@pytest.fixture
def run_workflow(run_dir):
    """Return a function that runs a WDL under miniwdl and returns its outputs dict."""

    runs = []

    def run(wdl_rel_path, inputs):
        runs.append(None)
        inputs_path = run_dir / f"inputs{len(runs)}.json"
        inputs_path.write_text(json.dumps(inputs, indent=2))
        # A trailing "/." makes miniwdl run directly in this directory instead of a timestamped subdirectory
        command = [
            "miniwdl",
            "run",
            str(REPO_ROOT / wdl_rel_path),
            "-i",
            str(inputs_path),
            "--dir",
            f"{run_dir}/run{len(runs)}/.",
            "--no-cache",
        ]
        result = subprocess.run(command, capture_output=True, text=True, env=_miniwdl_env())
        if result.returncode != 0:
            raise WorkflowError(f"miniwdl exited {result.returncode}\n{result.stderr[-4000:]}")
        return json.loads(result.stdout)["outputs"]

    return run


def _miniwdl_env():
    env = dict(os.environ)
    # miniwdl runs podman through sudo and always pulls, so use a wrapper that is rootless and keeps locally built images
    if env.get("MINIWDL__SCHEDULER__CONTAINER_BACKEND") == "podman":
        env.setdefault("MINIWDL__PODMAN__EXE", json.dumps([str(REPO_ROOT / "tests" / "podman_local.sh")]))
    return env


def load_records(path):
    """Read a VCF into plain dicts so assertions do not depend on pysam record lifetimes."""
    with pysam.VariantFile(str(path)) as vcf:
        samples = list(vcf.header.samples)
        records = [
            {
                "contig": r.contig,
                "pos": r.pos,
                "stop": r.stop,
                "id": r.id,
                "ref": r.ref,
                "alts": r.alts,
                "filter": tuple(r.filter.keys()),
                "info": dict(r.info.items()),
                "gts": {s: r.samples[s]["GT"] for s in samples},
                "phased": {s: r.samples[s].phased for s in samples},
                "format": {s: {k: v for k, v in r.samples[s].items() if k != "GT"} for s in samples},
            }
            for r in vcf
        ]
    return samples, records


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


def write_table(path, rows, header=None):
    """Write rows (and an optional header row) as tab-separated text and return the path."""
    lines = ([header] if header else []) + list(rows)
    path.write_text("".join("\t".join(str(v) for v in row) + "\n" for row in lines))
    return path


def index_table(path, preset):
    """Bgzip a tab-separated file in place, tabix-index it with preset ("bed" or "vcf") and return the .gz path."""
    pysam.tabix_index(str(path), preset=preset, force=True)
    return path.with_name(path.name + ".gz")


def write_sites_vcf(path, rows, info=(), contig=CONTIG, length=CONTIG_LENGTH):
    """Write a sites-only VCF from (pos, id, ref, alt, info text) rows, bgzipped and tabix-indexed, and return the .gz path.

    info is a sequence of (ID, Number, Type) INFO header definitions.
    """
    lines = ["##fileformat=VCFv4.2", f"##contig=<ID={contig},length={length}>"]
    lines += [f'##INFO=<ID={tag},Number={number},Type={kind},Description="{tag}">' for tag, number, kind in info]
    lines.append("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO")
    lines += [f"{contig}\t{pos}\t{rid}\t{ref}\t{alt}\t.\t.\t{fields or '.'}" for pos, rid, ref, alt, fields in rows]
    path.write_text("\n".join(lines) + "\n")
    return index_table(path, "vcf")


def write_mosdepth_bed(path, runs, contig=CONTIG):
    """Write run-length (start, end, depth) mosdepth rows as a bgzipped, tabix-indexed BED and return the .gz path."""
    return index_table(write_table(path, [[contig, start, end, depth] for start, end, depth in runs]), "bed")


def read_table(path):
    """Read a possibly gzipped TSV into a list of rows, keeping header lines."""
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt") as handle:
        return [line.rstrip("\n").split("\t") for line in handle if line.strip()]


def read_dicts(path):
    """Read a possibly gzipped TSV with a header row into a list of dicts."""
    rows = read_table(path)
    header = [h.lstrip("#") for h in rows[0]]
    return [dict(zip(header, row)) for row in rows[1:]]


@pytest.fixture(scope="session")
def source():
    return load_records(SUBSET_VCF)


@pytest.fixture
def vcf_inputs(subset_vcf, utils_docker):
    return {"vcf": str(subset_vcf), "vcf_idx": f"{subset_vcf}.tbi", "utils_docker": utils_docker}


WINDOW_FASTA = REPO_ROOT / "tests" / "data" / "chr22.windows.fa"
WINDOW_WIDTH, WINDOW_OFFSET = 1000, 400


def window_records(records, count=60):
    """Evenly spaced records of at most 200 bp, one per reference window of chr22.windows.fa, in window order."""
    eligible = [r for r in records if len(r["ref"]) <= 200]
    return eligible[:: len(eligible) // count][:count]


def window_position(index):
    """Position in the window FASTA of the record held by the index-th window."""
    return index * WINDOW_WIDTH + WINDOW_OFFSET + 1


def write_window_fasta(reference, path, records):
    """Cut a WINDOW_WIDTH bp reference window around each record and join them into one chr22 sequence."""
    with pysam.FastaFile(str(reference)) as fasta:
        windows = [
            fasta.fetch(r["contig"], r["pos"] - 1 - WINDOW_OFFSET, r["pos"] - 1 - WINDOW_OFFSET + WINDOW_WIDTH).upper()
            for r in window_records(records)
        ]
    sequence = "".join(windows)
    path.write_text(f">{CONTIG}\n" + "".join(sequence[i : i + 60] + "\n" for i in range(0, len(sequence), 60)))


@pytest.fixture
def window_vcf(run_dir, subset_vcf, source):
    """The window_records() of the fixture VCF moved into the coordinates of chr22.windows.fa, as a bgzipped, indexed VCF."""
    order = {site_key(r): i for i, r in enumerate(window_records(source[1]))}
    plain = run_dir / "windows.vcf"
    with pysam.VariantFile(str(subset_vcf)) as vcf_in:
        header = pysam.VariantHeader()
        for line in str(vcf_in.header).splitlines():
            if not line.startswith("#CHROM"):
                header.add_line(line.replace(f"length={CONTIG_LENGTH}>", f"length={WINDOW_WIDTH * len(order)}>"))
        for sample in vcf_in.header.samples:
            header.add_sample(sample)
        records = []
        for record in vcf_in:
            if (record.contig, record.pos, record.ref, record.alts) in order:
                new = record.copy()
                new.translate(header)
                new.pos = window_position(order[(record.contig, record.pos, record.ref, record.alts)])
                records.append(new)
        with pysam.VariantFile(str(plain), "w", header=header) as vcf_out:
            for record in sorted(records, key=lambda r: r.pos):
                vcf_out.write(record)
    pysam.tabix_index(str(plain), preset="vcf", force=True)
    return Path(f"{plain}.gz")
