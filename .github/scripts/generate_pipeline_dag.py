#!/usr/bin/env python3
"""Generate the Run Order graph in docs/pipeline.md from the method config snapshot in docs/method_configs.json.

One config feeds another when the second reads a data table or workspace attribute the first writes, and each
config is coloured by the numbered pipeline.md section that links its workflow. The graph is checked before it
is written: every edge must go to a higher step number, no attribute may have two writers, every attribute read
must be written by a config or be a cohort input, every workflow must be linked from exactly one section, and
the configs must match the Run Order table. Run after refreshing the snapshot with export_terra_method_configs.py or editing the
table; `--check` exits 1 when the committed graph differs from the regenerated one.
"""

import argparse
import collections
import json
import re
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
SNAPSHOT = REPO_ROOT / "docs" / "method_configs.json"
PIPELINE_MD = REPO_ROOT / "docs" / "pipeline.md"
START, END = "<!-- pipeline-dag:start -->", "<!-- pipeline-dag:end -->"

ATTRIBUTE_RE = re.compile(r"\b(?:this(?:\.\w+)+|workspace\.\w+)")
STEP_RE = re.compile(r"^(\d+)-")
# Data table each `this.<hop>` reference moves to
HOPS = {"samples": "sample", "sample_set": "sample_set", "contigs": "contig"}
SECTION_RE = re.compile(r"^## (\d)\. (.+)$", re.MULTILINE)
# Fill and border of each numbered section's boxes, in section order; white text is legible on every fill
SECTION_STYLES = ["fill:#1d4ed8,stroke:#1e3a8a", "fill:#15803d,stroke:#14532d", "fill:#c2410c,stroke:#7c2d12"]
INIT = '%%{init: {"themeVariables": {"fontSize": "18px"}, "flowchart": {"nodeSpacing": 30, "rankSpacing": 60}}}%%'

# The SV integration workflows hand their VCFs on through bucket directories rather than data table columns
DIRECTORY_EDGES = [
    ("03-MergeSampleSVCallsets", "04-MergeCohortSVCallsets"),
    ("03-MergeSampleSVCallsets", "04-MergeCohortLongSVCallsets"),
    ("04-MergeCohortSVCallsets", "05-RegenotypeSVCallsets"),
    ("05-RegenotypeSVCallsets", "06-MergeRegenotypedSVCallsets"),
    ("06-MergeRegenotypedSVCallsets", "07-PreprocessVcfs"),
    ("04-MergeCohortLongSVCallsets", "07-PreprocessVcfs"),
]

# Cohort inputs filled in by hand rather than written by a config; workspace attributes are all references
HAND_FILLED = {
    ("sample", "unaligned_bam"),
    ("sample", "sex"),
    ("sample", "sample_id"),
    ("sample_set", "sample_set_id"),
    ("sample_set", "sv_integration_dir"),
    ("sample_set", "sv_regenotype_dir"),
    ("sample_set", "sv_vcf"),
    ("sample_set", "sv_vcf_idx"),
    ("sample_set", "sv_long_vcf"),
    ("sample_set", "sv_long_vcf_idx"),
    ("contig", "contig_id"),
    ("contig", "dbsnp_vcf"),
    ("contig", "dbsnp_vcf_idx"),
    ("contig", "dbvar_vcf"),
    ("contig", "dbvar_vcf_idx"),
    ("contig", "drop_variant_ids"),
    ("contig", "gnomad_snv_indel_vcf"),
    ("contig", "gnomad_snv_indel_vcf_idx"),
}


def resolve(expression, root):
    """(data table, attribute) an attribute expression refers to, from a config rooted at `root`."""
    parts = expression.split(".")
    if parts[0] == "workspace":
        return "workspace", parts[1]
    table = root
    for hop in parts[1:-1]:
        table = HOPS[hop]
    return table, parts[-1]


def step(name):
    return int(STEP_RE.match(name).group(1))


def graph(configs):
    """Edges between config names, and a list of wiring errors."""
    writers = collections.defaultdict(set)
    for config in configs:
        for expression in config["outputs"].values():
            writers[resolve(expression, config["root_entity_type"])].add(config["name"])

    errors = [f"{t}.{a} is written by {', '.join(sorted(w))}" for (t, a), w in sorted(writers.items()) if len(w) > 1]
    edges = set(DIRECTORY_EDGES)
    for config in configs:
        for expression in config["inputs"].values():
            for reference in ATTRIBUTE_RE.findall(expression):
                key = resolve(reference, config["root_entity_type"])
                if key in writers:
                    edges.update((source, config["name"]) for source in writers[key] if source != config["name"])
                elif key[0] != "workspace" and key not in HAND_FILLED:
                    errors.append(f"{config['name']} reads {key[0]}.{key[1]}, which no config writes")
    errors += [f"{s} feeds {t} but does not have a lower step number" for s, t in sorted(edges) if step(s) >= step(t)]
    return sorted(edges), errors


def table_configs(markdown):
    """Config names listed in the Run Order table; TRGT stands for every `NN-TRGT-<catalog>` config."""
    names = set()
    for row in re.findall(r"^\| (\d\d) \|(.*)\|$", markdown, re.MULTILINE):
        number, cells = row
        for cell in cells.split("|"):
            for workflow in cell.split(","):
                workflow = re.sub(r"\s*\(.*\)", "", workflow).strip()
                if workflow and workflow != "-":
                    names.add(f"{number}-{workflow}")
    return names


def table_errors(configs, markdown):
    listed = table_configs(markdown)
    names = {re.sub(r"^(\d\d-TRGT)-.*", r"\1", c["name"]) for c in configs}
    return [f"{n} is in the Run Order table but not in the snapshot" for n in sorted(listed - names)] + [
        f"{n} is in the snapshot but not in the Run Order table" for n in sorted(names - listed)
    ]


def sections(configs, markdown):
    """Section number of each config, from the numbered section linking its workflow, and a list of errors."""
    starts = [(m.start(), int(m.group(1))) for m in SECTION_RE.finditer(markdown)]
    found, errors = {}, []
    for config in configs:
        links = re.finditer(r"\*\*\[" + re.escape(config["workflow"]) + r"\]\(", markdown)
        numbers = {max((n for start, n in starts if start < m.start()), default=None) for m in links} - {None}
        if len(numbers) == 1:
            found[config["name"]] = numbers.pop()
        else:
            errors.append(f"{config['workflow']} ({config['name']}) is linked from {len(numbers)} numbered sections, not 1")
    return found, errors


def node(name):
    return "s" + re.sub(r"\W", "_", name)


def mermaid(configs, edges, section_of):
    lines = ["```mermaid", INIT, "flowchart TB"]
    lines += [f"    classDef section{i} {style},color:#ffffff,stroke-width:2px" for i, style in enumerate(SECTION_STYLES, 1)]
    for config in configs:
        lines.append(f'    {node(config["name"])}["{config["name"]}"]:::section{section_of[config["name"]]}')
    lines += [f"    {node(s)} --> {node(t)}" for s, t in edges]
    # Mermaid numbers links in order of definition, so the real edges are the first len(edges)
    lines.append(f"    linkStyle {','.join(map(str, range(len(edges))))} stroke:#334155,stroke-width:2px")
    # Invisible links from each step's first config to every config of the next step keep the steps in rows
    steps = collections.defaultdict(list)
    for config in configs:
        steps[step(config["name"])].append(config["name"])
    numbers = sorted(steps)
    for previous, current in zip(numbers, numbers[1:]):
        lines += [f"    {node(steps[previous][0])} ~~~ {node(name)}" for name in steps[current]]
    return "\n".join(lines + ["```"])


def render(markdown, configs):
    """pipeline.md with the graph between the markers replaced; raises on wiring errors."""
    edges, errors = graph(configs)
    section_of, section_errors = sections(configs, markdown)
    errors += section_errors + table_errors(configs, markdown)
    if errors:
        raise ValueError("\n".join(errors))
    start, end = markdown.index(START) + len(START), markdown.index(END)
    return markdown[:start] + "\n" + mermaid(configs, edges, section_of) + "\n" + markdown[end:]


def main():
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--check", action="store_true", help="exit 1 if docs/pipeline.md differs from regenerated output")
    args = parser.parse_args()

    markdown = PIPELINE_MD.read_text()
    try:
        expected = render(markdown, json.loads(SNAPSHOT.read_text()))
    except ValueError as error:
        print(error)
        return 1
    if args.check:
        if expected != markdown:
            print("out of date: docs/pipeline.md\nRun: python .github/scripts/generate_pipeline_dag.py")
        return int(expected != markdown)
    PIPELINE_MD.write_text(expected)
    print(f"Wrote the Run Order graph to {PIPELINE_MD.relative_to(REPO_ROOT)}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
