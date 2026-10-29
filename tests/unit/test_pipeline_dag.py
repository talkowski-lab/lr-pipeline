import json

from conftest import load_script

dag = load_script(".github/scripts/generate_pipeline_dag.py")

TABLE = "| 01 | Align | - | - |\n| 02 | - | Merge | - |\n| 03 | - | - | Annotate |\n"


def config(name, root, inputs=None, outputs=None):
    return {"name": name, "workflow": name[3:], "root_entity_type": root, "inputs": inputs or {}, "outputs": outputs or {}}


def wiring():
    return [
        config("01-Align", "sample", {"W.bam": "this.unaligned_bam"}, {"W.out": "this.aligned_bam"}),
        config("02-Merge", "sample_set", {"W.bams": "this.samples.aligned_bam"}, {"W.out": "workspace.cohort_vcf"}),
        config("03-Annotate", "contig", {"W.vcf": "workspace.cohort_vcf", "W.ref": "workspace.ref_fa"}),
    ]


def test_committed_graph_matches_generator():
    markdown = dag.PIPELINE_MD.read_text()
    assert dag.render(markdown, json.loads(dag.SNAPSHOT.read_text())) == markdown, (
        "run python .github/scripts/generate_pipeline_dag.py"
    )


def test_resolve_follows_entity_hops():
    assert dag.resolve("this.samples.cpg_hap1_bed", "sample_set") == ("sample", "cpg_hap1_bed")
    assert dag.resolve("this.sample_set.failed_bins_bed", "contig") == ("sample_set", "failed_bins_bed")
    assert dag.resolve("workspace.integrated_vcf", "sample") == ("workspace", "integrated_vcf")


def test_edges_follow_written_attributes():
    edges, errors = dag.graph(wiring())
    assert errors == []
    assert ("01-Align", "02-Merge") in edges and ("02-Merge", "03-Annotate") in edges


def test_backward_edge_is_an_error():
    configs = wiring()
    configs[0]["inputs"]["W.vcf"] = "workspace.cohort_vcf"
    assert any("does not have a lower step number" in e for e in dag.graph(configs)[1])


def test_second_writer_is_an_error():
    configs = wiring() + [config("03-Other", "sample_set", outputs={"W.out": "workspace.cohort_vcf"})]
    assert any("is written by 02-Merge, 03-Other" in e for e in dag.graph(configs)[1])


def test_unwritten_read_is_an_error():
    configs = wiring()
    configs[2]["inputs"]["W.bed"] = "this.missing_bed"
    assert "03-Annotate reads contig.missing_bed, which no config writes" in dag.graph(configs)[1]


def test_table_must_match_snapshot():
    assert dag.table_errors(wiring(), TABLE) == []
    assert dag.table_errors(wiring()[:2], TABLE) == ["03-Annotate is in the Run Order table but not in the snapshot"]


def test_sections_come_from_workflow_links():
    markdown = "## 1. Calling\n**[Align](a.wdl)**\n## 2. Merging\n**[Merge](m.wdl)**\n## 3. Annotation\n"
    found, errors = dag.sections(wiring(), markdown)
    assert found == {"01-Align": 1, "02-Merge": 2}
    assert errors == ["Annotate (03-Annotate) is linked from 0 numbered sections, not 1"]
