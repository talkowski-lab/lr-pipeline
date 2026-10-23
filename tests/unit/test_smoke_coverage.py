import ast
import re

from conftest import REPO_ROOT

WDL_DIRS = ["wdl/annotation", "wdl/annotation_utils"]

# Workflows without a smoke test. Every WDL workflow in the two directories needs a test in tests/wdl or an entry here.
EXEMPT = {
    "AnnotateInSilicoPredictors": "runs Hail with external predictor tables (hail image, large resources)",
    "AnnotateVEPHail": "runs VEP through Hail with a VEP cache (vephail image, large resources)",
    "AnnotateVEPTRVs": "runs VEP with a VEP cache (vephail image, large resources)",
    "AnnotateIndelTRs": "needs the stranalysis image",
    "AnnotateL1MEAID": "needs the l1meaid, RepeatMasker and intact MEI images and resources",
    "AnnotateSVAN": "needs the SVAN image and repeat resources",
    "AnnotateSVAnnotate": "needs the GATK image and a GENCODE annotation",
    "AnnotateTruvariRemap": "needs the remap and bwa images with a reference index",
    "AnnotateVRS": "needs the VRS image and a SeqRepo database",
    "AnnotatePALMER": "needs PALMER output VCFs from the palmer image",
    "ConvertVcfToBed": "needs svtk from the gatk_sv_lr image",
    "CreateCohortDepthFiles": "needs mosdepth-sized inputs for all 24 contigs",
    "CreateCNVDepthPlots": "needs depth matrices for all 24 contigs and the plotting image",
    "CreateTRGTHistograms": "needs the stranalysis image",
    "SummarizeAnnotationsGATKSV": "needs SV_colors.txt in the utils image, absent until the file is committed",
    "CreateCohortPedigreeAncestryFilesAoUPhase2": "archived",
}


def workflow_files():
    files = {}
    for directory in WDL_DIRS:
        for path in sorted((REPO_ROOT / directory).glob("*.wdl")):
            if re.search(r"^workflow\s+\w+", path.read_text(), re.M):
                files[path.stem] = f"{directory}/{path.name}"
    return files


def string_constants(tree):
    return {
        target.id: node.value.value
        for node in tree.body
        if isinstance(node, ast.Assign) and isinstance(node.value, ast.Constant) and isinstance(node.value.value, str)
        for target in node.targets
        if isinstance(target, ast.Name)
    }


def tested_workflows():
    """Workflow files referenced as `{CONSTANT}/Name.wdl` or `wdl/<dir>/Name.wdl` in the tests/wdl modules."""
    constants = string_constants(ast.parse((REPO_ROOT / "tests" / "conftest.py").read_text()))
    tested = set()
    for path in sorted((REPO_ROOT / "tests" / "wdl").glob("*.py")):
        source = path.read_text()
        local = constants | string_constants(ast.parse(source))
        for name, stem in re.findall(r"\{(\w+)\}/(\w+)\.wdl", source):
            if local.get(name) in WDL_DIRS:
                tested.add(f"{local[name]}/{stem}.wdl")
        tested.update(
            f"{directory}/{stem}.wdl" for directory, stem in re.findall(r"(wdl/annotation(?:_utils)?)/(\w+)\.wdl", source)
        )
    return tested


def test_every_workflow_is_smoke_tested_or_exempt():
    workflows = workflow_files()
    tested = {stem for stem, rel in workflows.items() if rel in tested_workflows()}
    exempt = set(EXEMPT)

    unknown = sorted(exempt - set(workflows))
    assert not unknown, f"EXEMPT lists workflows that do not exist in {WDL_DIRS}: {unknown}"
    both = sorted(exempt & tested)
    assert not both, f"These workflows now have a smoke test; remove them from EXEMPT in {__file__}: {both}"
    missing = sorted(set(workflows) - tested - exempt)
    assert not missing, (
        f"These workflows have no smoke test in tests/wdl. Add a test, or add the workflow to EXEMPT in {__file__} "
        f"with a one-line reason: {missing}"
    )
    assert all(EXEMPT.values()), "Every EXEMPT entry needs a reason"
