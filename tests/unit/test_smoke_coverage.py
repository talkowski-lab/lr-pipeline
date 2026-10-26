import ast
import re

from conftest import REPO_ROOT

WDL_DIRS = ["wdl/annotation", "wdl/annotation_utils", "wdl/tools"]

# Workflows without a smoke test. Every WDL workflow in WDL_DIRS needs a test in tests/wdl or an entry here.
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
    "CreateCohortDepthFiles": "needs mosdepth-sized inputs for all 24 contigs",
    "CreateCNVDepthPlots": "needs depth matrices for all 24 contigs and the plotting image",
    "CreateTRGTHistograms": "needs the stranalysis image",
    "Automop": "cleans up a Terra workspace through FISS (automop image)",
    "DeepVariant": "needs the DeepVariant image and an aligned BAM",
    "GLNexus": "needs the GLnexus and hail images and per-sample gVCFs",
    "HiFiCNV": "needs the HiFiCNV image and an aligned BAM",
    "HiPhase": "needs the HiPhase images and an aligned BAM",
    "Hifiasm": "needs the hifiasm image and raw reads",
    "Kanpig": "needs the kanpig image and aligned BAMs",
    "LongReadCNVs": "needs the GATK and sv-pipeline images and depth from aligned BAMs",
    "MethylationProfiling": "needs the cpgtools image and an aligned BAM with methylation tags",
    "MinimapAlignment": "needs the minimap2 image and raw reads",
    "MinimapReadAlignment": "needs the minimap2 image and raw reads",
    "MosDepth": "needs the mosdepth images and an aligned BAM",
    "PALMERAssembly": "needs the PALMER image and aligned assembly BAMs",
    "PALMERDiploid": "needs the PALMER image and an aligned BAM",
    "PAV": "needs the PAV image and phased haplotype assemblies",
    "PBSV": "needs the pbsv image and an aligned BAM",
    "RepeatMasker": "needs the multi-GB Dfam TE Tools image",
    "ResetAlignedReads": "needs the minimap2 image for samtools reset and an aligned BAM",
    "Sawfish": "needs the sawfish image and aligned BAMs",
    "Sniffles": "needs the Sniffles image and an aligned BAM",
    "TRGT": "needs the TRGT image and an aligned BAM",
    "TRGTLPS": "needs the trgt-lps and stranalysis images",
    "Vamos": "needs the vamos image and an aligned BAM",
    "MergeCohortLongSVCallsets": "needs the collaborator SV-integration image",
    "MergeCohortSVCallsets": "needs the collaborator SV-integration image",
    "MergeRegenotypedSVCallsets": "needs the collaborator SV-integration image",
    "MergeSampleSVCallsets": "needs the collaborator SV-integration and XGBoost scoring images",
    "RegenotypeSVCallsets": "needs the collaborator SV-integration image",
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


def covered_workflows():
    """Workflow files referenced as `{CONSTANT}/Name.wdl` or `wdl/<dir>/Name.wdl` in the tests/wdl modules."""
    constants = string_constants(ast.parse((REPO_ROOT / "tests" / "conftest.py").read_text()))
    tested = set()
    for path in sorted((REPO_ROOT / "tests" / "wdl").glob("*.py")):
        source = path.read_text()
        local = constants | string_constants(ast.parse(source))
        for name, stem in re.findall(r"\{(\w+)\}/(\w+)\.wdl", source):
            if local.get(name) in WDL_DIRS:
                tested.add(f"{local[name]}/{stem}.wdl")
        literal = re.findall(r"(wdl/(?:annotation|annotation_utils|tools))/(\w+)\.wdl", source)
        tested.update(f"{directory}/{stem}.wdl" for directory, stem in literal)
    return tested


def test_every_workflow_is_smoke_tested_or_exempt():
    workflows = workflow_files()
    tested = {stem for stem, rel in workflows.items() if rel in covered_workflows()}
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
