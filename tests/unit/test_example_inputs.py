import json

import pytest
from conftest import REPO_ROOT, load_script

gen = load_script(".github/scripts/generate_example_inputs.py")
wdl_meta = load_script(".github/scripts/wdl_meta.py")

EXPECTED = gen.example_inputs()


def test_committed_inputs_match_generator():
    committed = {p: p.read_text() for p in gen.INPUTS_DIR.rglob("*.json")}
    assert committed == EXPECTED, "run python .github/scripts/generate_example_inputs.py"


def test_every_registered_workflow_has_a_file():
    assert len(EXPECTED) == len(gen.dockstore_workflows())
    for path in EXPECTED:
        assert path.is_file()
        assert not str(path.relative_to(REPO_ROOT)).startswith("inputs/archive")


@pytest.mark.parametrize("path", sorted(EXPECTED), ids=lambda p: p.stem)
def test_keys_are_declared_inputs(path):
    workflow = wdl_meta.parse_workflow(next(w for w in gen.dockstore_workflows() if w.stem == path.stem))
    declared = {f"{workflow.name}.{d.name}" for d in workflow.inputs}
    assert set(json.loads(path.read_text())) <= declared
