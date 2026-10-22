import ast

import pytest
from conftest import REPO_ROOT

# Hail is not installed outside its images, so these scripts are checked statically
HAIL_SCRIPTS = {
    "scripts/helper/split_vcf_hail.py": {"--vcf", "--n-shards", "--records-per-shard", "--prefix", "--cores", "--mem", "--build"},
    "scripts/helper/vep_annotate_hail.py": {"-i", "-o", "--cores", "--mem", "--build"},
    "scripts/annotation/annotate_insilico_predictors.py": {
        "--build",
        "--cores",
        "--mem",
        "--cadd_ht",
        "--pangolin_ht",
        "--phylop_ht",
        "--revel_ht",
        "--spliceai_ht",
        "--vcf",
        "--output_tsv",
    },
}


def parse(rel_path):
    return ast.parse((REPO_ROOT / rel_path).read_text(), filename=rel_path)


def required_arguments(tree):
    found = set()
    for node in ast.walk(tree):
        if isinstance(node, ast.Call) and getattr(node.func, "attr", None) == "add_argument":
            required = any(kw.arg == "required" and getattr(kw.value, "value", False) is True for kw in node.keywords)
            if required:
                found.add(node.args[0].value)
    return found


@pytest.mark.parametrize("rel_path", HAIL_SCRIPTS)
def test_script_defines_main_and_entry_point_guard(rel_path):
    tree = parse(rel_path)
    assert "main" in {node.name for node in tree.body if isinstance(node, ast.FunctionDef)}
    guards = [node for node in tree.body if isinstance(node, ast.If) and "__main__" in ast.dump(node.test)]
    assert guards, "missing if __name__ == '__main__' guard"


@pytest.mark.parametrize("rel_path,expected", HAIL_SCRIPTS.items())
def test_script_declares_required_arguments(rel_path, expected):
    assert required_arguments(parse(rel_path)) == expected
