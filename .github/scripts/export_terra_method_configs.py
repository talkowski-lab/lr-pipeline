#!/usr/bin/env python3
"""Export a template Terra workspace's method configs into its method config snapshot.

`--pipeline main` exports the main template workspace into docs/method_configs.json and `--pipeline cnvs` the
CNV template workspace into docs/method_configs_cnvs.json. Keeps each config's name, workflow, root entity type
(null for a config run without one), inputs and outputs - the wiring that generate_pipeline_dag.py turns into the
Run Order graph in docs/pipeline.md or docs/pipeline_cnvs.md. Method versions are left out, so repointing a config
between a branch and `main` never changes the snapshot. The unnumbered configs, for the optional analyses, are kept
so a clone can be rebuilt from the snapshot; step 00 holds the template-only downsampling and alignment steps and is
left out. Needs `gcloud` authenticated with read access to the workspace.
"""

import argparse
import json
import re
import subprocess
import sys
import urllib.parse
import urllib.request
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
# Template workspace and snapshot of each pipeline
PIPELINES = {
    "main": ("LR_GNOMAD-AoU_TALK_Pipeline-Template", REPO_ROOT / "docs" / "method_configs.json"),
    "cnvs": ("LR_GNOMAD-AoU_TALK_CNV-Pipeline-Template", REPO_ROOT / "docs" / "method_configs_cnvs.json"),
}
RAWLS = "https://api.firecloud.org/api/workspaces"
TEMPLATE_ONLY_RE = re.compile(r"^00-")


def get(url, token):
    request = urllib.request.Request(url, headers={"Authorization": f"Bearer {token}"})
    with urllib.request.urlopen(request) as response:
        return json.load(response)


def export(namespace, name):
    token = subprocess.run(["gcloud", "auth", "print-access-token"], capture_output=True, text=True, check=True).stdout.strip()
    base = f"{RAWLS}/{urllib.parse.quote(namespace)}/{urllib.parse.quote(name)}"
    configs = []
    for listed in get(f"{base}/methodconfigs?allRepos=true", token):
        if TEMPLATE_ONLY_RE.match(listed["name"]):
            continue
        config = get(
            f"{base}/method_configs/{urllib.parse.quote(listed['namespace'])}/{urllib.parse.quote(listed['name'])}", token
        )
        configs.append(
            {
                "name": config["name"],
                "workflow": config["methodRepoMethod"]["methodPath"].rsplit("/", 1)[-1],
                "root_entity_type": config.get("rootEntityType"),
                "inputs": dict(sorted(config["inputs"].items())),
                "outputs": dict(sorted(config["outputs"].items())),
            }
        )
    return sorted(configs, key=lambda c: c["name"])


def main():
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--pipeline", choices=sorted(PIPELINES), default="main")
    parser.add_argument("--workspace-namespace", default="LR_GNOMAD_1_CO-AoU_TALK")
    parser.add_argument("--workspace-name", help="overrides the pipeline's template workspace")
    parser.add_argument("--output", type=Path, help="overrides the pipeline's snapshot")
    args = parser.parse_args()

    workspace, snapshot = PIPELINES[args.pipeline]
    output = args.output or snapshot
    configs = export(args.workspace_namespace, args.workspace_name or workspace)
    output.write_text(json.dumps(configs, indent=4) + "\n")
    print(f"Wrote {len(configs)} method configs to {output}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
