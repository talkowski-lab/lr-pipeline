#!/usr/bin/env python3
"""Score breakpoint-sharpness tables (compute_breakpoint_sharpness.py output) with a gradient-boosted tree model
exported to JSON, using numpy only. Adds a p_gbm column: the probability that the call is real, from its two
breakpoints' edge_z, step and J2.

Features (as the model was trained): for each breakpoint, direction-signed log1p(edge_z), step and J2 (signed by
whether the step goes in the expected direction: down for DEL, up for DUP); the weaker and stronger of the two
breakpoints for each measure give the 6 features. Model JSON: {"features", "baseline", "trees": [{"feature",
"threshold", "missing_left", "left", "right", "leaf", "value"}]}; p = sigmoid(baseline + sum of leaf values), going
left when x <= threshold.

Usage:
  score_sharpness_gbm.py --model sharpness_gbm.json --in calls.bp_sharpness.tsv.gz --out calls.scored.tsv.gz
"""

import argparse
import json
import sys

import numpy as np
import pandas as pd


def features(d):
    exp = np.where(d["SVTYPE"] == "DEL", -1.0, 1.0)
    f = pd.DataFrame(index=d.index)
    per = {}
    for side in ("left", "right"):
        step = d[f"{side}_step"].fillna(0) * exp
        agree = np.where(step > 0, 1.0, -1.0)
        per[side] = {
            "log_edge_z": np.log1p(d[f"{side}_edge_z"].fillna(0).clip(lower=0)) * agree,
            "dir_step": step,
            "dir_J2": d[f"{side}_J2"].fillna(0).clip(-1.5, 1.5) * agree,
        }
    for k in ("log_edge_z", "dir_step", "dir_J2"):
        a, b = per["left"][k], per["right"][k]
        f[f"weak_{k}"], f[f"strong_{k}"] = np.minimum(a, b), np.maximum(a, b)
    return f


def predict(model, x):
    raw = np.full(len(x), model["baseline"], dtype=float)
    rows = np.arange(len(x))
    for t in model["trees"]:
        feat, thr = np.array(t["feature"]), np.array(t["threshold"])
        miss, left, right = np.array(t["missing_left"], bool), np.array(t["left"]), np.array(t["right"])
        leaf, value = np.array(t["leaf"], bool), np.array(t["value"])
        node = np.zeros(len(x), dtype=int)
        while True:
            active = ~leaf[node]
            if not active.any():
                break
            v = x[rows, feat[node]]
            go_left = np.where(np.isnan(v), miss[node], v <= thr[node])
            node = np.where(active, np.where(go_left, left[node], right[node]), node)
        raw += value[node]
    return 1 / (1 + np.exp(-raw))


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--model", required=True, help="model JSON")
    p.add_argument("--in", dest="input", required=True, help="compute_breakpoint_sharpness.py table")
    p.add_argument("--out", required=True, help="output table (.gz) with a p_gbm column")
    if len(sys.argv) == 1:
        p.print_help()
        sys.exit(1)
    a = p.parse_args()
    with open(a.model) as fh:
        model = json.load(fh)
    d = pd.read_csv(a.input, sep="\t")
    x = features(d)[model["features"]].to_numpy(dtype=float) if len(d) else np.zeros((0, len(model["features"])))
    d["p_gbm"] = predict(model, x) if len(d) else []
    d.to_csv(a.out, sep="\t", index=False, compression="gzip" if a.out.endswith(".gz") else None)
    print(f"scored {len(d)} calls", file=sys.stderr)


if __name__ == "__main__":
    main()
