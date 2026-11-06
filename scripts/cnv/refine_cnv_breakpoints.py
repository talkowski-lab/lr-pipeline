#!/usr/bin/env python3
"""Refine CNV breakpoints from 100 bp read-depth bins and score how sharp they are.

For every call with SVLEN >= --min-size (default 4 kb) and at least one carrier:

1. Normalisation (per carrier): CN = RD / median coverage x 2, then rescaled so the sample's median CN over its
   own copy-neutral flanks is 2. The flanks are taken beyond the search windows, one call-length each side:
   [start - window - L, start - window) and [end + window, end + window + L), L = SVLEN.
   No non-carrier reference is used (carriers / non-carriers are not known before genotyping).
2. Transition search at each breakpoint: candidate boundaries are bin edges within +/- --window (2 kb) of the
   called breakpoint. For a ramp width w in --ramps (0,100,200,400,600,800,1000 bp) the bins between the boundary
   and boundary + w (inward) are skipped, and the bins on the outside of the boundary are compared with the bins on
   the inside of the ramp (within the window; the inside never passes the middle of the call):
       step  = median(inside) - median(outside)
       noise = sqrt((MAD_out^2 + MAD_in^2) / 2) x 1.4826, floored at --noise-floor
       z     = step / noise
   A transition is significant if the step has the expected sign (DEL: down, DUP: up), |step| >= --min-step and
   |z| >= --min-z. The smallest w with a significant transition is used, and among boundaries the one with the
   largest |z|. New breakpoint = middle of the ramp; CI = ramp +/- half a bin.
3. Scores. Breakpoint quality: 3 = sharp (w <= 100 bp), 2 = w <= 400, 1 = w <= 1000, 0 = no significant
   transition within 1 kb. Call: bp_score = left + right (0-6), n_sharp = number of quality-3 breakpoints,
   status = PASS (both found) / PARTIAL (one) / FAIL (none).
   The call-level result uses the pooled profile (per-bin median CN across carriers); each carrier is also run on
   its own (per-sample table) and the spread of per-sample breakpoints is reported (flag MULTI_BREAKPOINT if the
   per-sample breakpoints differ by more than --multi-bp-spread).

Inputs: --bed (headered BED with #chrom start end SVID carriers SVTYPE, e.g. vcf_to_bed_with_carriers.py output;
SVLEN taken from an SVLEN column or end - start), --rd (bgzipped + tabix-indexed RD matrix, 100 bp bins, loaded one
contig at a time), --samples (IDs in RD column order if --rd has no '#' header), --medians (2-row file: IDs /
median coverage).

EXTENDED SEARCH (--extend-search): when one breakpoint is sharp (quality 3) and the other is not (quality < 3),
the profile is renormalised on the flank beyond the sharp breakpoint only, and the other end is searched from the
sharp breakpoint + --window out to sharp breakpoint + --extend-mult x SVLEN (or leftwards for a sharp right
breakpoint): the transition search above is run at centres every --window bp, keeping only quality-3 transitions
where (a) the CN between the sharp breakpoint and the candidate is on the CNV side of the flank level by >=
--min-step, and (b) the CN in the --window bp beyond the candidate is within 0.5 of the flank level (back to normal).
The candidate with the largest |z| replaces the weak end (flag EXT_RIGHT / EXT_LEFT).

Outputs: <prefix>.calls.tsv.gz (one row per call: original and refined coordinates, CIs, ramp widths, steps, z,
qualities, bp_score, n_sharp, status, flags), <prefix>.per_sample.tsv.gz (one row per call x carrier) and
<prefix>.refined.bed.gz (input columns with start/end replaced by the refined breakpoints where found, plus
orig_start, orig_end, bp_score, n_sharp, bp_status).

Usage:
  refine_cnv_breakpoints.py --bed chr6.gCNV_2k.merged_cnvs.per_sample_calls.bed.gz --rd hprc_hgsvc.chr6.RD.txt.gz \\
      --samples hprc_hgsvc.medianCov.transposed.bed --medians hprc_hgsvc.medianCov.transposed.bed \\
      --out-prefix chr6.gCNV_2k.per_sample_calls.bp_refined [--min-size 4000] [--window 2000]
"""

import argparse
import subprocess
import sys

import numpy as np
import pandas as pd


def parse_args():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--bed", required=True)
    p.add_argument("--rd", required=True)
    p.add_argument("--samples", help="sample IDs in RD column order (if --rd has no '#' header)")
    p.add_argument("--medians", required=True)
    p.add_argument("--out-prefix", required=True)
    p.add_argument("--min-size", type=int, default=4000, help="minimum SVLEN (default 4000)")
    p.add_argument("--window", type=int, default=2000, help="search +/- this many bp around each breakpoint (default 2000)")
    p.add_argument("--ramps", default="0,100,200,400,600,800,1000", help="ramp widths (bp) tried, smallest first")
    p.add_argument("--min-step", type=float, default=0.4, help="minimum |CN step| (default 0.4)")
    p.add_argument("--min-z", type=float, default=4.0, help="minimum |step| / noise (default 4)")
    p.add_argument("--noise-floor", type=float, default=0.05, help="minimum noise SD in CN units (default 0.05)")
    p.add_argument("--min-side-bins", type=int, default=3, help="minimum bins on each side of a boundary (default 3)")
    p.add_argument(
        "--multi-bp-spread",
        type=int,
        default=300,
        help="flag MULTI_BREAKPOINT if per-sample breakpoints span more than this (bp; default 300)",
    )
    p.add_argument("--min-flank-cn", type=float, default=0.3, help="skip a carrier whose flank median CN is below this")
    p.add_argument(
        "--extend-search",
        action="store_true",
        help="if one breakpoint is sharp (quality 3) and the other is not, search for the other end from the "
        "sharp breakpoint out to 2 x the call length (see EXTENDED SEARCH)",
    )
    p.add_argument("--extend-mult", type=float, default=2.0, help="extended search reaches this x call length (default 2)")
    if len(sys.argv) == 1:
        p.print_help()
        sys.exit(1)
    return p.parse_args()


def read_ids(path):
    lines = [x.rstrip("\n") for x in open(path) if x.strip()]
    if "\t" in lines[0]:
        return [x for x in lines[0].split("\t") if x]
    return [x.strip() for x in lines]


def f2(v):
    return None if v is None else round(float(v), 3)


def mad_sd(x):
    m = np.median(x)
    return 1.4826 * np.median(np.abs(x - m))


def find_transition(starts, cn, bp, side, expect_sign, args, call_mid):
    """Search a transition near breakpoint bp. side='left': outside = lower coordinates, inside = higher.
    Returns dict with new breakpoint, ci, ramp, step, z, quality."""
    binw = int(np.median(np.diff(starts))) if len(starts) > 1 else 100
    w0, w1 = bp - args.window, bp + args.window
    best = None
    for ramp in [int(r) for r in args.ramps.split(",")]:
        # 1) locate the boundary: single change-point with an L1 cost (sum of absolute deviations of each side
        #    from its own median; the ramp bins are left out), restricted to steps in the expected direction
        cand = None
        for b in range(w0 - w0 % binw, w1 + 1, binw):
            if side == "left":
                out_mask = (starts >= w0) & (starts + binw <= b)
                in_mask = (starts >= b + ramp) & (starts + binw <= min(w1, call_mid))
            else:
                # right breakpoint: inside = lower coordinates (within the call), outside = higher
                in_mask = (starts >= max(w0, call_mid)) & (starts + binw <= b - ramp)
                out_mask = (starts >= b) & (starts + binw <= w1)
            if out_mask.sum() < args.min_side_bins or in_mask.sum() < args.min_side_bins:
                continue
            o, i = cn[out_mask], cn[in_mask]
            mo, mi = np.median(o), np.median(i)
            if np.sign(mi - mo) != expect_sign:
                continue
            cost = np.abs(o - mo).sum() + np.abs(i - mi).sum()
            if cand is None or cost < cand["cost"]:
                cand = {"cost": cost, "b": b, "o": o, "i": i, "step": mi - mo}
        if cand is None:
            continue
        # 2) test the located boundary
        o, i, step, b = cand["o"], cand["i"], cand["step"], cand["b"]
        noise = max(args.noise_floor, np.sqrt((mad_sd(o) ** 2 + mad_sd(i) ** 2) / 2))
        z = step / noise
        if abs(step) < args.min_step or abs(z) < args.min_z:
            continue
        if side == "left":
            new_bp, ci = b + ramp // 2, (b - binw // 2, b + ramp + binw // 2)
        else:
            new_bp, ci = b - ramp // 2, (b - ramp - binw // 2, b + binw // 2)
        best = {"bp": new_bp, "ci_lo": ci[0], "ci_hi": ci[1], "ramp": ramp, "step": step, "z": z}
        break
    if best is None:
        return {"bp": None, "ci_lo": None, "ci_hi": None, "ramp": None, "step": None, "z": None, "q": 0}
    best["q"] = 3 if best["ramp"] <= 100 else 2 if best["ramp"] <= 400 else 1
    return best


def extend_end(st, x, anchor, side, sign, L, args):
    """search the other end from a sharp anchor breakpoint out to --extend-mult x L; see EXTENDED SEARCH"""
    reach = int(args.extend_mult * L)
    if side == "right":
        centres = range(anchor + args.window, anchor + reach + 1, args.window)
    else:
        centres = range(anchor - args.window, anchor - reach - 1, -args.window)
    best = None
    for c in centres:
        t = find_transition(st, x, c, side, sign, args, anchor)
        if t["q"] != 3:
            continue
        b = t["bp"]
        if side == "right":
            mid, beyond = x[(st >= anchor) & (st + 100 <= b)], x[(st >= b) & (st < b + args.window)]
        else:
            mid, beyond = x[(st >= b) & (st + 100 <= anchor)], x[(st >= b - args.window) & (st < b)]
        if len(mid) < args.min_side_bins or len(beyond) < args.min_side_bins:
            continue
        if (np.median(mid) - 2) * sign < args.min_step or abs(np.median(beyond) - 2) > 0.5:
            continue
        if best is None or abs(t["z"]) > abs(best["z"]):
            best = t
    return best


def main():
    args = parse_args()
    header = subprocess.run(["tabix", "-H", args.rd], capture_output=True, text=True).stdout.splitlines()
    if header and header[-1].startswith("#"):
        samples = header[-1].split("\t")[3:]
    elif args.samples:
        samples = read_ids(args.samples)
    else:
        sys.exit("ERROR: --rd has no '#' header; pass --samples")
    med_lines = [x.rstrip("\n").split("\t") for x in open(args.medians) if x.strip()]
    med = dict(zip([x for x in med_lines[0] if x], [float(x) for x in med_lines[1] if x.strip()]))
    scale = np.array([2.0 / med[s] for s in samples], dtype=np.float32)
    col = {s: i for i, s in enumerate(samples)}

    bed = pd.read_csv(args.bed, sep="\t", dtype={"carriers": str}, keep_default_na=False, low_memory=False)
    bed = bed.rename(columns={bed.columns[0]: "chrom"})
    if "SVLEN" not in bed.columns:
        bed["SVLEN"] = bed["end"] - bed["start"]
    bed["SVLEN"] = pd.to_numeric(bed["SVLEN"], errors="coerce").fillna(bed["end"] - bed["start"]).abs()
    keep = (bed["SVLEN"] >= args.min_size) & bed["SVTYPE"].isin(["DEL", "DUP"]) & ~bed["carriers"].isin(["", "NA", "."])
    todo = bed[keep]
    print(f"{len(bed)} calls; {len(todo)} DEL/DUP >= {args.min_size} bp with carriers to refine", file=sys.stderr)

    call_rows, samp_rows, refined = [], [], {}
    for chrom, sub in todo.groupby("chrom", sort=False):
        print(f"loading {chrom} RD ...", file=sys.stderr)
        proc = subprocess.Popen(["tabix", args.rd, chrom], stdout=subprocess.PIPE)
        rd = pd.read_csv(proc.stdout, sep="\t", header=None, dtype={0: str})
        proc.wait()
        starts = rd[1].to_numpy()
        cn_all = rd.iloc[:, 3:].to_numpy(dtype=np.float32) * scale
        del rd
        for n_done, r in enumerate(sub.itertuples(), 1):
            s, e = int(r.start), int(r.end)
            L = e - s
            sign = -1 if r.SVTYPE == "DEL" else 1
            reach = int(max(1.0, args.extend_mult) * L) if args.extend_search else L
            lo = np.searchsorted(starts, s - args.window - reach)
            hi = np.searchsorted(starts, e + args.window + reach)
            st = starts[lo:hi]
            fl = ((st >= s - args.window - L) & (st < s - args.window)) | ((st >= e + args.window) & (st < e + args.window + L))
            carriers = [c for c in r.carriers.split(",") if c in col]
            prof, raw = [], []
            for c in carriers:
                x = cn_all[lo:hi, col[c]].astype(float)
                raw.append(x)
                base = np.median(x[fl]) if fl.any() else np.nan
                if not np.isfinite(base) or base < args.min_flank_cn:
                    samp_rows.append([r.SVID, c, "skipped_low_flank", None, None, None, None, None, None])
                    continue
                x = x / base * 2
                prof.append(x)
                lb = find_transition(st, x, s, "left", sign, args, (s + e) / 2)
                rb = find_transition(st, x, e, "right", sign, args, (s + e) / 2)
                samp_rows.append(
                    [r.SVID, c, "ok", round(base, 3), lb["bp"], lb["q"], rb["bp"], rb["q"], (lb["q"] > 0) + (rb["q"] > 0)]
                )
            if not prof:
                call_rows.append(
                    [chrom, s, e, r.SVID, r.SVTYPE, L, len(carriers), 0] + [None] * 12 + [0, 0, 0, 0, "FAIL", "NO_USABLE_CARRIER"]
                )
                continue
            pooled = np.median(np.vstack(prof), axis=0)
            lb = find_transition(st, pooled, s, "left", sign, args, (s + e) / 2)
            rb = find_transition(st, pooled, e, "right", sign, args, (s + e) / 2)
            flags = []
            if args.extend_search and (lb["q"] == 3) != (rb["q"] == 3):
                side = "right" if lb["q"] == 3 else "left"
                anchor = lb["bp"] if side == "right" else rb["bp"]
                ofl = (
                    (st >= anchor - args.window - L) & (st < anchor - args.window)
                    if side == "right"
                    else (st >= anchor + args.window) & (st < anchor + args.window + L)
                )
                xs = []
                for x in raw:
                    base = np.median(x[ofl]) if ofl.sum() >= args.min_side_bins else np.nan
                    if np.isfinite(base) and base >= args.min_flank_cn:
                        xs.append(x / base * 2)
                ext = extend_end(st, np.median(np.vstack(xs), axis=0), anchor, side, sign, L, args) if xs else None
                if ext is not None:
                    if side == "right":
                        rb = ext
                    else:
                        lb = ext
                    flags.append("EXT_RIGHT" if side == "right" else "EXT_LEFT")
            n_found = (lb["q"] > 0) + (rb["q"] > 0)
            status = "PASS" if n_found == 2 else "PARTIAL" if n_found == 1 else "FAIL"
            ps = [x for x in samp_rows[-len(carriers) :] if x[0] == r.SVID and x[2] == "ok"]
            for k, nm in ((4, "LEFT"), (6, "RIGHT")):
                vals = [x[k] for x in ps if x[k] is not None]
                if len(vals) >= 2 and max(vals) - min(vals) > args.multi_bp_spread:
                    flags.append(f"MULTI_BREAKPOINT_{nm}")
            new_s = lb["bp"] if lb["bp"] is not None else s
            new_e = rb["bp"] if rb["bp"] is not None else e
            if new_e <= new_s:
                new_s, new_e = s, e
                flags.append("REFINED_INVERTED_KEPT_ORIGINAL")
            if lb["bp"] is not None and new_s != s:
                flags.append("START_MOVED")
            if rb["bp"] is not None and new_e != e:
                flags.append("END_MOVED")
            call_rows.append(
                [
                    chrom,
                    s,
                    e,
                    r.SVID,
                    r.SVTYPE,
                    L,
                    len(carriers),
                    len(prof),
                    new_s,
                    new_e,
                    lb["ci_lo"],
                    lb["ci_hi"],
                    rb["ci_lo"],
                    rb["ci_hi"],
                    lb["ramp"],
                    rb["ramp"],
                    f2(lb["step"]),
                    f2(rb["step"]),
                    f2(lb["z"]),
                    f2(rb["z"]),
                    lb["q"],
                    rb["q"],
                    lb["q"] + rb["q"],
                    (lb["q"] == 3) + (rb["q"] == 3),
                    status,
                    ",".join(flags) or ".",
                ]
            )
            refined[r.SVID] = (new_s, new_e, lb["q"] + rb["q"], (lb["q"] == 3) + (rb["q"] == 3), status)
            if n_done % 1000 == 0:
                print(f"  {chrom}: {n_done}/{len(sub)} calls", file=sys.stderr)
        del cn_all

    cols = [
        "chrom",
        "orig_start",
        "orig_end",
        "SVID",
        "SVTYPE",
        "SVLEN",
        "n_carriers",
        "n_carriers_used",
        "new_start",
        "new_end",
        "start_ci_lo",
        "start_ci_hi",
        "end_ci_lo",
        "end_ci_hi",
        "start_ramp_bp",
        "end_ramp_bp",
        "start_step",
        "end_step",
        "start_z",
        "end_z",
        "start_q",
        "end_q",
        "bp_score",
        "n_sharp",
        "status",
        "flags",
    ]
    calls = pd.DataFrame(call_rows, columns=cols)
    calls.to_csv(args.out_prefix + ".calls.tsv.gz", sep="\t", index=False, compression="gzip")
    pd.DataFrame(
        samp_rows,
        columns=["SVID", "sample", "state", "flank_median_cn", "new_start", "start_q", "new_end", "end_q", "n_bp_found"],
    ).to_csv(args.out_prefix + ".per_sample.tsv.gz", sep="\t", index=False, compression="gzip")
    out = bed.copy()
    out["orig_start"], out["orig_end"] = out["start"], out["end"]
    out["bp_score"], out["n_sharp"], out["bp_status"] = "NA", "NA", "NOT_TESTED"
    for i, sid in enumerate(out["SVID"]):
        if sid in refined:
            ns, ne, sc, nsh, stt = refined[sid]
            out.iat[i, out.columns.get_loc("start")], out.iat[i, out.columns.get_loc("end")] = ns, ne
            out.iat[i, out.columns.get_loc("bp_score")] = sc
            out.iat[i, out.columns.get_loc("n_sharp")] = nsh
            out.iat[i, out.columns.get_loc("bp_status")] = stt
    out = out.rename(columns={"chrom": "#chrom"})
    out.to_csv(args.out_prefix + ".refined.bed.gz", sep="\t", index=False, compression="gzip")
    print(
        f"{len(calls)} calls refined: {calls.status.value_counts().to_dict()}; n_sharp {calls.n_sharp.value_counts().to_dict()}",
        file=sys.stderr,
    )


if __name__ == "__main__":
    main()
