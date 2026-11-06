#!/usr/bin/env python3
"""SD-aware breakpoint refinement and scoring for CNV calls with one or both breakpoints inside segmental duplications
(SDs), or lying within an SD. Depth-only; no family / inheritance information is used.

1. Classify each call (original coordinates) against merged SD intervals (--segdup):
     within_SD   >= --within-frac of the call covered by SDs
     2bp_SD      start and end both inside SDs (and not within_SD)
     1bp_SD      exactly one end inside an SD (and not within_SD)
   calls with no end in an SD and < --within-frac coverage are skipped (use refine_cnv_breakpoints.py for those).
2. Merge calls of one event (same SVTYPE):
     2bp_SD     same SD at the start and the same SD at the end
     1bp_SD     same SD at the SD end, and unique ends within --unique-merge-dist bp (single linkage)
     within_SD  same SD block and reciprocal overlap >= --within-ro (single linkage)
   Event carriers = union of the members' carriers; the event profile = per-bin median of the carriers' CN.
3. Per carrier, CN = RD / median coverage x 2, rescaled so its median over the outer flanks is 2. Outer flank =
   --flank bp of sequence immediately outside the SD (for an SD end) or beyond end + 2 kb (for a unique end).
4. Each SD end: plateau = median CN over the unique sequence between the SD's inner edge and the other end (or, for
   2bp_SD, between the two SDs); step = plateau - outer flank level. The breakpoint is the 50% point between the
   flank and plateau levels inside the SD, taken as the edge after the last bin (outside -> inside, 1 kb running
   median) still on the flank side, so later bins stay on the plateau side; the 25% / 75% points give the interval.
   The end passes if the step is in the expected direction with |step| >= --min-step and >= --min-plateau-frac of
   the plateau bins lie on the plateau side of the 50% level.
   A unique end (1bp_SD) is scored per member call (pooled over that call's carriers, at that call's end; members
   of one event can have different unique ends). The end is first moved to the 100 bp bin edge within +/- 2 kb that
   maximises the expected-direction step between the 1 kb outside and inside, then scored as in
   compute_breakpoint_sharpness.py: step, J2 and edge_z (background = this end's outer-flank bins outside SDs);
   sharp if edge_z >= --unique-edge-z, |step| >= --min-step and J2 >= --unique-J2. The end passes if the sharp
   members carry >= 50% of the event's carriers (<side>_ok_carrier_frac); reported values = best member.
   within_SD: plateau = median CN over the call; flank = outside the SD block
   (--local-flank: 2 kb outside the call); pass if the step is in the expected
   direction with |step| >= --min-step and the plateau-fraction test holds. Breakpoints stay at the merged call.
5. Status: SD_PASS (all ends pass), SD_PARTIAL (one of two passes), SD_FAIL. Flags: NEAR_CHROM_END (within --end-mask
   of a chromosome end), NO_PLATEAU (< 20 unique plateau bins), SR_SITE_MATCH (any member matches a short-read call).

--require-sharp: SD ends and within_SD ends also need a sharp transition (unique-end test, within +/- 2 kb of the
   call end); the breakpoint is then that sharp edge.
--sharp-rule majority|score: combine the sharp checks of both ends (see --help); default all = unchanged.
--no-merge: no step 2; every call is its own event (per-sample mode: give a BED with one carrier per row).
Outputs: <prefix>.events.tsv.gz (one row per merged event) and <prefix>.members.tsv.gz (call -> event).

Usage:
  refine_sd_cnv_breakpoints.py --bed chr6.gCNV_2k.merged_cnvs.per_sample_calls.bed.gz --rd hprc_hgsvc.chr6.RD.txt.gz \\
      --samples hprc_hgsvc.medianCov.transposed.bed --medians hprc_hgsvc.medianCov.transposed.bed \\
      --segdup hg38.SegDup.sorted.merged.bed.gz --fai GRCh38.fa.fai --sr-bed 1KGP_sr_hg38.chr6.DEL_DUP_CNV.bed.gz \\
      --out-prefix chr6.gCNV_2k.sd_refined [--min-size 4000]
"""

import argparse
import bisect
import gzip
import subprocess
import sys
from collections import defaultdict

import numpy as np
import pandas as pd

TYPE_MATCH = {"DEL": {"DEL", "CNV"}, "DUP": {"DUP", "CNV"}}


def parse_args():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--bed", required=True, help="per-sample calls: #chrom start end SVID carriers ... SVTYPE")
    p.add_argument("--rd", required=True)
    p.add_argument("--samples", help="sample IDs in RD column order (if --rd has no '#' header)")
    p.add_argument("--medians", required=True)
    p.add_argument("--segdup", required=True)
    p.add_argument("--fai", help="reference .fai, for chromosome lengths (NEAR_CHROM_END flag)")
    p.add_argument("--sr-bed", help="short-read SV BED (#chrom start end ID SVTYPE FILTER) for the SR_SITE_MATCH flag")
    p.add_argument("--out-prefix", required=True)
    p.add_argument("--min-size", type=int, default=4000)
    p.add_argument("--within-frac", type=float, default=0.9)
    p.add_argument("--unique-merge-dist", type=int, default=4000)
    p.add_argument("--within-ro", type=float, default=0.5)
    p.add_argument("--flank", type=int, default=10000, help="outer flank length (bp; default 10000)")
    p.add_argument("--min-step", type=float, default=0.7)
    p.add_argument("--unique-edge-z", type=float, default=4.0, help="edge_z cutoff at unique ends (default 4)")
    p.add_argument("--unique-J2", type=float, default=0.7, help="J2 cutoff at unique ends (default 0.7)")
    p.add_argument("--min-plateau-frac", type=float, default=0.8)
    p.add_argument("--end-mask", type=int, default=500000, help="flag events within this distance of a chromosome end")
    p.add_argument(
        "--no-merge",
        action="store_true",
        help="score every call as its own event (no merging of calls; use with a one-sample-per-row BED)",
    )
    p.add_argument(
        "--require-sharp",
        action="store_true",
        help="SD ends and within_SD ends must also show a sharp transition: the end is moved to the best "
        "100 bp edge within +/- 2 kb of the call end and must have edge_z >= --unique-edge-z, "
        "|step| >= --min-step and J2 >= --unique-J2 (as for unique ends); breakpoint = that edge",
    )
    p.add_argument(
        "--sharp-rule",
        choices=["all", "majority", "score"],
        default="all",
        help="with --require-sharp, how the 3 sharp checks (edge_z, J2, |step|) of the 2 ends are combined: "
        "all = every check at every end (default); majority = each end passes >= 2 of 3 and each check "
        "passes at >= 1 end; score = mean of the 6 value/cutoff ratios (capped at 2) >= 1 and none < 0.8. "
        "Cutoffs are unchanged; the level test (SD / within_SD ends) and the step direction must pass "
        "at both ends in every rule",
    )
    p.add_argument(
        "--local-flank",
        action="store_true",
        help="within_SD: outer flanks = --flank bp starting 2 kb outside the call, instead of outside the SD block",
    )
    if len(sys.argv) == 1:
        p.print_help()
        sys.exit(1)
    return p.parse_args()


def open_text(path):
    return gzip.open(path, "rt") if path.endswith((".gz", ".bgz")) else open(path)


def read_ids(path):
    lines = [x.rstrip("\n") for x in open(path) if x.strip()]
    return [x for x in lines[0].split("\t") if x] if "\t" in lines[0] else [x.strip() for x in lines]


def load_merged(path):
    iv = defaultdict(list)
    with open_text(path) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) >= 3 and f[1].isdigit():
                iv[f[0]].append((int(f[1]), int(f[2])))
    out = {}
    for c, v in iv.items():
        v.sort()
        m = [list(v[0])]
        for a, b in v[1:]:
            if a <= m[-1][1]:
                m[-1][1] = max(m[-1][1], b)
            else:
                m.append([a, b])
        out[c] = [tuple(x) for x in m]
    return out


def sd_index(sd, c, pos):
    """index of the SD interval containing pos, or None"""
    v = sd.get(c, [])
    i = bisect.bisect_right(v, (pos, float("inf"))) - 1
    return i if i >= 0 and v[i][0] <= pos < v[i][1] else None


def sd_cov(sd, c, s, e):
    v = sd.get(c, [])
    i = max(0, bisect.bisect_right(v, (s, float("inf"))) - 1)
    t = 0
    while i < len(v) and v[i][0] < e:
        t += max(0, min(e, v[i][1]) - max(s, v[i][0]))
        i += 1
    return t


def running_median(v, half=5):
    return np.array([np.median(v[max(0, i - half) : i + half]) for i in range(len(v))])


def sd_crossing(st, x, sd_iv, side, lv_out, lv_in, q):
    """breakpoint inside the SD: edge after the last bin (walking outside -> inside) still on the flank side of the
    q point between lv_out and lv_in"""
    m = (st >= sd_iv[0]) & (st < sd_iv[1])
    pos, v = st[m], x[m]
    if len(v) == 0:
        return None
    v = running_median(v)
    if side == "R":
        pos, v = pos[::-1], v[::-1]
    thr = lv_out + q * (lv_in - lv_out)
    idx = np.where((v - thr) * np.sign(lv_out - lv_in) > 0)[0]
    if len(idx) == 0:
        return int(pos[0]) if side == "L" else int(pos[0]) + 100
    k = idx[-1]
    return int(pos[k]) + 100 if side == "L" else int(pos[k])


def refine_unique_end(st, x, bp0, side, sign, search=2000, seg=1000):
    """bin edge within +/- search of bp0 maximising the expected-direction step between the seg bp outside and inside"""
    best, best_bp = -np.inf, bp0
    for b in range(bp0 - search, bp0 + search + 1, 100):
        lo, hi = x[(st >= b - seg) & (st < b)], x[(st >= b) & (st < b + seg)]
        if len(lo) < 3 or len(hi) < 3:
            continue
        d = (np.median(hi) - np.median(lo)) * (sign if side == "L" else -sign)
        if d > best:
            best, best_bp = d, b
    return best_bp


def unique_end_metrics(st, x, bp, side, bg, window=2000):
    if side == "L":
        out, inn = x[(st >= bp - window) & (st < bp - window // 2)], x[(st >= bp + window // 2) & (st < bp + window)]
        seq = x[(st >= bp - window) & (st < bp + window)]
    else:
        out, inn = x[(st >= bp + window // 2) & (st < bp + window)], x[(st >= bp - window) & (st < bp - window // 2)]
        seq = x[(st >= bp - window) & (st < bp + window)][::-1]
    if len(out) < 3 or len(inn) < 3:
        return np.nan, np.nan, np.nan
    step = np.median(inn) - np.median(out)
    if abs(step) < 1e-6:
        return 0.0, 0.0, 0.0
    sg = np.sign(step)
    d1 = np.diff(seq) * sg
    d2 = (seq[2:] - seq[:-2]) * sg
    return step, d2.max() / abs(step), d1.max() / bg


def sharp_end(st, xp, in_sd, side, bp_call, sign, args):
    """move an end to the best edge within +/- 2 kb and score it with the unique-end sharp test"""
    bp0 = refine_unique_end(st, xp, int(bp_call), side, sign)
    if side == "L":
        own = (st >= bp0 - 2000 - args.flank) & (st < bp0 - 2000)
    else:
        own = (st >= bp0 + 2000) & (st < bp0 + 2000 + args.flank)
    fl_bg = own & ~in_sd if (own & ~in_sd).sum() >= 20 else own
    if fl_bg.sum() < 5:
        return bp0, np.nan, np.nan, np.nan, False
    bg = max(np.percentile(np.abs(np.diff(xp[fl_bg])), 99), 0.02)
    stp_u, j2, ez = unique_end_metrics(st, xp, bp0, side, bg)
    ok_u = bool(np.sign(stp_u) == sign and abs(stp_u) >= args.min_step and j2 >= args.unique_J2 and ez >= args.unique_edge_z)
    return bp0, stp_u, j2, ez, ok_u


def joint_pass(ends, sign, args):
    """evaluate the sharp checks of both ends together (--sharp-rule majority / score); cutoffs unchanged"""
    vals = []
    for d in ends:
        v = (d.get("edge_z"), d.get("J2"), d.get("sharp_step"))
        if not d.get("level_ok") or any(x is None or not np.isfinite(x) for x in v) or np.sign(v[2]) != sign:
            return False
        vals.append((v[0] / args.unique_edge_z, v[1] / args.unique_J2, abs(v[2]) / args.min_step))
    r = np.array(vals)
    if args.sharp_rule == "majority":
        return bool(((r >= 1).sum(axis=1) >= 2).all() and (r >= 1).any(axis=0).all())
    rc = np.minimum(r, 2.0)
    return bool(rc.mean() >= 1 and r.min() >= 0.8)


def main():
    args = parse_args()
    header = subprocess.run(["tabix", "-H", args.rd], capture_output=True, text=True).stdout.splitlines()
    samples = header[-1].split("\t")[3:] if header and header[-1].startswith("#") else read_ids(args.samples)
    med_lines = [x.rstrip("\n").split("\t") for x in open(args.medians) if x.strip()]
    med = dict(zip([x for x in med_lines[0] if x], [float(x) for x in med_lines[1] if x.strip()]))
    scale = np.array([2.0 / med[s] for s in samples], dtype=np.float32)
    col = {s: i for i, s in enumerate(samples)}
    sd = load_merged(args.segdup)
    clen = {}
    if args.fai:
        for line in open(args.fai):
            f = line.split("\t")
            clen[f[0]] = int(f[1])
    sr = defaultdict(list)
    if args.sr_bed:
        with open_text(args.sr_bed) as fh:
            for line in fh:
                if line.startswith("#"):
                    continue
                c, s, e, i, t = line.rstrip("\n").split("\t")[:5]
                sr[(c, t)].append((int(s), int(e), i))
        for k in sr:
            sr[k].sort()

    def sr_ids(c, t, s, e):
        out = []
        for tt in TYPE_MATCH.get(t, ()):
            for ss, se, i in sr.get((c, tt), []):
                if ss >= e:
                    break
                ov = min(e, se) - max(s, ss)
                if ov > 0 and min(ov / (e - s), ov / (se - ss)) >= 0.5:
                    out.append(i)
        return out

    bed = pd.read_csv(args.bed, sep="\t", dtype={"carriers": str}, keep_default_na=False, low_memory=False)
    bed = bed.rename(columns={bed.columns[0]: "chrom"})
    bed = bed[
        bed["SVTYPE"].isin(["DEL", "DUP"])
        & ((bed["end"] - bed["start"]) >= args.min_size)
        & ~bed["carriers"].isin(["", "NA", "."])
    ]
    rows = []
    for r in bed.itertuples():
        L = r.end - r.start
        frac = sd_cov(sd, r.chrom, r.start, r.end) / L
        li, ri = sd_index(sd, r.chrom, r.start), sd_index(sd, r.chrom, r.end - 1)
        if frac >= args.within_frac:
            g = "within_SD"
        elif li is not None and ri is not None:
            g = "2bp_SD"
        elif li is not None or ri is not None:
            g = "1bp_SD"
        else:
            continue
        rows.append(
            dict(
                chrom=r.chrom,
                start=r.start,
                end=r.end,
                SVID=r.SVID,
                SVTYPE=r.SVTYPE,
                carriers=r.carriers,
                sd_frac=frac,
                sd_left=li,
                sd_right=ri,
                group=g,
                sr=sr_ids(r.chrom, r.SVTYPE, r.start, r.end),
            )
        )
    member_cols = ["chrom", "start", "end", "SVID", "SVTYPE", "carriers", "sd_frac", "sd_left", "sd_right", "group", "sr"]
    calls = pd.DataFrame(rows, columns=member_cols)
    if calls.empty:
        pd.DataFrame(columns=["event", "chrom", "group", "SVTYPE", "status"]).to_csv(
            args.out_prefix + ".events.tsv.gz", sep="\t", index=False, compression="gzip"
        )
        pd.DataFrame(columns=member_cols[:-1] + ["event", "sr_site_match"]).to_csv(
            args.out_prefix + ".members.tsv.gz", sep="\t", index=False, compression="gzip"
        )
        print("0 calls with SD involvement", file=sys.stderr)
        return
    print(f"{len(calls)} calls with SD involvement: {calls.group.value_counts().to_dict()}", file=sys.stderr)

    # ---- merge into events
    eid = {}
    n_ev = 0
    if args.no_merge:
        eid = {sv: i + 1 for i, sv in enumerate(calls.SVID)}
    for (c, _t, g), sub in [] if args.no_merge else calls.groupby(["chrom", "SVTYPE", "group"]):
        if g == "2bp_SD":
            for _, x in sub.groupby(["sd_left", "sd_right"]):
                n_ev += 1
                for s in x.SVID:
                    eid[s] = n_ev
        elif g == "1bp_SD":
            sub = sub.assign(
                sd_end=np.where(sub.sd_left.notna(), "L", "R"),
                sd_id=sub.sd_left.fillna(sub.sd_right),
                uniq=np.where(sub.sd_left.notna(), sub.end, sub.start),
            )
            for _, x in sub.groupby(["sd_end", "sd_id"]):
                x = x.sort_values("uniq")
                last = None
                for s, u in zip(x.SVID, x.uniq):
                    if last is None or u - last > args.unique_merge_dist:
                        n_ev += 1
                    eid[s] = n_ev
                    last = u
        else:
            sub = sub.assign(block=[sd_index(sd, c, (a + b) // 2) for a, b in zip(sub.start, sub.end)]).sort_values("start")
            for _, x in sub.groupby("block", dropna=False):
                cur = []
                for r in x.itertuples():
                    if cur:
                        ov = min(r.end, cur[-1][1]) - max(r.start, cur[-1][0])
                        if ov > 0 and min(ov / (r.end - r.start), ov / (cur[-1][1] - cur[-1][0])) >= args.within_ro:
                            eid[r.SVID] = n_ev
                            cur[-1] = (min(r.start, cur[-1][0]), max(r.end, cur[-1][1]))
                            continue
                    n_ev += 1
                    eid[r.SVID] = n_ev
                    cur.append((r.start, r.end))
    calls["event"] = calls.SVID.map(eid)

    # ---- score events
    out = []
    for chrom, csub in calls.groupby("chrom", sort=False):
        print(f"loading {chrom} RD ...", file=sys.stderr)
        proc = subprocess.Popen(["tabix", args.rd, chrom], stdout=subprocess.PIPE)
        rd = pd.read_csv(proc.stdout, sep="\t", header=None, dtype={0: str})
        proc.wait()
        starts = rd[1].to_numpy()
        cn_all = rd.iloc[:, 3:].to_numpy(dtype=np.float32) * scale
        del rd
        sdc = sd.get(chrom, [])
        for ev, x in csub.groupby("event"):
            g, t = x.group.iloc[0], x.SVTYPE.iloc[0]
            sign = -1.0 if t == "DEL" else 1.0
            s, e = int(x.start.min()), int(x.end.max())
            carriers = sorted({c for cs in x.carriers for c in cs.split(",") if c in col})
            flags = []
            if g == "2bp_SD":
                SL, SR = sdc[int(x.sd_left.iloc[0])], sdc[int(x.sd_right.iloc[0])]
                plat = (SL[1], SR[0])
                oL, oR = (SL[0] - args.flank, SL[0]), (SR[1], SR[1] + args.flank)
            elif g == "1bp_SD":
                if x.sd_left.notna().iloc[0]:
                    SL, SR = sdc[int(x.sd_left.iloc[0])], None
                    plat = (SL[1], int(x.end.median()) - 2000)
                    oL, oR = (SL[0] - args.flank, SL[0]), (int(x.end.median()) + 2000, int(x.end.median()) + 2000 + args.flank)
                else:
                    SL, SR = None, sdc[int(x.sd_right.iloc[0])]
                    plat = (int(x.start.median()) + 2000, SR[0])
                    oL = (int(x.start.median()) - 2000 - args.flank, int(x.start.median()) - 2000)
                    oR = (SR[1], SR[1] + args.flank)
            else:
                blk = sd_index(sd, chrom, (s + e) // 2)
                B = sdc[blk] if blk is not None else (s, e)
                SL = SR = None
                plat = (s, e)
                oL, oR = (B[0] - args.flank, B[0]), (B[1], B[1] + args.flank)
                if args.local_flank:
                    oL, oR = (s - 2000 - args.flank, s - 2000), (e + 2000, e + 2000 + args.flank)
            lo = np.searchsorted(starts, oL[0])
            hi = np.searchsorted(starts, oR[1])
            st = starts[lo:hi]
            flm = ((st >= oL[0]) & (st < oL[1])) | ((st >= oR[0]) & (st < oR[1]))
            profc = {}
            for c in carriers:
                v = cn_all[lo:hi, col[c]].astype(float)
                base = np.median(v[flm]) if flm.any() else np.nan
                if np.isfinite(base) and base >= 0.3:
                    profc[c] = v / base * 2
            prof = list(profc.values())
            row = dict(
                event=ev,
                chrom=chrom,
                group=g,
                SVTYPE=t,
                n_members=len(x),
                members=",".join(x.SVID),
                n_carriers=len(carriers),
                carriers_used=len(prof),
                call_start=s,
                call_end=e,
                sr_site_match=bool(any(len(v) for v in x.sr)),
                sr_ids=",".join(sorted({i for v in x.sr for i in v})) or ".",
            )
            if args.fai and chrom in clen and (s < args.end_mask or clen[chrom] - e < args.end_mask):
                flags.append("NEAR_CHROM_END")
            if not prof or flm.sum() < 20:
                row.update(status="SD_FAIL", flags=",".join(flags + ["NO_USABLE_FLANK"]))
                out.append(row)
                continue
            xp = np.median(np.vstack(prof), axis=0)
            pm = (st >= plat[0]) & (st < plat[1])
            if pm.sum() < 20:
                flags.append("NO_PLATEAU")
                pm = (st >= s) & (st < e)
            lvp = np.median(xp[pm])
            lvl = np.median(xp[(st >= oL[0]) & (st < oL[1])])
            lvr = np.median(xp[(st >= oR[0]) & (st < oR[1])])
            in_sd = np.zeros(len(st), bool)
            for a, b in sdc[max(0, bisect.bisect_right(sdc, (int(st[0]), 0)) - 1) :]:
                if a > st[-1]:
                    break
                in_sd |= (st >= a) & (st < b)

            ends = []
            for side, S_iv, lv_out in (("L", SL, lvl), ("R", SR, lvr)):
                stp = lvp - lv_out
                mid = lv_out + 0.5 * stp
                pf = np.mean((xp[pm] - mid) * np.sign(stp) > 0) if stp != 0 else 0.0
                if g == "within_SD" or S_iv is not None:
                    ok = np.sign(stp) == sign and abs(stp) >= args.min_step and pf >= args.min_plateau_frac
                    if S_iv is not None:
                        bp = sd_crossing(st, xp, S_iv, side, lv_out, lvp, 0.5)
                        ci = (
                            sd_crossing(st, xp, S_iv, side, lv_out, lvp, 0.25 if side == "L" else 0.75),
                            sd_crossing(st, xp, S_iv, side, lv_out, lvp, 0.75 if side == "L" else 0.25),
                        )
                        kind = "SD"
                    else:
                        bp, ci, kind = (s if side == "L" else e), (None, None), "within_SD"
                    j2_e = ez_e = st_e = np.nan
                    level_ok = bool(ok)
                    if args.require_sharp:
                        bp_s, st_e, j2_e, ez_e, ok_s = sharp_end(st, xp, in_sd, side, s if side == "L" else e, sign, args)
                        ok = bool(ok and ok_s)
                        bp, ci = bp_s, (None, None)
                    ends.append(
                        dict(
                            side=side,
                            kind=kind,
                            bp=bp,
                            ci=ci,
                            step=stp,
                            plateau_frac=pf,
                            ok=ok,
                            J2=j2_e,
                            edge_z=ez_e,
                            sharp_step=st_e,
                            level_ok=level_ok,
                        )
                    )
                else:
                    # score the unique end per member call (its own carriers and end), so that members with
                    # different unique ends are not smeared together
                    res = []
                    for mr in x.itertuples():
                        mc = [c for c in mr.carriers.split(",") if c in profc]
                        if not mc:
                            continue
                        xm = np.median(np.vstack([profc[c] for c in mc]), axis=0)
                        bp0 = refine_unique_end(st, xm, int(mr.start if side == "L" else mr.end), side, sign)
                        own = (
                            (st >= bp0 - 2000 - args.flank) & (st < bp0 - 2000)
                            if side == "L"
                            else (st >= bp0 + 2000) & (st < bp0 + 2000 + args.flank)
                        )
                        fl_bg = own & ~in_sd if (own & ~in_sd).sum() >= 20 else own
                        if fl_bg.sum() < 5:
                            continue
                        bg = max(np.percentile(np.abs(np.diff(xm[fl_bg])), 99), 0.02)
                        stp_u, j2, ez = unique_end_metrics(st, xm, bp0, side, bg)
                        ok = bool(
                            np.sign(stp_u) == sign
                            and abs(stp_u) >= args.min_step
                            and j2 >= args.unique_J2
                            and ez >= args.unique_edge_z
                        )
                        res.append(dict(bp=bp0, step=stp_u, J2=j2, edge_z=ez, ok=ok, n=len(mc)))
                    n_ok_car = sum(r_["n"] for r_ in res if r_["ok"])
                    ok = bool(res) and n_ok_car >= 0.5 * len(profc)
                    if res:
                        best = max(res, key=lambda r_: (r_["ok"], r_["n"], np.nan_to_num(r_["edge_z"])))
                    else:
                        best = dict(
                            bp=int(x.end.median()) if side == "R" else int(x.start.median()),
                            step=np.nan,
                            J2=np.nan,
                            edge_z=np.nan,
                        )
                    ends.append(
                        dict(
                            side=side,
                            kind="unique",
                            bp=best["bp"],
                            ci=(None, None),
                            step=best["step"],
                            plateau_frac=pf,
                            ok=ok,
                            J2=best["J2"],
                            edge_z=best["edge_z"],
                            ok_carrier_frac=n_ok_car / len(profc),
                            sharp_step=best["step"],
                            level_ok=True,
                        )
                    )
            n_ok = sum(bool(d["ok"]) for d in ends)
            status = "SD_PASS" if n_ok == 2 else "SD_PARTIAL" if n_ok == 1 else "SD_FAIL"
            if g == "within_SD":
                status = "SD_PASS" if n_ok == 2 else "SD_FAIL"
            if args.require_sharp and args.sharp_rule != "all" and len(ends) == 2 and joint_pass(ends, sign, args):
                status = "SD_PASS"
            row.update(
                plateau_cn=round(lvp, 3),
                flank_cn_left=round(lvl, 3),
                flank_cn_right=round(lvr, 3),
                status=status,
                flags=",".join(flags) or ".",
            )
            for d in ends:
                k = d["side"]
                row.update(
                    {
                        f"{k}_kind": d["kind"],
                        f"{k}_bp": d["bp"],
                        f"{k}_ci_lo": d["ci"][0],
                        f"{k}_ci_hi": d["ci"][1],
                        f"{k}_step": round(float(d["step"]), 3) if np.isfinite(d["step"]) else np.nan,
                        f"{k}_plateau_frac": round(float(d["plateau_frac"]), 3),
                        f"{k}_J2": d["J2"],
                        f"{k}_edge_z": d["edge_z"],
                        f"{k}_pass": bool(d["ok"]),
                        f"{k}_ok_carrier_frac": d.get("ok_carrier_frac", np.nan),
                        f"{k}_sharp_step": d.get("sharp_step", np.nan),
                        f"{k}_level_ok": d.get("level_ok", np.nan),
                    }
                )
            out.append(row)
        del cn_all
    ev = pd.DataFrame(out)
    ev.to_csv(args.out_prefix + ".events.tsv.gz", sep="\t", index=False, compression="gzip")
    calls.drop(columns=["sr"]).assign(sr_site_match=[bool(v) for v in calls.sr]).to_csv(
        args.out_prefix + ".members.tsv.gz", sep="\t", index=False, compression="gzip"
    )
    print(f"{len(calls)} calls -> {len(ev)} events; {ev.groupby(['group', 'status']).size().to_dict()}", file=sys.stderr)


if __name__ == "__main__":
    main()
