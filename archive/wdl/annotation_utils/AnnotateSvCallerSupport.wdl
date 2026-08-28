version 1.0

import "../utils/Helpers.wdl"
import "../utils/Structs.wdl"

workflow AnnotateSvCallerSupport {
    input {
        File sv_vcf
        File sv_vcf_idx
        Array[File] kanpig_vcfs
        Array[File] kanpig_vcf_idxs
        Array[String] sample_ids
        String prefix

        Array[File?]? sample_sv_stats
        Array[File?]? cutesv_vcfs
        Array[File?]? cutesv_vcf_idxs
        Array[File?]? sniffles_vcfs
        Array[File?]? sniffles_vcf_idxs
        Array[File?]? delly_vcfs
        Array[File?]? delly_vcf_idxs
        Array[File?]? pbsv_vcfs
        Array[File?]? pbsv_vcf_idxs
        Array[File?]? sawfish_vcfs
        Array[File?]? sawfish_vcf_idxs
        Array[File?]? dipcall_vcfs
        Array[File?]? dipcall_vcf_idxs
        Array[File?]? hapdiff_vcfs
        Array[File?]? hapdiff_vcf_idxs
        
        Int truvari_breakpoint_window = 500
        Float truvari_reciprocal_overlap = 0.0
        Float truvari_sequence_similarity = 0.7
        Float truvari_size_similarity = 0.7

        Boolean fuzzy_match_vcf_to_stats = true
        Int fuzzy_match_breakpoint_window = 500

        Boolean match_gt_kanpig = true
        Boolean match_gt_non_kanpig = true

        File? swap_samples
        File? null_file

        String utils_docker

        RuntimeAttr? runtime_attr_swap_samples
        RuntimeAttr? runtime_attr_subset_cohort_to_samples
        RuntimeAttr? runtime_attr_extract_sample
        RuntimeAttr? runtime_attr_process_sample
        RuntimeAttr? runtime_attr_merge_vcfs
        RuntimeAttr? runtime_attr_concat_tsvs
    }

    if (defined(swap_samples)) {
        call Helpers.SwapSampleIds {
            input:
                vcf = sv_vcf,
                vcf_idx = sv_vcf_idx,
                sample_swap_list = select_first([swap_samples]),
                prefix = "~{prefix}.swapped",
                docker = utils_docker,
                runtime_attr_override = runtime_attr_swap_samples
        }
    }

    File final_sv_vcf = select_first([SwapSampleIds.swapped_vcf, sv_vcf])
    File final_sv_vcf_idx = select_first([SwapSampleIds.swapped_vcf_idx, sv_vcf_idx])

    call Helpers.SubsetVcfToSamples as SubsetCohortToSamples {
        input:
            vcf = final_sv_vcf,
            vcf_idx = final_sv_vcf_idx,
            samples = sample_ids,
            prefix = "~{prefix}.subset",
            docker = utils_docker,
            runtime_attr_override = runtime_attr_subset_cohort_to_samples
    }

    scatter (i in range(length(sample_ids))) {
        call Helpers.SubsetVcfToSamples as ExtractSample {
            input:
                vcf = SubsetCohortToSamples.subset_vcf,
                vcf_idx = SubsetCohortToSamples.subset_vcf_idx,
                samples = [sample_ids[i]],
                filter_to_sample = false,
                prefix = "~{prefix}.~{sample_ids[i]}.subset",
                docker = utils_docker,
                runtime_attr_override = runtime_attr_extract_sample
        }

        File? sv_stats_i      = if defined(sample_sv_stats)   then select_first([sample_sv_stats])[i]   else null_file
        File? cutesv_vcf_i    = if defined(cutesv_vcfs)       then select_first([cutesv_vcfs])[i]       else null_file
        File? cutesv_idx_i    = if defined(cutesv_vcf_idxs)   then select_first([cutesv_vcf_idxs])[i]   else null_file
        File? sniffles_vcf_i  = if defined(sniffles_vcfs)     then select_first([sniffles_vcfs])[i]     else null_file
        File? sniffles_idx_i  = if defined(sniffles_vcf_idxs) then select_first([sniffles_vcf_idxs])[i] else null_file
        File? delly_vcf_i     = if defined(delly_vcfs)        then select_first([delly_vcfs])[i]        else null_file
        File? delly_idx_i     = if defined(delly_vcf_idxs)    then select_first([delly_vcf_idxs])[i]    else null_file
        File? pbsv_vcf_i      = if defined(pbsv_vcfs)         then select_first([pbsv_vcfs])[i]         else null_file
        File? pbsv_idx_i      = if defined(pbsv_vcf_idxs)     then select_first([pbsv_vcf_idxs])[i]     else null_file
        File? sawfish_vcf_i   = if defined(sawfish_vcfs)      then select_first([sawfish_vcfs])[i]      else null_file
        File? sawfish_idx_i   = if defined(sawfish_vcf_idxs)  then select_first([sawfish_vcf_idxs])[i]  else null_file
        File? dipcall_vcf_i   = if defined(dipcall_vcfs)      then select_first([dipcall_vcfs])[i]      else null_file
        File? dipcall_idx_i   = if defined(dipcall_vcf_idxs)  then select_first([dipcall_vcf_idxs])[i]  else null_file
        File? hapdiff_vcf_i   = if defined(hapdiff_vcfs)      then select_first([hapdiff_vcfs])[i]      else null_file
        File? hapdiff_idx_i   = if defined(hapdiff_vcf_idxs)  then select_first([hapdiff_vcf_idxs])[i]  else null_file

        call ProcessSample {
            input:
                sample_id = sample_ids[i],
                subset_vcf = ExtractSample.subset_vcf,
                subset_vcf_idx = ExtractSample.subset_vcf_idx,
                kanpig_vcf = kanpig_vcfs[i],
                kanpig_vcf_idx = kanpig_vcf_idxs[i],
                sv_stats = sv_stats_i,
                cutesv_vcf = cutesv_vcf_i,
                cutesv_vcf_idx = cutesv_idx_i,
                sniffles_vcf = sniffles_vcf_i,
                sniffles_vcf_idx = sniffles_idx_i,
                delly_vcf = delly_vcf_i,
                delly_vcf_idx = delly_idx_i,
                pbsv_vcf = pbsv_vcf_i,
                pbsv_vcf_idx = pbsv_idx_i,
                sawfish_vcf = sawfish_vcf_i,
                sawfish_vcf_idx = sawfish_idx_i,
                dipcall_vcf = dipcall_vcf_i,
                dipcall_vcf_idx = dipcall_idx_i,
                hapdiff_vcf = hapdiff_vcf_i,
                hapdiff_vcf_idx = hapdiff_idx_i,
                truvari_breakpoint_window = truvari_breakpoint_window,
                truvari_reciprocal_overlap = truvari_reciprocal_overlap,
                truvari_sequence_similarity = truvari_sequence_similarity,
                truvari_size_similarity = truvari_size_similarity,
                fuzzy_match_vcf_to_stats = fuzzy_match_vcf_to_stats,
                fuzzy_match_breakpoint_window = fuzzy_match_breakpoint_window,
                match_gt_kanpig = match_gt_kanpig,
                match_gt_non_kanpig = match_gt_non_kanpig,
                prefix = "~{prefix}.~{sample_ids[i]}.added",
                docker = utils_docker,
                runtime_attr_override = runtime_attr_process_sample
        }
    }

    call Helpers.MergeVcfs {
        input:
            vcfs = ProcessSample.processed_vcf,
            vcf_idxs = ProcessSample.processed_vcf_idx,
            prefix = "~{prefix}.added",
            extra_args = "--merge id",
            docker = utils_docker,
            runtime_attr_override = runtime_attr_merge_vcfs
    }

    call Helpers.ConcatTsvs {
        input:
            tsvs = ProcessSample.match_counts_tsv,
            sort_output = false,
            preserve_header = true,
            prefix = "~{prefix}.match_counts",
            docker = utils_docker,
            runtime_attr_override = runtime_attr_concat_tsvs
    }

    output {
        File sv_added_vcf = MergeVcfs.merged_vcf
        File sv_added_vcf_idx = MergeVcfs.merged_vcf_idx
        File sv_match_counts_tsv = ConcatTsvs.concatenated_tsv
    }
}

task ProcessSample {
    input {
        String sample_id
        File subset_vcf
        File subset_vcf_idx
        File kanpig_vcf
        File kanpig_vcf_idx
        File? sv_stats
        File? cutesv_vcf
        File? cutesv_vcf_idx
        File? sniffles_vcf
        File? sniffles_vcf_idx
        File? delly_vcf
        File? delly_vcf_idx
        File? pbsv_vcf
        File? pbsv_vcf_idx
        File? sawfish_vcf
        File? sawfish_vcf_idx
        File? dipcall_vcf
        File? dipcall_vcf_idx
        File? hapdiff_vcf
        File? hapdiff_vcf_idx
        Int truvari_breakpoint_window
        Float truvari_reciprocal_overlap
        Float truvari_sequence_similarity
        Float truvari_size_similarity
        Boolean fuzzy_match_vcf_to_stats
        Int fuzzy_match_breakpoint_window
        Boolean match_gt_kanpig
        Boolean match_gt_non_kanpig
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        python3 <<CODE
import bisect
import gzip
import math
import os
from collections import defaultdict
from math import comb

import pysam
import truvari

SKIP_AD_CALLERS = {"dipcall", "hapdiff"}
PER_CALLER_NAMES = ["cutesv", "sniffles", "delly", "pbsv", "sawfish", "dipcall", "hapdiff"]
EXCLUDED_SUPP = {"pav", "kanpig"}

sample_id = "~{sample_id}"
reciprocal_overlap = float(~{truvari_reciprocal_overlap})
size_similarity = float(~{truvari_size_similarity})
sequence_similarity = float(~{truvari_sequence_similarity})
breakpoint_window = int(~{truvari_breakpoint_window})
sv_stats_path = "~{sv_stats}"
fuzzy_match = ~{true="True" false="False" fuzzy_match_vcf_to_stats}
fuzzy_match_breakpoint_window = int(~{fuzzy_match_breakpoint_window})
match_gt_kanpig = ~{true="True" false="False" match_gt_kanpig}
match_gt_non_kanpig = ~{true="True" false="False" match_gt_non_kanpig}

caller_files = {
    "cutesv":  ("~{cutesv_vcf}",  "~{cutesv_vcf_idx}"),
    "sniffles": ("~{sniffles_vcf}", "~{sniffles_vcf_idx}"),
    "delly":   ("~{delly_vcf}",   "~{delly_vcf_idx}"),
    "pbsv":    ("~{pbsv_vcf}",    "~{pbsv_vcf_idx}"),
    "sawfish": ("~{sawfish_vcf}", "~{sawfish_vcf_idx}"),
    "dipcall": ("~{dipcall_vcf}", "~{dipcall_vcf_idx}"),
    "hapdiff": ("~{hapdiff_vcf}", "~{hapdiff_vcf_idx}"),
}

def is_called(gt):
    return any(a is not None and a > 0 for a in gt)

def is_missing(gt):
    return all(a is None for a in gt)

def count_non_ref(gt):
    return sum(1 for a in gt if a is not None and a > 0)

def calculate_pl(ref_reads, alt_reads):
    if ref_reads + alt_reads == 0:
        return (0, 0, 0)
    means = [0.05, 0.50, 0.95]
    log10 = math.log(10)
    ll = [(alt_reads * math.log(means[i]) + ref_reads * math.log(1.0 - means[i])) / log10 for i in range(3)]
    max_ll = max(ll)
    return tuple(int(round(-10 * (x - max_ll))) for x in ll)

def calculate_gq(pls):
    return min(sorted(pls)[1], 99)

def clear_format(rec, sample, n_alleles, clear_gt=True):
    rec.samples[sample]["GQ"] = None
    rec.samples[sample]["AD"] = tuple(None for _ in range(n_alleles))
    rec.samples[sample]["PL"] = tuple(None for _ in range(comb(n_alleles + 1, 2)))
    rec.samples[sample]["DP"] = None
    if clear_gt:
        rec.samples[sample]["GT"] = tuple(None for _ in range(n_alleles))

def kanpig_key(rec):
    ref = rec.ref.upper() if rec.ref else rec.ref
    alts = tuple(a.upper() for a in rec.alts) if rec.alts else ()
    return (rec.chrom, rec.pos, ref, alts, rec.id)

def get_ad_from_record(record, caller, sample_name, target_svlen=None):
    fmt = record.samples[sample_name]
    if caller in ("cutesv", "sniffles"):
        dr = fmt.get("DR")
        dv = fmt.get("DV")
        if dr is not None and dv is not None:
            return (int(dr), int(dv))
    elif caller in ("pbsv", "sawfish"):
        ad = fmt.get("AD")
        if ad is None or any(x is None for x in ad):
            return None
        ad = tuple(int(x) for x in ad)
        if len(ad) == 2:
            return ad
        if len(ad) > 2 and record.alts and target_svlen is not None:
            ref_len = len(record.ref) if record.ref else 0
            best_idx, best_diff = 0, float("inf")
            for i, alt in enumerate(record.alts):
                diff = abs(abs(len(str(alt)) - ref_len) - abs(int(target_svlen)))
                if diff < best_diff:
                    best_diff, best_idx = diff, i
            return (ad[0], ad[best_idx + 1])
    elif caller == "delly":
        rr = fmt.get("RR")
        rv = fmt.get("RV")
        if rr is not None and rv is not None:
            return (int(rr), int(rv))
    return None

def best_truvari_match(cohort_rec, caller_vf):
    chrom, pos = cohort_rec.chrom, cohort_rec.pos
    best = None
    best_score = -1.0
    try:
        for cand in caller_vf.fetch(chrom, max(0, pos - breakpoint_window), pos + breakpoint_window):
            try:
                m = cohort_rec.match(cand)
            except Exception:
                continue
            if not m.state:
                continue
            score = m.score if m.score is not None else 0.0
            if score > best_score:
                best_score, best = score, cand
    except (ValueError, OSError):
        pass
    return best

tv_params = truvari.VariantParams(
    pctsize=size_similarity,
    pctseq=sequence_similarity,
    pctovl=reciprocal_overlap,
    refdist=breakpoint_window,
    sizemin=0,
    sizefilt=0,
    passonly=False,
    skip_gt=True,
    no_ref="a",
)

# sv_stats: optional per-sample BED listing callers supporting each variant.
sv_stats_map = {}
sv_stats_spatial = defaultdict(list)
if sv_stats_path:
    with gzip.open(sv_stats_path, "rt") as f:
        header = None
        for line in f:
            line = line.rstrip("\n")
            if not line:
                continue
            if header is None:
                header = line.lstrip("#").split("\t")
                continue
            fields = dict(zip(header, line.split("\t")))
            callers = tuple(sorted({c.strip() for c in fields["SUPP"].split(",") if c.strip() and c.strip() not in EXCLUDED_SUPP}))
            stat = {
                "id": fields["ID"],
                "chrom": fields["CHROM"],
                "pos": int(fields["POS"]),
                "end": int(fields["END"]),
                "svlen": abs(int(fields["SVLEN"])),
                "svtype": fields["SVTYPE"],
                "callers": callers,
            }
            sv_stats_map[fields["ID"]] = stat
            sv_stats_spatial[fields["CHROM"]].append((stat["pos"], stat))
    for chrom in sv_stats_spatial:
        sv_stats_spatial[chrom].sort(key=lambda x: x[0])

def find_matching_stat(chrom, pos, end, svlen, svtype):
    # Score each candidate by size similarity + reciprocal overlap + breakpoint closeness, pick the max.
    entries = sv_stats_spatial.get(chrom, [])
    if not entries:
        return None
    lo = bisect.bisect_left(entries, (pos - fuzzy_match_breakpoint_window,))
    rec_span = max(1, end - pos)
    rec_len = max(1, abs(svlen))
    best = None
    best_score = -1.0
    for i in range(lo, len(entries)):
        entry_pos, stat = entries[i]
        if entry_pos > pos + fuzzy_match_breakpoint_window:
            break
        if stat["svtype"] != svtype or not stat["callers"]:
            continue
        stat_span = max(1, stat["end"] - stat["pos"])
        stat_len = max(1, stat["svlen"])
        size_sim = min(rec_len, stat_len) / max(rec_len, stat_len)
        overlap = max(0, min(end, stat["end"]) - max(pos, stat["pos"]))
        recip_overlap = overlap / max(rec_span, stat_span)
        bp_closeness = max(0.0, 1.0 - abs(entry_pos - pos) / fuzzy_match_breakpoint_window)
        score = size_sim + recip_overlap + bp_closeness
        if score > best_score:
            best_score = score
            best = stat
    return best

def callers_for_record(rec):
    if not sv_stats_map:
        return None
    stat = sv_stats_map.get(rec.id)
    if stat is None and fuzzy_match:
        rec_svtype = rec.info.get("SVTYPE", None)
        if isinstance(rec_svtype, tuple):
            rec_svtype = rec_svtype[0] if rec_svtype else None
        if rec_svtype:
            rec_svlen = rec.info.get("SVLEN")
            if isinstance(rec_svlen, tuple):
                rec_svlen = rec_svlen[0] if rec_svlen else None
            rec_svlen = abs(int(rec_svlen)) if rec_svlen is not None else abs(rec.stop - rec.pos)
            stat = find_matching_stat(rec.chrom, rec.pos, rec.stop, rec_svlen, rec_svtype)
    if stat is None:
        return None
    return set(stat["callers"])

# Kanpig: load all records into a dict keyed by (CHROM, POS, REF.upper, ALT.upper, ID).
kp_in = pysam.VariantFile("~{kanpig_vcf}")
kp_sample = sample_id if sample_id in kp_in.header.samples else list(kp_in.header.samples)[0]
kp_index = {}
for r in kp_in:
    kp_index[kanpig_key(r)] = r

# Per-caller: open via truvari with shared params.
caller_handles = {}
for c, (path, idx) in caller_files.items():
    if not path:
        continue
    local_vcf = f"caller_{c}.vcf.gz"
    local_idx = f"caller_{c}.vcf.gz.tbi"
    os.symlink(path, local_vcf)
    os.symlink(idx, local_idx)
    vf = truvari.VariantFile(local_vcf, params=tv_params)
    samples = list(vf.header.samples)
    if not samples:
        continue
    csample = sample_id if sample_id in samples else samples[0]
    caller_handles[c] = (vf, csample)

vcf_in = truvari.VariantFile("~{subset_vcf}", params=tv_params)
header = vcf_in.header
for tag, line in [
    ("AD", '##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths">'),
    ("GQ", '##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype Quality">'),
    ("PL", '##FORMAT=<ID=PL,Number=G,Type=Integer,Description="Phred-scaled genotype likelihoods">'),
    ("DP", '##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Approximate read depth">'),
    ("EV", '##FORMAT=<ID=EV,Number=.,Type=String,Description="Callers supporting this sample call; entries as caller_(refAD_altAD); assembly-based callers (dipcall, hapdiff) omit the AD fields">'),
    ("BEV", '##FORMAT=<ID=BEV,Number=1,Type=String,Description="Representative caller for this sample call: kanpig if in EV, else caller from EV with AD that provides the highest GQ">'),
]:
    if tag not in header.formats:
        header.add_line(line)

out = truvari.VariantFile("~{prefix}.vcf.gz", "w", header=header)

v_count = v_match = v_gt_match = v_gt_ad = v_kp_site = v_kp_gt = v_gt_ad_kp = v_nonkp_site = v_nonkp_gt = v_gt_ad_nonkp = 0

for rec in vcf_in:
    rec.translate(out.header)
    sample_name = list(rec.samples.keys())[0]
    s = rec.samples[sample_name]
    gt = s["GT"]
    n_alleles = len(rec.alts) + 1
    chrom = rec.chrom

    if not gt or not is_called(gt):
        # Case 1: Base ref or missing
        kp_rec = kp_index.get(kanpig_key(rec))
        if kp_rec is not None:
            kp_s = kp_rec.samples[kp_sample]
            kp_gt = kp_s["GT"]
            if not is_called(kp_gt) and not is_missing(kp_gt):
                # Case 1a: Kanpig ref --> set FORMAT fields
                kp_ad = kp_s["AD"]
                s["DP"] = kp_s.get("DP")
                s["AD"] = tuple(int(x) for x in kp_ad)
                pls = calculate_pl(int(kp_ad[0]), int(kp_ad[1]))
                s["PL"] = pls
                s["GQ"] = calculate_gq(pls)
                s["GT"] = tuple(0 for _ in range(n_alleles))
            else:
                # Case 1b: Kanpig missing or non-ref --> clear
                clear_format(rec, sample_name, n_alleles)
        else:
            # No kanpig match --> clear
            clear_format(rec, sample_name, n_alleles, clear_gt=False)

        cur_gt = s["GT"]
        if cur_gt is not None:
            s["GT"] = tuple(sorted(cur_gt, key=lambda a: (a is None, a if a is not None else 0)))
        s.phased = False
    else:
        # Case 2: Base non-ref
        v_count += 1
        base_state = count_non_ref(gt)
        supporters = []
        had_any_match = False
        had_gt_match = False
        had_gt_ad = False
        had_kp_site = False
        had_kp_gt = False
        had_gt_ad_kp = False
        had_nonkp_site = False
        had_nonkp_gt = False
        had_gt_ad_nonkp = False

        # Assess Kanpig match
        kp_dp_val = None
        kp_rec = kp_index.get(kanpig_key(rec))
        if kp_rec is not None:
            had_any_match = True
            had_kp_site = True
            kp_s = kp_rec.samples[kp_sample]
            kp_dp_val = kp_s.get("DP")
            kp_gt = kp_s["GT"]
            kp_gt_match = (count_non_ref(kp_gt) == base_state)
            keep_kp = kp_gt_match if match_gt_kanpig else is_called(kp_gt)
            if keep_kp:
                kp_ad = kp_s["AD"]
                had_gt_match = True
                had_kp_gt = True
                had_gt_ad = True
                had_gt_ad_kp = True
                supporters.append(("kanpig", (int(kp_ad[0]), int(kp_ad[1]))))

        # Assess non-Kanpig match
        target_svlen = rec.info.get("SVLEN")
        if isinstance(target_svlen, tuple):
            target_svlen = target_svlen[0] if target_svlen else None
        allowed_callers = callers_for_record(rec)
        for caller, (vf, csample) in caller_handles.items():
            if allowed_callers is not None and caller not in allowed_callers:
                continue
            match = best_truvari_match(rec, vf)
            if match is None:
                continue
            had_any_match = True
            had_nonkp_site = True
            m_gt = match.samples[csample]["GT"]
            m_gt_match = (count_non_ref(m_gt) == base_state)
            keep_match = m_gt_match if match_gt_non_kanpig else is_called(m_gt)
            if not keep_match:
                continue
            had_gt_match = True
            had_nonkp_gt = True
            if caller in SKIP_AD_CALLERS:
                supporters.append((caller, None))
            else:
                ad = get_ad_from_record(match, caller, csample, target_svlen=target_svlen)
                if ad is not None:
                    had_gt_ad = True
                    had_gt_ad_nonkp = True
                supporters.append((caller, ad))

        # Update counts for TSV
        if had_any_match:
            v_match += 1
        if had_gt_match:
            v_gt_match += 1
        if had_gt_ad:
            v_gt_ad += 1
        if had_kp_site:
            v_kp_site += 1
        if had_kp_gt:
            v_kp_gt += 1
        if had_gt_ad_kp:
            v_gt_ad_kp += 1
        if had_nonkp_site:
            v_nonkp_site += 1
        if had_nonkp_gt:
            v_nonkp_gt += 1
        if had_gt_ad_nonkp:
            v_gt_ad_nonkp += 1

        # Select BEV
        bev = None
        best_ad = None
        best_pls = None
        if supporters:
            supporters.sort(key=lambda x: x[0])
            kanpig_sup = next((sup for sup in supporters if sup[0] == "kanpig"), None)
            if kanpig_sup is not None:
                bev = "kanpig"
                best_ad = kanpig_sup[1]
            else:
                best_gq = -1
                for name, ad in supporters:
                    if name in SKIP_AD_CALLERS or ad is None:
                        continue
                    pls_candidate = calculate_pl(ad[0], ad[1])
                    gq_candidate = calculate_gq(pls_candidate)
                    if gq_candidate > best_gq:
                        best_gq = gq_candidate
                        bev = name
                        best_ad = ad
                        best_pls = pls_candidate

        # Set FORMAT fields
        if supporters:
            ev_entries = []
            for name, ad in supporters:
                if name in SKIP_AD_CALLERS or ad is None:
                    ev_entries.append(name)
                else:
                    ev_entries.append(f"{name}_({ad[0]}_{ad[1]})")
            s["EV"] = tuple(ev_entries)
        
        if bev is not None:
            s["BEV"] = bev
            if best_pls is None:
                best_pls = calculate_pl(best_ad[0], best_ad[1])
            s["AD"] = best_ad
            s["PL"] = best_pls
            s["GQ"] = calculate_gq(best_pls)
            if bev == "kanpig" and kp_dp_val is not None:
                s["DP"] = int(kp_dp_val)
        else:
            clear_format(rec, sample_name, n_alleles, clear_gt=False)

    out.write(rec)

vcf_in.close()
out.close()
kp_in.close()
for vf, _ in caller_handles.values():
    vf.close()

with open("~{prefix}.match_counts.tsv", "w") as fh:
    fh.write("sample_id\tcount\tcount_match_site\tcount_match_gt\tcount_match_gt_ad\tcount_match_site_kp\tcount_match_gt_kp\tcount_match_gt_ad_kp\tcount_match_site_nonkp\tcount_match_gt_nonkp\tcount_match_gt_ad_nonkp\n")
    fh.write(f"{sample_id}\t{v_count}\t{v_match}\t{v_gt_match}\t{v_gt_ad}\t{v_kp_site}\t{v_kp_gt}\t{v_gt_ad_kp}\t{v_nonkp_site}\t{v_nonkp_gt}\t{v_gt_ad_nonkp}\n")
CODE

        tabix -p vcf ~{prefix}.vcf.gz
    >>>

    output {
        File processed_vcf = "~{prefix}.vcf.gz"
        File processed_vcf_idx = "~{prefix}.vcf.gz.tbi"
        File match_counts_tsv = "~{prefix}.match_counts.tsv"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 8,
        disk_gb: ceil(size(subset_vcf, "GB") + size(kanpig_vcf, "GB") + size(select_all([cutesv_vcf, sniffles_vcf, delly_vcf, pbsv_vcf, sawfish_vcf, dipcall_vcf, hapdiff_vcf]), "GB")) + 20,
        boot_disk_gb: 10,
        preemptible_tries: 1,
        max_retries: 0
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
        memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
        disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
        bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
        docker: docker
        preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
    }
}
