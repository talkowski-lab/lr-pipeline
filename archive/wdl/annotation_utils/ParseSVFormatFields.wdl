version 1.0

import "../utils/Helpers.wdl"
import "../utils/Structs.wdl"

workflow ParseSVFormatFields {
    input {
        File cohort_vcf
        File cohort_vcf_idx
        Array[String] sample_ids
        Array[File] sample_sv_stats
        Array[File?] cutesv_vcfs
        Array[File?] cutesv_vcf_idxs
        Array[File?] sniffles_vcfs
        Array[File?] sniffles_vcf_idxs
        Array[File?] delly_vcfs
        Array[File?] delly_vcf_idxs
        Array[File?] pbsv_vcfs
        Array[File?] pbsv_vcf_idxs
        Array[File?] sawfish_vcfs
        Array[File?] sawfish_vcf_idxs
        Array[File?] dipcall_vcfs
        Array[File?] dipcall_vcf_idxs
        Array[File?] hapdiff_vcfs
        Array[File?] hapdiff_vcf_idxs
        String contig
        String prefix

        File? swap_samples

        String utils_docker

        RuntimeAttr? runtime_attr_swap_samples
        RuntimeAttr? runtime_attr_subset_cohort_to_samples
        RuntimeAttr? runtime_attr_extract_sample
        RuntimeAttr? runtime_attr_parse_format_fields
        RuntimeAttr? runtime_attr_concat_tsvs
    }

    if (defined(swap_samples)) {
        call Helpers.SwapSampleIds {
            input:
                vcf = cohort_vcf,
                vcf_idx = cohort_vcf_idx,
                sample_swap_list = select_first([swap_samples]),
                prefix = "~{prefix}.swapped",
                docker = utils_docker,
                runtime_attr_override = runtime_attr_swap_samples
        }
    }

    File final_cohort_vcf = select_first([SwapSampleIds.swapped_vcf, cohort_vcf])
    File final_cohort_vcf_idx = select_first([SwapSampleIds.swapped_vcf_idx, cohort_vcf_idx])

    call Helpers.SubsetVcfToSamples as SubsetCohortToSamples {
        input:
            vcf = final_cohort_vcf,
            vcf_idx = final_cohort_vcf_idx,
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

        call ParseFormatFields {
            input:
                sample_id = sample_ids[i],
                subset_vcf = ExtractSample.subset_vcf,
                subset_vcf_idx = ExtractSample.subset_vcf_idx,
                sv_stats = sample_sv_stats[i],
                cutesv_vcf = cutesv_vcfs[i],
                cutesv_vcf_idx = cutesv_vcf_idxs[i],
                sniffles_vcf = sniffles_vcfs[i],
                sniffles_vcf_idx = sniffles_vcf_idxs[i],
                delly_vcf = delly_vcfs[i],
                delly_vcf_idx = delly_vcf_idxs[i],
                pbsv_vcf = pbsv_vcfs[i],
                pbsv_vcf_idx = pbsv_vcf_idxs[i],
                sawfish_vcf = sawfish_vcfs[i],
                sawfish_vcf_idx = sawfish_vcf_idxs[i],
                prefix = "~{prefix}.~{sample_ids[i]}.gq_calls",
                docker = utils_docker,
                runtime_attr_override = runtime_attr_parse_format_fields
        }
    }

    call Helpers.ConcatTsvs {
        input:
            tsvs = ParseFormatFields.gq_calls_tsv,
            sort_output = false,
            preserve_header = true,
            prefix = "~{prefix}.gq_calls",
            docker = utils_docker,
            runtime_attr_override = runtime_attr_concat_tsvs
    }

    output {
        File gq_calls_tsv = ConcatTsvs.concatenated_tsv
    }
}

task ParseFormatFields {
    input {
        String sample_id
        File subset_vcf
        File subset_vcf_idx
        File sv_stats
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
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        python3 <<CODE
import bisect
import gzip
import os
import pysam
from collections import defaultdict

def get_ad_from_record(record, caller, sample_name, target_svlen=None):
    fmt = record.samples[sample_name]
    if caller in ('cutesv', 'sniffles'):
        dr = fmt.get('DR')
        dv = fmt.get('DV')
        if dr is not None and dv is not None:
            return (int(dr), int(dv))
    elif caller in ('pbsv', 'sawfish'):
        ad = fmt.get('AD')
        if ad is not None:
            if any(x is None for x in ad):
                return None
            ad = tuple(int(x) for x in ad)
            if len(ad) == 2:
                return ad
            if len(ad) > 2 and record.alts and target_svlen is not None:
                ref_len = len(record.ref) if record.ref else 0
                best_idx = 0
                best_diff = float('inf')
                for i, alt in enumerate(record.alts):
                    allele_svlen = abs(len(str(alt)) - ref_len)
                    diff = abs(allele_svlen - abs(int(target_svlen)))
                    if diff < best_diff:
                        best_diff = diff
                        best_idx = i
                return (ad[0], ad[best_idx + 1])
    elif caller == 'delly':
        rr = fmt.get('RR')
        rv = fmt.get('RV')
        if rr is not None and rv is not None:
            return (int(rr), int(rv))
    return None

def rec_has_svtype(rec, svtype):
    try:
        svt = rec.info.get('SVTYPE', '')
        if isinstance(svt, tuple):
            svt = svt[0]
        return svt == svtype
    except ValueError:
        ref_len = len(rec.ref) if rec.ref else 0
        if rec.alts:
            for alt in rec.alts:
                alt_len = len(str(alt))
                if svtype == 'INS' and alt_len > ref_len:
                    return True
                elif svtype == 'DEL' and ref_len > alt_len:
                    return True
        return False

def get_rec_svlen(rec):
    try:
        svl = rec.info.get('SVLEN', None)
        if svl is not None:
            if isinstance(svl, tuple):
                svl = svl[0]
            return abs(svl)
    except (ValueError, TypeError):
        pass
    ref_len = len(rec.ref) if rec.ref else 0
    alt_len = max((len(str(a)) for a in rec.alts), default=0) if rec.alts else 0
    return abs(ref_len - alt_len)

def find_matching_variant(vcf, chrom, pos, svtype, window=500):
    start = max(0, pos - window)
    end = pos + window
    best = None
    best_dist = float('inf')
    try:
        for rec in vcf.fetch(chrom, start, end):
            if not rec_has_svtype(rec, svtype):
                continue
            if get_rec_svlen(rec) < 50:
                continue
            dist = abs(rec.pos - pos)
            if dist < best_dist:
                best_dist = dist
                best = rec
    except ValueError:
        pass
    return best

SKIP_CALLERS = {'hapdiff', 'dipcall'}
sample_id = "~{sample_id}"

# Load sv_stats by variant ID and build a spatial index for fuzzy matching.
sv_stats_map = {}
sv_stats_spatial = defaultdict(list)
with gzip.open("~{sv_stats}", 'rt') as f:
    header = None
    for line in f:
        line = line.strip()
        if not line:
            continue
        if header is None:
            header = line.lstrip('#').split('\t')
            continue
        fields = dict(zip(header, line.split('\t')))
        callers = [c.strip() for c in fields['SUPP'].split(',') if c.strip() != 'pav']
        stat = {
            'stat_id': fields['ID'],
            'stat_pos': int(fields['POS']),
            'svtype': fields['SVTYPE'],
            'svlen': fields['SVLEN'],
            'callers': callers,
        }
        sv_stats_map[fields['ID']] = stat
        sv_stats_spatial[fields['CHROM']].append((int(fields['POS']), stat))
for chrom in sv_stats_spatial:
    sv_stats_spatial[chrom].sort(key=lambda x: x[0])

def find_matching_stat(chrom, pos, svtype, window=500):
    entries = sv_stats_spatial.get(chrom, [])
    if not entries:
        return None
    lo = bisect.bisect_left(entries, (pos - window,))
    best = None
    best_dist = float('inf')
    for i in range(lo, len(entries)):
        entry_pos, stat = entries[i]
        if entry_pos > pos + window:
            break
        if stat['svtype'] != svtype or not stat['callers']:
            continue
        dist = abs(entry_pos - pos)
        if dist < best_dist:
            best_dist = dist
            best = stat
    return best

caller_files = {
    'cutesv':   ("~{cutesv_vcf}",   "~{cutesv_vcf_idx}"),
    'sniffles': ("~{sniffles_vcf}", "~{sniffles_vcf_idx}"),
    'delly':    ("~{delly_vcf}",    "~{delly_vcf_idx}"),
    'pbsv':     ("~{pbsv_vcf}",     "~{pbsv_vcf_idx}"),
    'sawfish':  ("~{sawfish_vcf}",  "~{sawfish_vcf_idx}"),
}

# Open each caller VCF once so pysam reuses its tabix index across all
# per-record fetches (the Populate task re-opened on every call).
caller_vcf_map = {}
for caller, (vcf_path, idx_path) in caller_files.items():
    if not vcf_path:
        continue
    local_vcf = f"caller_{caller}.vcf.gz"
    local_idx = f"caller_{caller}.vcf.gz.tbi"
    os.symlink(vcf_path, local_vcf)
    os.symlink(idx_path, local_idx)
    handle = pysam.VariantFile(local_vcf)
    samples = list(handle.header.samples)
    if not samples:
        handle.close()
        continue
    caller_sample = sample_id if sample_id in samples else samples[0]
    caller_vcf_map[caller] = (handle, caller_sample)

vcf_in = pysam.VariantFile("~{subset_vcf}")
out_path = "~{prefix}.tsv"
with open(out_path, 'w') as out:
    out.write("variant\tsample\tcaller\tref_depth\talt_depth\n")
    for record in vcf_in:
        sample_name = list(record.samples.keys())[0]
        gt = record.samples[sample_name]['GT']
        # Skip uncalled / ref-only genotypes.
        if not gt or not any(a not in (0, None) for a in gt):
            continue

        stat = sv_stats_map.get(record.id)
        if stat is None:
            rec_svtype = record.info.get('SVTYPE', None)
            if isinstance(rec_svtype, tuple):
                rec_svtype = rec_svtype[0]
            if rec_svtype:
                stat = find_matching_stat(record.chrom, record.pos, rec_svtype)
        if stat is None or not stat['callers']:
            continue

        variant_id = record.id if record.id else '.'
        for caller in stat['callers']:
            if caller in SKIP_CALLERS:
                continue
            entry = caller_vcf_map.get(caller)
            if entry is None:
                continue
            vcf_handle, caller_sample = entry
            match = find_matching_variant(vcf_handle, record.chrom, record.pos, stat['svtype'])
            if match is None:
                continue
            ad = get_ad_from_record(match, caller, caller_sample, target_svlen=stat['svlen'])
            if ad is None:
                continue
            out.write(f"{variant_id}\t{sample_id}\t{caller}\t{ad[0]}\t{ad[1]}\n")

vcf_in.close()
for handle, _ in caller_vcf_map.values():
    handle.close()
CODE
    >>>

    output {
        File gq_calls_tsv = "~{prefix}.tsv"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 4,
        disk_gb: ceil(size(subset_vcf, "GB") + size(sv_stats, "GB") + size(select_all([cutesv_vcf, sniffles_vcf, delly_vcf, pbsv_vcf, sawfish_vcf]), "GB")) + 10,
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
