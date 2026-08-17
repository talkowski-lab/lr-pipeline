version 1.0

import "../utils/Helpers.wdl"
import "../utils/Structs.wdl"

workflow PopulateSVFormatFields {
    input {
        File cohort_vcf
        File cohort_vcf_idx
        String prefix

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
        String merge_args = "--merge id"
        Boolean fuzzy_match_vcf_to_stats = true
        File? swap_samples

        String utils_docker

        RuntimeAttr? runtime_attr_compute_counts
        RuntimeAttr? runtime_attr_aggregate_counts
        RuntimeAttr? runtime_attr_swap_samples
        RuntimeAttr? runtime_attr_subset_cohort_to_samples
        RuntimeAttr? runtime_attr_extract_sample
        RuntimeAttr? runtime_attr_fill_format_fields
        RuntimeAttr? runtime_attr_merge_vcfs
        RuntimeAttr? runtime_attr_concat_missing
    }

    scatter (i in range(length(sample_ids))) {
        call ComputePerSampleCallerCounts {
            input:
                sample_id = sample_ids[i],
                sv_stats = sample_sv_stats[i],
                prefix = "~{prefix}.~{sample_ids[i]}.caller_counts",
                docker = utils_docker,
                runtime_attr_override = runtime_attr_compute_counts
        }
    }

    call AggregateCallerCounts {
        input:
            per_sample_counts = ComputePerSampleCallerCounts.counts_tsv,
            prefix = "~{prefix}.caller_counts",
            docker = utils_docker,
            runtime_attr_override = runtime_attr_aggregate_counts
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

        call FillFormatFields {
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
                dipcall_vcf = dipcall_vcfs[i],
                dipcall_vcf_idx = dipcall_vcf_idxs[i],
                hapdiff_vcf = hapdiff_vcfs[i],
                hapdiff_vcf_idx = hapdiff_vcf_idxs[i],
                caller_counts_tsv = AggregateCallerCounts.counts_tsv,
                fuzzy_match_vcf_to_stats = fuzzy_match_vcf_to_stats,
                prefix = "~{prefix}.~{sample_ids[i]}.filled",
                docker = utils_docker,
                runtime_attr_override = runtime_attr_fill_format_fields
        }
    }

    call Helpers.MergeVcfs {
        input:
            vcfs = FillFormatFields.filled_vcf,
            vcf_idxs = FillFormatFields.filled_vcf_idx,
            prefix = "~{prefix}.filled",
            extra_args = merge_args,
            docker = utils_docker,
            runtime_attr_override = runtime_attr_merge_vcfs
    }

    call Helpers.ConcatTsvs as ConcatCallerSource {
        input:
            tsvs = FillFormatFields.caller_source_tsv,
            sort_output = false,
            preserve_header = true,
            prefix = "~{prefix}.caller_source",
            docker = utils_docker,
            runtime_attr_override = runtime_attr_concat_missing
    }

    output {
        File sv_filled_vcf = MergeVcfs.merged_vcf
        File sv_filled_vcf_idx = MergeVcfs.merged_vcf_idx
        File sv_caller_counts_tsv = AggregateCallerCounts.counts_tsv
        File sv_caller_source_tsv = ConcatCallerSource.concatenated_tsv
    }
}

task ComputePerSampleCallerCounts {
    input {
        String sample_id
        File sv_stats
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        python3 <<CODE
import gzip

def get_len_bucket(svlen):
    l = abs(int(svlen))
    if l < 50:
        return "<50"
    elif l < 500:
        return "50-500"
    elif l < 5000:
        return "500-5000"
    else:
        return ">5000"

counts = {}

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
        svtype = fields['SVTYPE']
        if svtype not in ('INS', 'DEL'):
            continue
        len_bucket = get_len_bucket(fields['SVLEN'])
        for caller in fields['SUPP'].split(','):
            caller = caller.strip()
            key = (svtype, len_bucket, caller)
            counts[key] = counts.get(key, 0) + 1

with open("~{prefix}.tsv", 'w') as out:
    out.write("TYPE\tLEN\tCALLER\tCOUNT\n")
    for (svtype, len_bucket, caller), count in sorted(counts.items()):
        out.write(f"{svtype}\t{len_bucket}\t{caller}\t{count}\n")
CODE
    >>>

    output {
        File counts_tsv = "~{prefix}.tsv"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 4,
        disk_gb: ceil(size(sv_stats, "GB")) + 5,
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

task AggregateCallerCounts {
    input {
        Array[File] per_sample_counts
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        python3 <<CODE
import csv
from collections import defaultdict

input_files = "~{sep=',' per_sample_counts}".split(',')

agg = defaultdict(int)
for f in input_files:
    with open(f, 'r') as fh:
        reader = csv.DictReader(fh, delimiter='\t')
        for row in reader:
            key = (row['TYPE'], row['LEN'], row['CALLER'])
            agg[key] += int(row['COUNT'])

callers = sorted(set(k[2] for k in agg.keys()))
len_order = ["<50", "50-500", "500-5000", ">5000"]
type_order = ["INS", "DEL"]

with open("~{prefix}.tsv", 'w') as out:
    header = ["TYPE", "LEN"] + [f"COUNT_{c}" for c in callers]
    out.write("\t".join(header) + "\n")
    for svtype in type_order:
        for len_bucket in len_order:
            row = [svtype, len_bucket] + [str(agg.get((svtype, len_bucket, c), 0)) for c in callers]
            out.write("\t".join(row) + "\n")
CODE
    >>>

    output {
        File counts_tsv = "~{prefix}.tsv"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 4,
        disk_gb: ceil(size(per_sample_counts, "GB")) + 5,
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

task FillFormatFields {
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
        File? dipcall_vcf
        File? dipcall_vcf_idx
        File? hapdiff_vcf
        File? hapdiff_vcf_idx
        File caller_counts_tsv
        Boolean fuzzy_match_vcf_to_stats
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
import csv
import pysam
from collections import defaultdict
from math import lgamma, log

NO_AD_CALLERS = {'hapdiff'}

def get_len_bucket(svlen):
    l = abs(int(svlen))
    if l < 50:
        return "<50"
    elif l < 500:
        return "50-500"
    elif l < 5000:
        return "500-5000"
    else:
        return ">5000"

def calc_pls(ref_reads, alt_reads):
    total = ref_reads + alt_reads
    if total == 0:
        return [0, 0, 0]
    
    means = [0.03, 0.50, 0.97]
    priors = [0.33, 0.34, 0.33]
    
    ll_raw = []
    for i in range(3):
        ll = log(priors[i]) + alt_reads * log(means[i]) + ref_reads * log(1.0 - means[i])
        ll_10 = ll / log(10)
        ll_raw.append(ll_10)
    
    max_ll = max(ll_raw)
    return [int(round(-10 * (ll - max_ll))) for ll in ll_raw]

def calc_gq(norm_pls):
    return min(sorted(norm_pls)[1], 99)

def get_ad_from_record(record, caller, sample_name, target_svlen=None):
    fmt = record.samples[sample_name]
    if caller in ('cutesv', 'sniffles'):
        dr = fmt.get('DR')
        dv = fmt.get('DV')
        if dr is not None and dv is not None:
            return (int(dr), int(dv))
    elif caller in ('pbsv', 'sawfish', 'dipcall'):
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

def find_matching_variant(vcf_path, chrom, pos, svtype, window=500):
    vcf = pysam.VariantFile(vcf_path)
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
    finally:
        vcf.close()
    return best

sample_id = "~{sample_id}"

# Load caller counts: (TYPE, LEN_BUCKET) -> {caller: count}
caller_counts = {}
with open("~{caller_counts_tsv}", 'r') as f:
    reader = csv.DictReader(f, delimiter='\t')
    for row in reader:
        key = (row['TYPE'], row['LEN'])
        caller_counts[key] = {k[6:]: int(v) for k, v in row.items() if k.startswith('COUNT_')}

# Load sv_stats by variant ID (excluding pav) and build spatial index for fuzzy matching
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
        stat = {            'stat_id': fields['ID'],
            'stat_pos': int(fields['POS']),            'svtype': fields['SVTYPE'],
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

fuzzy_match = ~{true="True" false="False" fuzzy_match_vcf_to_stats}

caller_files = {
    'cutesv':  ("~{cutesv_vcf}",  "~{cutesv_vcf_idx}"),
    'sniffles': ("~{sniffles_vcf}", "~{sniffles_vcf_idx}"),
    'delly':   ("~{delly_vcf}",   "~{delly_vcf_idx}"),
    'pbsv':    ("~{pbsv_vcf}",    "~{pbsv_vcf_idx}"),
    'sawfish': ("~{sawfish_vcf}", "~{sawfish_vcf_idx}"),
    'dipcall': ("~{dipcall_vcf}", "~{dipcall_vcf_idx}"),
    'hapdiff': ("~{hapdiff_vcf}", "~{hapdiff_vcf_idx}"),
}

caller_vcf_map = {}
for vcf_caller, (vcf_path, idx_path) in caller_files.items():
    if not vcf_path:
        continue
    local_vcf = f"caller_{vcf_caller}.vcf.gz"
    local_idx = f"caller_{vcf_caller}.vcf.gz.tbi"
    os.symlink(vcf_path, local_vcf)
    os.symlink(idx_path, local_idx)
    with pysam.VariantFile(local_vcf) as tmp:
        samples = list(tmp.header.samples)
        if not samples:
            continue
        caller_sample = sample_id if sample_id in samples else samples[0]
    caller_vcf_map[vcf_caller] = (local_vcf, caller_sample)

vcf_in = pysam.VariantFile("~{subset_vcf}")
header = vcf_in.header
if 'AD' not in header.formats:
    header.add_line('##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths for the ref and alt alleles">')
if 'GQ' not in header.formats:
    header.add_line('##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype quality">')
if 'PL' not in header.formats:
    header.add_line('##FORMAT=<ID=PL,Number=G,Type=Integer,Description="Phred-scaled genotype likelihoods">')
if 'EV' not in header.formats:
    header.add_line('##FORMAT=<ID=EV,Number=.,Type=String,Description="Callers supporting this variant in this sample">')
if 'BEV' not in header.formats:
    header.add_line('##FORMAT=<ID=BEV,Number=1,Type=String,Description="Best caller for this variant type and size bucket">')

vcf_out = pysam.VariantFile("~{prefix}.vcf.gz", 'w', header=header)
caller_source_rows = []

for record in vcf_in:
    sample_name = list(record.samples.keys())[0]
    sample_fmt = record.samples[sample_name]

    # pysam re-encodes missing/phased alleles on write (e.g. 0|. → 0/.), so save and restore GT
    orig_gt = sample_fmt['GT']
    orig_phased = sample_fmt.phased

    def write_record():
        sample_fmt['GT'] = orig_gt
        sample_fmt.phased = orig_phased
        vcf_out.write(record)

    # Skip uncalled / ref-only genotypes — pass through untouched
    if not orig_gt or not any(a not in (0, None) for a in orig_gt):
        write_record()
        continue

    stat = sv_stats_map.get(record.id)

    if stat is None and fuzzy_match:
        rec_svtype = record.info.get('SVTYPE', None)
        if isinstance(rec_svtype, tuple):
            rec_svtype = rec_svtype[0]
        if rec_svtype:
            stat = find_matching_stat(record.chrom, record.pos, rec_svtype)

    if stat is None or not stat['callers']:
        caller_source_rows.append((
            sample_id, record.id, '.', '.', '.', '.', '.', '.', '.', 'NO_STATS_MATCH'
        ))
        write_record()
        continue

    stat_id = stat['stat_id']
    stat_pos = stat['stat_pos']
    svtype = stat['svtype']
    len_bucket = get_len_bucket(stat['svlen'])
    callers = stat['callers']
    ev = ','.join(callers)

    bucket_counts = caller_counts.get((svtype, len_bucket), {})
    rankable = [c for c in callers if c not in NO_AD_CALLERS]
    ranked_callers = sorted(rankable, key=lambda c: bucket_counts.get(c, 0), reverse=True)

    ad = None
    bev = None
    raw_id = '.'
    raw_pos = None
    for caller in ranked_callers:
        entry = caller_vcf_map.get(caller)
        if entry is None:
            continue
        vcf_path, caller_sample = entry

        match = find_matching_variant(vcf_path, record.chrom, record.pos, svtype)
        if match is None:
            continue

        ad = get_ad_from_record(match, caller, caller_sample, target_svlen=stat['svlen'])
        if ad is not None:
            bev = caller
            raw_id = match.id if match.id else '.'
            raw_pos = match.pos
            break

    sample_fmt['EV'] = ev
    sample_fmt['BEV'] = bev if bev else '.'

    match_dist_vcf_stats = str(abs(record.pos - stat_pos))
    match_dist_vcf_raw = str(abs(record.pos - raw_pos)) if raw_pos is not None else '.'
    match_dist_stats_raw = str(abs(stat_pos - raw_pos)) if raw_pos is not None else '.'

    if ad is None:
        caller_source_rows.append((
            sample_id, record.id, stat_id, raw_id, ev, '.', match_dist_vcf_stats, match_dist_vcf_raw, match_dist_stats_raw, 'NO_AD'
        ))
    else:
        sample_fmt['AD'] = ad
        pls = calc_pls(ad[0], ad[1])
        sample_fmt['PL'] = tuple(pls)
        sample_fmt['GQ'] = calc_gq(pls)
        caller_source_rows.append((
            sample_id, record.id, stat_id, raw_id, ev, bev, match_dist_vcf_stats, match_dist_vcf_raw, match_dist_stats_raw, 'MATCHED'
        ))

    write_record()

vcf_in.close()
vcf_out.close()

with open("~{prefix}.caller_source.tsv", 'w') as f:
    f.write("SAMPLE\tVARIANT_ID_VCF\tVARIANT_ID_STATS\tVARIANT_ID_RAW\tEV\tBEV\tMATCH_DIST_VCF_STATS\tMATCH_DIST_VCF_RAW\tMATCH_DIST_STATS_RAW\tSTATUS\n")
    for row in caller_source_rows:
        f.write("\t".join(row) + "\n")
CODE

        tabix -p vcf ~{prefix}.vcf.gz
    >>>

    output {
        File filled_vcf = "~{prefix}.vcf.gz"
        File filled_vcf_idx = "~{prefix}.vcf.gz.tbi"
        File caller_source_tsv = "~{prefix}.caller_source.tsv"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 4,
        disk_gb: ceil(size(subset_vcf, "GB") + size(sv_stats, "GB") + size(select_all([cutesv_vcf, sniffles_vcf, delly_vcf, pbsv_vcf, sawfish_vcf, dipcall_vcf, hapdiff_vcf]), "GB")) + 10,
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
