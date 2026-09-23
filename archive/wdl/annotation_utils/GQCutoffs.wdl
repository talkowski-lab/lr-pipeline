version 1.0

import "../utils/Helpers.wdl"
import "../utils/Structs.wdl"

workflow GQCutoffs {
    meta {
        description: [
            "This utility derives genotype-quality cutoffs by comparing a callset against trio and truth-set expectations. Trio children are identified from a PED, the callset and truth VCFs are subset to those samples, and precision and recall are tabulated across genotype-quality thresholds and variant length bins to produce a cutoff table."
        ]
    }

    parameter_meta {
        vcfs: "Per-contig callset VCFs to evaluate."
        vcf_idxs: "Index for `vcfs`."
        truth_vcfs: "Truth VCFs the callset is compared against."
        truth_vcf_idxs: "Index for `truth_vcfs`."
        subset_vcf_string: "`bcftools view` arguments used to pre-subset the callset."
        ped: "Six-column PED used to identify trio children."
        swap_samples_truth: "Optional sample-swap list applied to the truth VCFs."
        skip_trv: "Whether to skip tandem-repeat variants."
        length_bins: "Allele-length bin boundaries defining the size buckets."
        min_length_heuristic_comparison: "Minimum variant length at which heuristic matching replaces exact matching."
        del_size_similarity: "Minimum size similarity for matching deletions."
        del_reciprocal_overlap: "Minimum reciprocal overlap for matching deletions."
        del_breakpoint_window: "Breakpoint window, in bp, for matching deletions."
        ins_size_similarity: "Minimum size similarity for matching insertions."
        ins_breakpoint_window: "Breakpoint window, in bp, for matching insertions."
        gq_cutoffs_tsv: "Table of genotype-quality cutoffs per variant class and length bin."
    }

    input {
        Array[File] vcfs
        Array[File] vcf_idxs
        Array[File] truth_vcfs
        Array[File] truth_vcf_idxs
        String prefix

        String? subset_vcf_string
        File ped
        File? swap_samples_truth
        Boolean skip_trv = true

        Array[Int] length_bins = [0, 1, 2, 6, 10, 30, 50, 100, 500, 5000, 50000]
        Int min_length_heuristic_comparison = 30

        Float del_size_similarity = 0.8
        Float del_reciprocal_overlap = 0.8
        Int del_breakpoint_window = 500
        Float ins_size_similarity = 0.8
        Int ins_breakpoint_window = 100

        String utils_docker

        RuntimeAttr? runtime_attr_find_trio_children
        RuntimeAttr? runtime_attr_swap_sample_ids
        RuntimeAttr? runtime_attr_subset_vcf
        RuntimeAttr? runtime_attr_subset_to_children
        RuntimeAttr? runtime_attr_cutoff_analysis
        RuntimeAttr? runtime_attr_merge
    }

    call FindTrioChildren {
        input:
            vcf = vcfs[0],
            vcf_idx = vcf_idxs[0],
            ped = ped,
            prefix = prefix,
            docker = utils_docker,
            runtime_attr_override = runtime_attr_find_trio_children
    }

    scatter (i in range(length(vcfs))) {
        if (defined(swap_samples_truth)) {
            call Helpers.SwapSampleIds as SwapTruthSampleIds {
                input:
                    vcf = truth_vcfs[i],
                    vcf_idx = truth_vcf_idxs[i],
                    sample_swap_list = select_first([swap_samples_truth]),
                    prefix = "~{prefix}.truth_swapped.~{i}",
                    docker = utils_docker,
                    runtime_attr_override = runtime_attr_swap_sample_ids
            }
        }

        if (defined(subset_vcf_string)) {
            call Helpers.SubsetVcfByArgs as SubsetVcf {
                input:
                    vcf = vcfs[i],
                    vcf_idx = vcf_idxs[i],
                    extra_args = select_first([subset_vcf_string]),
                    prefix = "~{prefix}.presolved.~{i}",
                    docker = utils_docker,
                    runtime_attr_override = runtime_attr_subset_vcf
            }
        }

        call Helpers.SubsetVcfToSamples as SubsetToTrioSamples {
            input:
                vcf = select_first([SubsetVcf.subset_vcf, vcfs[i]]),
                vcf_idx = select_first([SubsetVcf.subset_vcf_idx, vcf_idxs[i]]),
                samples = read_lines(FindTrioChildren.trio_sample_ids_file),
                prefix = "~{prefix}.trio_subset.~{i}",
                docker = utils_docker,
                runtime_attr_override = runtime_attr_subset_to_children
        }

        call Helpers.SubsetVcfToSamples as SubsetTruthToTrioSamples {
            input:
                vcf = select_first([SwapTruthSampleIds.swapped_vcf, truth_vcfs[i]]),
                vcf_idx = select_first([SwapTruthSampleIds.swapped_vcf_idx, truth_vcf_idxs[i]]),
                samples = read_lines(FindTrioChildren.trio_sample_ids_file),
                extra_args = "--force-samples",
                prefix = "~{prefix}.truth_trio_subset.~{i}",
                docker = utils_docker,
                runtime_attr_override = runtime_attr_subset_to_children
        }

        call CutoffAnalysis {
            input:
                vcf = SubsetToTrioSamples.subset_vcf,
                vcf_idx = SubsetToTrioSamples.subset_vcf_idx,
                truth_vcf = SubsetTruthToTrioSamples.subset_vcf,
                truth_vcf_idx = SubsetTruthToTrioSamples.subset_vcf_idx,
                trio_definitions = FindTrioChildren.trio_definitions,
                length_bins = length_bins,
                skip_trv = skip_trv,
                del_size_similarity = del_size_similarity,
                del_reciprocal_overlap = del_reciprocal_overlap,
                del_breakpoint_window = del_breakpoint_window,
                ins_size_similarity = ins_size_similarity,
                ins_breakpoint_window = ins_breakpoint_window,
                min_length_heuristic_comparison = min_length_heuristic_comparison,
                prefix = "~{prefix}.gq_cutoffs.~{i}",
                docker = utils_docker,
                runtime_attr_override = runtime_attr_cutoff_analysis
        }
    }

    call MergeResults {
        input:
            tsvs = CutoffAnalysis.cutoff_tsv,
            length_bins = length_bins,
            prefix = "~{prefix}.gq_cutoffs",
            docker = utils_docker,
            runtime_attr_override = runtime_attr_merge
    }

    output {
        File gq_cutoffs_tsv = MergeResults.merged_tsv
    }
}

task FindTrioChildren {
    input {
        File vcf
        File vcf_idx
        File ped
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        python3 <<CODE
import pysam

vcf = pysam.VariantFile("~{vcf}")
vcf_samples = set(vcf.header.samples)
vcf.close()

trios = []
with open("~{ped}") as f:
    for line in f:
        if line.startswith("#"):
            continue
        fields = line.strip().split("\t")
        if len(fields) < 4:
            continue
        sample, father, mother = fields[1], fields[2], fields[3]
        if father != "0" and mother != "0":
            if sample in vcf_samples and father in vcf_samples and mother in vcf_samples:
                trios.append((sample, father, mother))

with open("~{prefix}.trio_definitions.tsv", "w") as out:
    for child, father, mother in trios:
        out.write(f"{child}\t{father}\t{mother}\n")

all_samples = set()
for child, father, mother in trios:
    all_samples.update([child, father, mother])

with open("~{prefix}.trio_sample_ids.txt", "w") as out:
    for sample in sorted(all_samples):
        out.write(sample + "\n")
CODE
    >>>

    output {
        File trio_definitions = "~{prefix}.trio_definitions.tsv"
        File trio_sample_ids_file = "~{prefix}.trio_sample_ids.txt"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 4,
        disk_gb: 10,
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

task CutoffAnalysis {
    input {
        File vcf
        File vcf_idx
        File truth_vcf
        File truth_vcf_idx
        File trio_definitions
        Array[Int] length_bins
        Boolean skip_trv = true
        Float del_size_similarity
        Float del_reciprocal_overlap
        Int del_breakpoint_window
        Float ins_size_similarity
        Int ins_breakpoint_window
        Int min_length_heuristic_comparison
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        bcftools norm -m-any --no-version -Oz -o input.split.vcf.gz ~{vcf}
        tabix -p vcf input.split.vcf.gz
        bcftools norm -m-any --no-version -Oz -o truth.split.vcf.gz ~{truth_vcf}
        tabix -p vcf truth.split.vcf.gz

        python3 <<CODE
import pysam
from bisect import bisect_left, bisect_right
from collections import defaultdict

LENGTH_BINS = [~{sep=", " length_bins}]
if not LENGTH_BINS:
    raise ValueError("length_bins must not be empty")
if any(left >= right for left, right in zip(LENGTH_BINS, LENGTH_BINS[1:])):
    raise ValueError("length_bins must be strictly increasing")

SIZE_LABELS = [f"{start}-{end - 1}" for start, end in zip(LENGTH_BINS, LENGTH_BINS[1:])] + [f"{LENGTH_BINS[-1]}+"]
SIZE_ORDER = {label: index for index, label in enumerate(SIZE_LABELS)}

DEL_SIZE_SIM = float(~{del_size_similarity})
DEL_REC_OVL = float(~{del_reciprocal_overlap})
DEL_BP_WIN = int(~{del_breakpoint_window})
INS_SIZE_SIM = float(~{ins_size_similarity})
INS_BP_WIN = int(~{ins_breakpoint_window})
MIN_LEN_HEUR = int(~{min_length_heuristic_comparison})
SKIP_TRV = ~{if skip_trv then "True" else "False"}

def get_type(variant_id):
    vid = (variant_id or "").upper()
    if "INS" in vid:
        return "INS"
    elif "DEL" in vid:
        return "DEL"
    elif "TRV" in vid:
        return "TRV"
    return "SNV"

def get_size_bucket(allele_length):
    size = abs(allele_length)
    for index in range(len(LENGTH_BINS)):
        if index + 1 == len(LENGTH_BINS) or size < LENGTH_BINS[index + 1]:
            return SIZE_LABELS[index]

def is_nonref(gt):
    return gt is not None and any(a is not None and a != 0 for a in gt)

def is_ref(gt):
    return gt is not None and all(a == 0 for a in gt)

def get_info_int(rec, field):
    raw = rec.info.get(field)
    if raw is None:
        return None
    if isinstance(raw, (list, tuple)):
        raw = raw[0]
    try:
        return int(raw)
    except (ValueError, TypeError):
        return None

def passes_del(qs, qe, ql, cs, ce, cl):
    if ql == 0 or cl == 0:
        return False
    if min(ql, cl) / max(ql, cl) < DEL_SIZE_SIM:
        return False
    overlap = max(0, min(qe, ce) - max(qs, cs))
    if overlap / min(ql, cl) < DEL_REC_OVL:
        return False
    if abs(qs - cs) > DEL_BP_WIN and abs(qe - ce) > DEL_BP_WIN:
        return False
    return True

def passes_ins(qs, ql, cs, cl):
    if ql == 0 or cl == 0:
        return False
    if min(ql, cl) / max(ql, cl) < INS_SIZE_SIM:
        return False
    if abs(qs - cs) > INS_BP_WIN:
        return False
    return True

trios = []
with open("~{trio_definitions}") as f:
    for line in f:
        child, father, mother = line.strip().split("\t")
        trios.append((child, father, mother))
children = [c for c, _, _ in trios]
parents_of = {c: (f, m) for c, f, m in trios}

# Truth side
truth_in = pysam.VariantFile("truth.split.vcf.gz")
truth_samples = set(truth_in.header.samples)
common_children = [c for c in children if c in truth_samples]
common_children_set = set(common_children)

exact_lookup = {}                                  # (chrom, pos, ref, alt) -> frozenset(nonref_children)
sv_index = defaultdict(lambda: defaultdict(list))  # type -> chrom -> [(start, stop, abslen, nonref_children)]
truth_total_counts = defaultdict(int)              # (bucket_type, bucket_size) -> total truth child-variants

for record in truth_in:
    svtype = record.info.get("SVTYPE")
    if isinstance(svtype, (list, tuple)):
        svtype = svtype[0]
    svlen = get_info_int(record, "SVLEN") or 0
    if svtype == "SNV":
        svlen = 0
    nonref = []
    for s in common_children:
        if is_nonref(record.samples[s].get("GT")):
            nonref.append(s)
    if not nonref and svtype not in ("DEL", "INS"):
        # still need exact lookup empty for precision, but skip storing empty sets
        pass

    abslen = abs(svlen)
    if svtype in ("DEL", "INS") and abslen >= MIN_LEN_HEUR:
        sv_index[svtype][record.chrom].append((record.start, record.stop, abslen, frozenset(nonref)))
        size_bucket = get_size_bucket(svlen)
        truth_total_counts[(svtype, size_bucket)] += len(nonref)
    else:
        for alt in (record.alts or []):
            key = (record.chrom, record.pos, record.ref, alt)
            prev = exact_lookup.get(key)
            new_set = frozenset(nonref)
            if prev is None:
                exact_lookup[key] = new_set
            else:
                exact_lookup[key] = prev | new_set
        # Bucket assignment for truth-total: type from SVTYPE if available else SNV
        truth_type = svtype if svtype in ("DEL", "INS", "SNV") else "SNV"
        size_bucket = get_size_bucket(svlen)
        truth_total_counts[(truth_type, size_bucket)] += len(nonref)

truth_in.close()

# Sort SV index entries by start for binary search
for t in sv_index:
    for chrom in sv_index[t]:
        sv_index[t][chrom].sort(key=lambda x: x[0])

def heuristic_match(svtype, chrom, qs, qe, ql):
    if chrom not in sv_index[svtype]:
        return frozenset()
    entries = sv_index[svtype][chrom]
    if svtype == "INS":
        win = INS_BP_WIN
        lo = qs - win
        hi = qs + win + 1
    else:
        margin = DEL_BP_WIN + int(ql / DEL_SIZE_SIM) + 1
        lo = qs - margin
        hi = qe + margin
    # Binary search by start; widen with a margin to catch overlapping intervals
    starts = [e[0] for e in entries]
    left = bisect_left(starts, lo - 1)
    # Linear scan from left until start > hi
    matched = set()
    for i in range(left, len(entries)):
        cs, ce, cl, nonref = entries[i]
        if cs > hi:
            break
        if svtype == "INS":
            ok = passes_ins(qs, ql, cs, cl)
        else:
            ok = passes_del(qs, qe, ql, cs, ce, cl)
        if ok:
            matched |= nonref
    return frozenset(matched)

# Per (type, size, gq) -> [n_calls, n_truth, n_denovo]
counts = defaultdict(lambda: [0, 0, 0])

vcf_in = pysam.VariantFile("input.split.vcf.gz")
for record in vcf_in:
    vtype = get_type(record.id)
    if SKIP_TRV and vtype == "TRV":
        continue
    svlen = get_info_int(record, "SVLEN")
    if svlen is None:
        svlen = get_info_int(record, "allele_length") or 0
    abslen = abs(svlen)
    size_bucket = get_size_bucket(svlen)

    # Determine truth-supported child set
    if vtype in ("DEL", "INS") and abslen >= MIN_LEN_HEUR:
        truth_children = heuristic_match(vtype, record.chrom, record.start, record.stop, abslen)
    else:
        truth_children = frozenset()
        for alt in (record.alts or []):
            key = (record.chrom, record.pos, record.ref, alt)
            entry = exact_lookup.get(key)
            if entry:
                truth_children = truth_children | entry

    for child, father, mother in trios:
        child_gt = record.samples[child].get("GT")
        if not is_nonref(child_gt):
            continue
        gq = record.samples[child].get("GQ")
        if gq is None:
            continue
        key = (vtype, size_bucket, int(gq))
        counts[key][0] += 1
        if child in truth_children:
            counts[key][1] += 1
        f_gt = record.samples[father].get("GT")
        m_gt = record.samples[mother].get("GT")
        if is_ref(f_gt) and is_ref(m_gt):
            counts[key][2] += 1

vcf_in.close()

n_children = len(children)

with open("~{prefix}.tsv", "w") as out:
    out.write("BUCKET_TYPE\tBUCKET_SIZE\tGQ\tN_CALLS\tN_CALLS_TRUTH_SUPPORTED\tN_DENOVO\tN_CHILDREN\tN_TRUTH_TOTAL\n")
    for key in sorted(counts.keys(), key=lambda k: (k[0], SIZE_ORDER.get(k[1], 99), k[2])):
        bt, bs, gq = key
        ncalls, ntruth, ndenovo = counts[key]
        ntruth_total = truth_total_counts.get((bt, bs), 0)
        out.write(f"{bt}\t{bs}\t{gq}\t{ncalls}\t{ntruth}\t{ndenovo}\t{n_children}\t{ntruth_total}\n")
CODE
    >>>

    output {
        File cutoff_tsv = "~{prefix}.tsv"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 16,
        disk_gb: 2 * ceil(size(vcf, "GB")) + 2 * ceil(size(truth_vcf, "GB")) + 5,
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

task MergeResults {
    input {
        Array[File] tsvs
        Array[Int] length_bins
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        python3 <<CODE
from collections import defaultdict

LENGTH_BINS = [~{sep=", " length_bins}]
SIZE_LABELS = [f"{start}-{end - 1}" for start, end in zip(LENGTH_BINS, LENGTH_BINS[1:])] + [f"{LENGTH_BINS[-1]}+"]
SIZE_ORDER = {label: index for index, label in enumerate(SIZE_LABELS)}

# (bucket_type, bucket_size, gq) -> [n_calls, n_truth, n_denovo]
counts = defaultdict(lambda: [0, 0, 0])
# (bucket_type, bucket_size) -> total truth (summed across shards)
truth_totals = defaultdict(int)
truth_totals_seen = defaultdict(set)
n_children = 0

input_files = "~{sep=',' tsvs}".split(",")
for f in input_files:
    with open(f) as fh:
        next(fh)
        for line in fh:
            fields = line.rstrip("\n").split("\t")
            key = (fields[0], fields[1], int(fields[2]))
            counts[key][0] += int(fields[3])
            counts[key][1] += int(fields[4])
            counts[key][2] += int(fields[5])
            n_children = max(n_children, int(fields[6]))
            tt_key = (fields[0], fields[1])
            if f not in truth_totals_seen[tt_key]:
                truth_totals[tt_key] += int(fields[7])
                truth_totals_seen[tt_key].add(f)

with open("~{prefix}.tsv", "w") as out:
    out.write("BUCKET_TYPE\tBUCKET_SIZE\tGQ\tN_CALLS\tN_CALLS_TRUTH_SUPPORTED\tN_DENOVO\tN_CHILDREN\tN_TRUTH_TOTAL\n")
    for key in sorted(counts.keys(), key=lambda k: (k[0], SIZE_ORDER.get(k[1], 99), k[2])):
        bt, bs, gq = key
        ncalls, ntruth, ndenovo = counts[key]
        ntruth_total = truth_totals.get((bt, bs), 0)
        out.write(f"{bt}\t{bs}\t{gq}\t{ncalls}\t{ntruth}\t{ndenovo}\t{n_children}\t{ntruth_total}\n")
CODE
    >>>

    output {
        File merged_tsv = "~{prefix}.tsv"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 4,
        disk_gb: 10,
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
