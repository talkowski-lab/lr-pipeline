version 1.0

import "../utils/Helpers.wdl"
import "../utils/Structs.wdl"

workflow GQCalculateCounts {
    meta {
        description: [
            "This utility computes GQ-stratified count tables used to derive GQ filtering cutoffs, from both a trio de novo analysis and a truth-set concordance analysis. Counts are bucketed by variant type, allele-length bin and supporting caller. For structural variants (`abs(allele_length) >= 50`), the `CALLER` column expands each call by its supporting callers using the `EV`/`BEV` FORMAT fields written by `AnnotateSvCallerSupport`: `kanpig`-backed calls are recorded under `CALLER=kanpig` with their own GQ, calls backed by other callers are split into one row per caller carrying an allelic depth in `EV` (with a per-caller GQ recomputed from that depth), and calls with no `BEV` are recorded with a blank `CALLER`. It outputs one TSV per analysis."
        ]
    }

    parameter_meta {
        vcfs: "Cohort VCFs to analyze."
        vcf_idxs: "Indexes for the cohort VCFs."
        truth_vcfs: "Truth-set VCFs, one per input VCF, for the concordance analysis."
        truth_vcf_idxs: "Indexes for the truth-set VCFs."
        length_bins: "Allele-length bin boundaries defining the size buckets."
        subset_vcf_string: "Optional `bcftools view` argument string to pre-subset each VCF."
        ped: "Pedigree used to identify trios for the de novo analysis."
        swap_samples_truth: "Optional sample-swap list applied to the truth VCFs."
        run_trio_qc: "Whether to run the trio de novo analysis."
        run_truth_qc: "Whether to run the truth-set concordance analysis."
        skip_trv: "Whether to skip tandem-repeat variants."
        drop_kanpig_supported_gq: "Whether a Kanpig-supported call also expands its `EV` callers, so co-supporting callers contribute their own genotype-quality rows."
        min_fuzzy_match: "Minimum variant length to perform fuzzy matching for truth concordance."
        del_breakpoint_window: "Breakpoint window, in bp, for matching deletions during truth concordance."
        del_reciprocal_overlap: "Minimum reciprocal overlap for matching deletions during truth concordance."
        del_size_similarity: "Minimum size similarity for matching deletions during truth concordance."
        ins_breakpoint_window: "Breakpoint window, in bp, for matching insertions during truth concordance."
        ins_reciprocal_overlap: "Minimum reciprocal overlap for matching insertions during truth concordance."
        ins_size_similarity: "Minimum size similarity for matching insertions during truth concordance."
        trio_denovo_tsv: "GQ-stratified trio de novo count table."
        truth_concordance_tsv: "GQ-stratified truth-set concordance count table."
    }

    input {
        Array[File] vcfs
        Array[File] vcf_idxs
        Array[File]? truth_vcfs
        Array[File]? truth_vcf_idxs
        String prefix

        Array[Int] length_bins = [0, 1, 10, 30, 50, 100, 500, 5000, 50000]
        String? subset_vcf_string
        File? ped
        File? swap_samples_truth
        Boolean run_trio_qc = true
        Boolean run_truth_qc = true
        Boolean skip_trv = true
        Boolean drop_kanpig_supported_gq = false

        Int min_fuzzy_match = 20

        Int del_breakpoint_window = 500
        Float del_reciprocal_overlap = 0.7
        Float del_size_similarity = 0.7

        Int ins_breakpoint_window = 200
        Float ins_reciprocal_overlap = 0.0
        Float ins_size_similarity = 0.5

        String utils_docker

        RuntimeAttr? runtime_attr_find_trios
        RuntimeAttr? runtime_attr_swap_sample_ids
        RuntimeAttr? runtime_attr_subset_vcf
        RuntimeAttr? runtime_attr_trio_analysis
        RuntimeAttr? runtime_attr_truth_analysis
        RuntimeAttr? runtime_attr_merge_trio
        RuntimeAttr? runtime_attr_merge_truth
    }

    if (run_trio_qc) {
        # Trio de novo analysis
        call FindTrios {
            input:
                vcf = vcfs[0],
                vcf_idx = vcf_idxs[0],
                ped = select_first([ped]),
                prefix = prefix,
                docker = utils_docker,
                runtime_attr_override = runtime_attr_find_trios
        }

        scatter (i in range(length(vcfs))) {
            if (defined(subset_vcf_string)) {
                call Helpers.SubsetVcfByArgs as SubsetTrioVcf {
                    input:
                        vcf = vcfs[i],
                        vcf_idx = vcf_idxs[i],
                        extra_args = select_first([subset_vcf_string]),
                        prefix = "~{prefix}.trio_presolved.~{i}",
                        docker = utils_docker,
                        runtime_attr_override = runtime_attr_subset_vcf
                }
            }

            call Helpers.SubsetVcfToSamples as SubsetToTrioSamples {
                input:
                    vcf = select_first([SubsetTrioVcf.subset_vcf, vcfs[i]]),
                    vcf_idx = select_first([SubsetTrioVcf.subset_vcf_idx, vcf_idxs[i]]),
                    samples = read_lines(FindTrios.trio_sample_ids_file),
                    prefix = "~{prefix}.trio_subset.~{i}",
                    docker = utils_docker,
                    runtime_attr_override = runtime_attr_subset_vcf
            }

            call TrioDeNovoAnalysis {
                input:
                    vcf = SubsetToTrioSamples.subset_vcf,
                    vcf_idx = SubsetToTrioSamples.subset_vcf_idx,
                    trio_definitions = FindTrios.trio_definitions,
                    length_bins = length_bins,
                    skip_trv = skip_trv,
                    drop_kanpig_supported_gq = drop_kanpig_supported_gq,
                    prefix = "~{prefix}.trio_denovo.~{i}",
                    docker = utils_docker,
                    runtime_attr_override = runtime_attr_trio_analysis
            }
        }

        call MergeTrioResults {
            input:
                tsvs = TrioDeNovoAnalysis.trio_denovo_tsv,
                length_bins = length_bins,
                prefix = "~{prefix}.trio_denovo",
                docker = utils_docker,
                runtime_attr_override = runtime_attr_merge_trio
        }
    }

    if (run_truth_qc) {
        # Truth set concordance analysis
        scatter (i in range(length(vcfs))) {
            if (defined(swap_samples_truth)) {
                call Helpers.SwapSampleIds as SwapTruthSampleIds {
                    input:
                        vcf = select_first([truth_vcfs])[i],
                        vcf_idx = select_first([truth_vcf_idxs])[i],
                        sample_swap_list = select_first([swap_samples_truth]),
                        prefix = "~{prefix}.truth_swapped.~{i}",
                        docker = utils_docker,
                        runtime_attr_override = runtime_attr_swap_sample_ids
                }
            }

            if (defined(subset_vcf_string)) {
                call Helpers.SubsetVcfByArgs as SubsetTruthEvalVcf {
                    input:
                        vcf = vcfs[i],
                        vcf_idx = vcf_idxs[i],
                        extra_args = select_first([subset_vcf_string]),
                        prefix = "~{prefix}.truth_presolved.~{i}",
                        docker = utils_docker,
                        runtime_attr_override = runtime_attr_subset_vcf
                }
            }

            call TruthSetAnalysis {
                input:
                    vcf = select_first([SubsetTruthEvalVcf.subset_vcf, vcfs[i]]),
                    vcf_idx = select_first([SubsetTruthEvalVcf.subset_vcf_idx, vcf_idxs[i]]),
                    truth_vcf = select_first([SwapTruthSampleIds.swapped_vcf, select_first([truth_vcfs])[i]]),
                    truth_vcf_idx = select_first([SwapTruthSampleIds.swapped_vcf_idx, select_first([truth_vcf_idxs])[i]]),
                    length_bins = length_bins,
                    skip_trv = skip_trv,
                    drop_kanpig_supported_gq = drop_kanpig_supported_gq,
                    min_fuzzy_match = min_fuzzy_match,
                    del_breakpoint_window = del_breakpoint_window,
                    del_reciprocal_overlap = del_reciprocal_overlap,
                    del_size_similarity = del_size_similarity,
                    ins_breakpoint_window = ins_breakpoint_window,
                    ins_reciprocal_overlap = ins_reciprocal_overlap,
                    ins_size_similarity = ins_size_similarity,
                    prefix = "~{prefix}.truth.~{i}",
                    docker = utils_docker,
                    runtime_attr_override = runtime_attr_truth_analysis
            }
        }

        call MergeTruthResults {
            input:
                tsvs = TruthSetAnalysis.truth_concordance_tsv,
                length_bins = length_bins,
                prefix = "~{prefix}.truth_concordance",
                docker = utils_docker,
                runtime_attr_override = runtime_attr_merge_truth
        }
    }

    output {
        File? trio_denovo_tsv = MergeTrioResults.merged_tsv
        File? truth_concordance_tsv = MergeTruthResults.merged_tsv
    }
}

task FindTrios {
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

task TrioDeNovoAnalysis {
    input {
        File vcf
        File vcf_idx
        File trio_definitions
        Array[Int] length_bins
        Boolean skip_trv = true
        Boolean drop_kanpig_supported_gq = false
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        python3 <<CODE
import math
import re
import pysam
from collections import defaultdict

drop_kanpig_supported_gq = ~{true="True" false="False" drop_kanpig_supported_gq}

trios = []
with open("~{trio_definitions}") as f:
    for line in f:
        child, father, mother = line.strip().split("\t")
        trios.append((child, father, mother))

LENGTH_BINS = [~{sep=", " length_bins}]

if not LENGTH_BINS:
    raise ValueError("length_bins must not be empty")
if any(left >= right for left, right in zip(LENGTH_BINS, LENGTH_BINS[1:])):
    raise ValueError("length_bins must be strictly increasing")

SIZE_LABELS = [f"{start}-{end - 1}" for start, end in zip(LENGTH_BINS, LENGTH_BINS[1:])] + [f"{LENGTH_BINS[-1]}+"]
SIZE_ORDER = {label: index for index, label in enumerate(SIZE_LABELS)}

EV_RE = re.compile(r"^(.+)_\((\d+)_(\d+)\)$")

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
    for index, start in enumerate(LENGTH_BINS):
        if index + 1 == len(LENGTH_BINS) or size < LENGTH_BINS[index + 1]:
            return SIZE_LABELS[index]

def is_nonref_unphased(gt):
    return gt is not None and any(a is not None and a != 0 for a in gt)

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

def caller_gq_emissions(sample):
    # (caller, gq) rows for one non-ref call: no_caller if no BEV, kanpig if BEV=kanpig,
    # otherwise one row per EV caller carrying an AD, with a per-caller GQ from that AD.
    # When drop_kanpig_supported_gq, kanpig-BEV calls also expand their EV so co-supporting
    # callers feed the non-kanpig plots (kanpig itself stays via its EV entry).
    bev = sample.get("BEV")
    if bev is None:
        gq = sample.get("GQ")
        return [] if gq is None else [("no_caller", gq)]
    if bev == "kanpig" and not drop_kanpig_supported_gq:
        gq = sample.get("GQ")
        return [] if gq is None else [("kanpig", gq)]
    emissions = []
    for entry in (sample.get("EV") or ()):
        match = EV_RE.match(entry)
        if not match:
            continue
        emissions.append((match.group(1), calculate_gq(calculate_pl(int(match.group(2)), int(match.group(3))))))
    return emissions

# (bucket_type, bucket_size, caller, gq, trv_status) -> [count_inherited, count_de_novo]
counts = defaultdict(lambda: [0, 0])

vcf = pysam.VariantFile("~{vcf}")
for record in vcf:
    if ~{if skip_trv then "True" else "False"} and "TRV" in (record.id or "").upper():
        continue

    al = record.info.get("allele_length")
    if al is None:
        al = 0
    elif isinstance(al, (list, tuple)):
        al = al[0]

    variant_type = get_type(record.id)
    size_bucket = get_size_bucket(al)
    trv_status = "TR_ENVELOPED" in record.info

    for child, father, mother in trios:
        child_sample = record.samples[child]
        if not is_nonref_unphased(child_sample["GT"]):
            continue

        emissions = caller_gq_emissions(child_sample)
        if not emissions:
            continue

        inherited = is_nonref_unphased(record.samples[father]["GT"]) or is_nonref_unphased(record.samples[mother]["GT"])
        slot = 0 if inherited else 1
        for caller, gq in emissions:
            counts[(variant_type, size_bucket, caller, gq, trv_status)][slot] += 1

vcf.close()

with open("~{prefix}.tsv", "w") as out:
    out.write("BUCKET_TYPE\tBUCKET_SIZE\tCALLER\tGQ\tTRV_STATUS\tCOUNT\tCOUNT_INHERITED\tCOUNT_DE_NOVO\n")
    for key in sorted(counts.keys(), key=lambda k: (k[0], SIZE_ORDER.get(k[1], 99), k[2], k[3], k[4])):
        bt, bs, caller, gq, trv_status = key
        inherited, de_novo = counts[key]
        count = inherited + de_novo
        out.write(f"{bt}\t{bs}\t{caller}\t{gq}\t{trv_status}\t{count}\t{inherited}\t{de_novo}\n")
CODE
    >>>

    output {
        File trio_denovo_tsv = "~{prefix}.tsv"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 8,
        disk_gb: 2 * ceil(size(vcf, "GB")) + 5,
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

task MergeTrioResults {
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

counts = defaultdict(lambda: [0, 0])
input_files = "~{sep=',' tsvs}".split(",")

for f in input_files:
    with open(f) as fh:
        next(fh)
        for line in fh:
            fields = line.rstrip("\n").split("\t")
            key = (fields[0], fields[1], fields[2], int(fields[3]), fields[4] == "True")
            counts[key][0] += int(fields[6])
            counts[key][1] += int(fields[7])

with open("~{prefix}.tsv", "w") as out:
    out.write("BUCKET_TYPE\tBUCKET_SIZE\tCALLER\tGQ\tTRV_STATUS\tCOUNT\tCOUNT_INHERITED\tCOUNT_DE_NOVO\n")
    for key in sorted(counts.keys(), key=lambda k: (k[0], SIZE_ORDER.get(k[1], 99), k[2], k[3], k[4])):
        bt, bs, caller, gq, trv_status = key
        inherited, de_novo = counts[key]
        count = inherited + de_novo
        out.write(f"{bt}\t{bs}\t{caller}\t{gq}\t{trv_status}\t{count}\t{inherited}\t{de_novo}\n")
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

task TruthSetAnalysis {
    input {
        File vcf
        File vcf_idx
        File truth_vcf
        File truth_vcf_idx
        Array[Int] length_bins
        Boolean skip_trv = true
        Boolean drop_kanpig_supported_gq = false
        Int min_fuzzy_match
        Int del_breakpoint_window
        Float del_reciprocal_overlap
        Float del_size_similarity
        Int ins_breakpoint_window
        Float ins_reciprocal_overlap
        Float ins_size_similarity
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        # Subset base VCF to common samples
        bcftools query -l ~{vcf} | sort > vcf_samples.txt
        bcftools query -l ~{truth_vcf} | sort > truth_samples.txt
        comm -12 vcf_samples.txt truth_samples.txt > common_samples.txt

        bcftools view -S common_samples.txt --min-ac 1 -Oz -o subset.vcf.gz ~{vcf}
        tabix -p vcf subset.vcf.gz

        if [[ "~{truth_vcf_idx}" != "~{truth_vcf}.tbi" ]]; then
            ln -sf "~{truth_vcf_idx}" "~{truth_vcf}.tbi"
        fi

        python3 <<CODE
import math
import re
import pysam
from collections import defaultdict

LENGTH_BINS = [~{sep=", " length_bins}]

if not LENGTH_BINS:
    raise ValueError("length_bins must not be empty")
if any(left >= right for left, right in zip(LENGTH_BINS, LENGTH_BINS[1:])):
    raise ValueError("length_bins must be strictly increasing")

SIZE_LABELS = [f"{start}-{end - 1}" for start, end in zip(LENGTH_BINS, LENGTH_BINS[1:])] + [f"{LENGTH_BINS[-1]}+"]
SIZE_ORDER = {label: index for index, label in enumerate(SIZE_LABELS)}

EV_RE = re.compile(r"^(.+)_\((\d+)_(\d+)\)$")

drop_kanpig_supported_gq = ~{true="True" false="False" drop_kanpig_supported_gq}
MIN_FUZZY_MATCH = int(~{min_fuzzy_match})
DEL_SIZE_SIM = float(~{del_size_similarity})
DEL_REC_OVL = float(~{del_reciprocal_overlap})
DEL_BP_WIN = int(~{del_breakpoint_window})
INS_SIZE_SIM = float(~{ins_size_similarity})
INS_REC_OVL = float(~{ins_reciprocal_overlap})
INS_BP_WIN = int(~{ins_breakpoint_window})

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
    for index, start in enumerate(LENGTH_BINS):
        if index + 1 == len(LENGTH_BINS) or size < LENGTH_BINS[index + 1]:
            return SIZE_LABELS[index]

def is_nonref_unphased(gt):
    return gt is not None and any(a is not None and a != 0 for a in gt)

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

def caller_gq_emissions(sample):
    # (caller, gq) rows for one non-ref call: no_caller if no BEV, kanpig if BEV=kanpig,
    # otherwise one row per EV caller carrying an AD, with a per-caller GQ from that AD.
    # When drop_kanpig_supported_gq, kanpig-BEV calls also expand their EV so co-supporting
    # callers feed the non-kanpig plots (kanpig itself stays via its EV entry).
    bev = sample.get("BEV")
    if bev is None:
        gq = sample.get("GQ")
        return [] if gq is None else [("no_caller", gq)]
    if bev == "kanpig" and not drop_kanpig_supported_gq:
        gq = sample.get("GQ")
        return [] if gq is None else [("kanpig", gq)]
    emissions = []
    for entry in (sample.get("EV") or ()):
        match = EV_RE.match(entry)
        if not match:
            continue
        emissions.append((match.group(1), calculate_gq(calculate_pl(int(match.group(2)), int(match.group(3))))))
    return emissions

# DEL/INS use fuzzy size/overlap/breakpoint matching; SNVs match exactly on CHROM+POS+REF+ALT.
def passes_del(qs, qe, ql, cs, ce, cl):
    if ql == 0 or cl == 0:
        return False
    if min(ql, cl) / max(ql, cl) < DEL_SIZE_SIM:
        return False
    overlap = max(0, min(qe, ce) - max(qs, cs))
    if overlap / min(ql, cl) < DEL_REC_OVL:
        return False
    if abs(qs - cs) + abs(qe - ce) > DEL_BP_WIN:
        return False
    return True

def passes_ins(qs, qe, ql, cs, ce, cl):
    if ql == 0 or cl == 0:
        return False
    if min(ql, cl) / max(ql, cl) < INS_SIZE_SIM:
        return False
    overlap = max(0, min(qe, ce) - max(qs, cs))
    if overlap / min(ql, cl) < INS_REC_OVL:
        return False
    if abs(qs - cs) + abs(qe - ce) > INS_BP_WIN:
        return False
    return True

with open("common_samples.txt") as f:
    common_samples = [line.strip() for line in f if line.strip()]

truth_in = pysam.VariantFile("~{truth_vcf}")
truth_samples = set(truth_in.header.samples)
truth_contigs = set(truth_in.header.contigs)
common_in_truth = [s for s in common_samples if s in truth_samples]

# Non-ref sample set per truth record, cached so the same candidate is scanned once.
nonref_cache = {}
def truth_nonref_set(cand):
    cache_key = (cand.chrom, cand.pos, cand.ref, cand.alts, cand.id)
    cached = nonref_cache.get(cache_key)
    if cached is not None:
        return cached
    nonref = set()
    for sample in common_in_truth:
        if is_nonref_unphased(cand.samples[sample]["GT"]):
            nonref.add(sample)
    nonref_cache[cache_key] = nonref
    return nonref

def find_truth_match(record, variant_type, ql):
    chrom = record.chrom
    if chrom not in truth_contigs:
        return False, set()
    has_match = False
    truth_nonref = set()

    # Fuzzy matching only for DEL/INS at or above the size threshold; everything else (SNVs and
    # sub-threshold indels) is matched exactly on CHROM+POS+REF+ALT.
    if not (variant_type in ("DEL", "INS") and ql >= MIN_FUZZY_MATCH):
        ref_upper = record.ref.upper() if record.ref else record.ref
        alts_upper = set(a.upper() for a in (record.alts or ()))
        for cand in truth_in.fetch(chrom, record.start, record.start + 1):
            if cand.pos != record.pos:
                continue
            cand_ref = cand.ref.upper() if cand.ref else cand.ref
            if cand_ref != ref_upper:
                continue
            if alts_upper & set(a.upper() for a in (cand.alts or ())):
                has_match = True
                truth_nonref |= truth_nonref_set(cand)
        return has_match, truth_nonref

    qs = record.start
    qe = record.stop
    if variant_type == "DEL":
        margin = DEL_BP_WIN + (int(ql / DEL_SIZE_SIM) if DEL_SIZE_SIM > 0 else ql) + 1
        region_start = max(0, qs - margin)
        region_end = qe + margin
    else:
        region_start = max(0, qs - INS_BP_WIN)
        region_end = qs + INS_BP_WIN + 1
    for cand in truth_in.fetch(chrom, region_start, region_end):
        cand_ref = cand.ref or ""
        cs = cand.start
        ce = cand.start + len(cand_ref)
        for alt in (cand.alts or ()):
            signed = len(alt) - len(cand_ref)
            cl = abs(signed)
            if variant_type == "DEL":
                if signed >= 0:
                    continue
                ok = passes_del(qs, qe, ql, cs, ce, cl)
            else:
                if signed <= 0:
                    continue
                ok = passes_ins(qs, qe, ql, cs, ce, cl)
            if ok:
                has_match = True
                truth_nonref |= truth_nonref_set(cand)
                break
    return has_match, truth_nonref

# Per (bucket_type, bucket_size, caller, gq) tracking
variant_sites = defaultdict(set)
variant_match_sites = defaultdict(set)
match_call_count = defaultdict(int)
match_concordant_count = defaultdict(int)

vcf_in = pysam.VariantFile("subset.vcf.gz")
for record in vcf_in:
    if ~{if skip_trv then "True" else "False"} and "TRV" in (record.id or "").upper():
        continue

    al = record.info.get("allele_length")
    if al is None:
        al = 0
    elif isinstance(al, (list, tuple)):
        al = al[0]

    variant_type = get_type(record.id)
    size_bucket = get_size_bucket(al)

    has_match, truth_nonref = find_truth_match(record, variant_type, abs(al))
    trv_status = "TR_ENVELOPED" in record.info

    for s in common_samples:
        sample = record.samples[s]
        if not is_nonref_unphased(sample["GT"]):
            continue

        emissions = caller_gq_emissions(sample)
        if not emissions:
            continue

        for caller, gq in emissions:
            bucket_key = (variant_type, size_bucket, caller, gq, trv_status)
            variant_sites[bucket_key].add(record.id)

            if has_match:
                variant_match_sites[bucket_key].add(record.id)
                match_call_count[bucket_key] += 1
                if s in truth_nonref:
                    match_concordant_count[bucket_key] += 1

vcf_in.close()
truth_in.close()

all_keys = sorted(variant_sites.keys(), key=lambda k: (k[0], SIZE_ORDER.get(k[1], 99), k[2], k[3], k[4]))
with open("~{prefix}.tsv", "w") as out:
    out.write("BUCKET_TYPE\tBUCKET_SIZE\tCALLER\tGQ\tTRV_STATUS\tCOUNT_VARIANT\tCOUNT_VARIANT_MATCH\tCOUNT_CALL_MATCH\tCOUNT_CALL_MATCH_CONCORDANT\tCOUNT_CALL_MATCH_DISCONCORDANT\n")
    for key in all_keys:
        bt, bs, caller, gq, trv_status = key
        vc = len(variant_sites[key])
        vmc = len(variant_match_sites.get(key, set()))
        mcc = match_call_count.get(key, 0)
        mccc = match_concordant_count.get(key, 0)
        out.write(f"{bt}\t{bs}\t{caller}\t{gq}\t{trv_status}\t{vc}\t{vmc}\t{mcc}\t{mccc}\t{mcc - mccc}\n")
CODE
    >>>

    output {
        File truth_concordance_tsv = "~{prefix}.tsv"
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
        memory: 20 + " GiB"
        disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
        bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
        docker: docker
        preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
    }
}

task MergeTruthResults {
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

# (bucket_type, bucket_size, caller, gq, trv_status) -> [variant_count, variant_match_count, match_call_count, match_concordant_count]
counts = defaultdict(lambda: [0, 0, 0, 0])
input_files = "~{sep=',' tsvs}".split(",")

for f in input_files:
    with open(f) as fh:
        next(fh)
        for line in fh:
            fields = line.rstrip("\n").split("\t")
            key = (fields[0], fields[1], fields[2], int(fields[3]), fields[4] == "True")
            counts[key][0] += int(fields[5])
            counts[key][1] += int(fields[6])
            counts[key][2] += int(fields[7])
            counts[key][3] += int(fields[8])

with open("~{prefix}.tsv", "w") as out:
    out.write("BUCKET_TYPE\tBUCKET_SIZE\tCALLER\tGQ\tTRV_STATUS\tCOUNT_VARIANT\tCOUNT_VARIANT_MATCH\tCOUNT_CALL_MATCH\tCOUNT_CALL_MATCH_CONCORDANT\tCOUNT_CALL_MATCH_DISCONCORDANT\n")
    for key in sorted(counts.keys(), key=lambda k: (k[0], SIZE_ORDER.get(k[1], 99), k[2], k[3], k[4])):
        bt, bs, caller, gq, trv_status = key
        vc, vmc, mcc, mccc = counts[key]
        out.write(f"{bt}\t{bs}\t{caller}\t{gq}\t{trv_status}\t{vc}\t{vmc}\t{mcc}\t{mccc}\t{mcc - mccc}\n")
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
