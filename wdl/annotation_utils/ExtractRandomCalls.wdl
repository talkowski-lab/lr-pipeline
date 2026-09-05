version 1.0

import "../utils/Structs.wdl"

workflow ExtractRandomCalls {
    input {
        Array[File] vcfs
        Array[File] vcf_idxs
        String prefix

        Int count
        Float? min_af
        Float? max_af
        Int? min_ac
        Int? max_ac
        Boolean singleton = false
        Array[String] filter_values = []
        Array[String] allele_types = []
        Int? min_allele_length
        Int? max_allele_length
        Array[String] variant_types = []
        Array[String] include_samples = []
        Array[String] exclude_samples = []
        Int random_seed = 42

        String utils_docker

        RuntimeAttr? runtime_attr_sample
        RuntimeAttr? runtime_attr_merge
    }

    scatter (i in range(length(vcfs))) {
        call SampleShardCalls {
            input:
                vcf = vcfs[i],
                vcf_idx = vcf_idxs[i],
                count = count,
                min_af = min_af,
                max_af = max_af,
                min_ac = min_ac,
                max_ac = max_ac,
                singleton = singleton,
                filter_values = filter_values,
                allele_types = allele_types,
                min_allele_length = min_allele_length,
                max_allele_length = max_allele_length,
                variant_types = variant_types,
                include_samples = include_samples,
                exclude_samples = exclude_samples,
                random_seed = random_seed + i,
                prefix = "~{prefix}.shard_~{i}",
                docker = utils_docker,
                runtime_attr_override = runtime_attr_sample
        }
    }

    call MergeRandomCalls {
        input:
            pair_tsvs = SampleShardCalls.pairs_tsv,
            candidate_counts = SampleShardCalls.candidate_count,
            count = count,
            random_seed = random_seed,
            prefix = prefix,
            docker = utils_docker,
            runtime_attr_override = runtime_attr_merge
    }

    output {
        File variant_sample_pairs = MergeRandomCalls.pairs_tsv
        File candidate_summary = MergeRandomCalls.summary_tsv
    }
}

task SampleShardCalls {
    input {
        File vcf
        File vcf_idx
        Int count
        Float? min_af
        Float? max_af
        Int? min_ac
        Int? max_ac
        Boolean singleton
        Array[String] filter_values
        Array[String] allele_types
        Int? min_allele_length
        Int? max_allele_length
        Array[String] variant_types
        Array[String] include_samples
        Array[String] exclude_samples
        Int random_seed
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        join_by() {
            local sep="$1"
            shift
            local first="$1"
            shift
            printf '%s' "$first" "${@/#/$sep}"
        }

        EXPR_PARTS=()

        MIN_AF="~{default='' min_af}"
        [ -n "$MIN_AF" ] && EXPR_PARTS+=("INFO/AF>=$MIN_AF")
        MAX_AF="~{default='' max_af}"
        [ -n "$MAX_AF" ] && EXPR_PARTS+=("INFO/AF<=$MAX_AF")

        if [ "~{true='1' false='0' singleton}" = "1" ]; then
            EXPR_PARTS+=("INFO/AC=1")
        else
            MIN_AC="~{default='' min_ac}"
            [ -n "$MIN_AC" ] && EXPR_PARTS+=("INFO/AC>=$MIN_AC")
            MAX_AC="~{default='' max_ac}"
            [ -n "$MAX_AC" ] && EXPR_PARTS+=("INFO/AC<=$MAX_AC")
        fi

        FILTER_VALUES="~{sep=',' filter_values}"
        if [ -n "$FILTER_VALUES" ]; then
            IFS=',' read -ra FVALS <<< "$FILTER_VALUES"
            FSUB=()
            for v in "${FVALS[@]}"; do FSUB+=("FILTER=\"$v\""); done
            EXPR_PARTS+=("($(join_by '||' "${FSUB[@]}"))")
        fi

        ALLELE_TYPES="~{sep=',' allele_types}"
        if [ -n "$ALLELE_TYPES" ]; then
            IFS=',' read -ra AVALS <<< "$ALLELE_TYPES"
            ASUB=()
            for v in "${AVALS[@]}"; do ASUB+=("INFO/allele_type=\"$v\""); done
            EXPR_PARTS+=("($(join_by '||' "${ASUB[@]}"))")
        fi

        MIN_ALLELE_LENGTH="~{default='' min_allele_length}"
        [ -n "$MIN_ALLELE_LENGTH" ] && EXPR_PARTS+=("abs(INFO/allele_length)>=$MIN_ALLELE_LENGTH")
        MAX_ALLELE_LENGTH="~{default='' max_allele_length}"
        [ -n "$MAX_ALLELE_LENGTH" ] && EXPR_PARTS+=("abs(INFO/allele_length)<=$MAX_ALLELE_LENGTH")

        VARIANT_TYPES="~{sep=',' variant_types}"
        if [ -n "$VARIANT_TYPES" ]; then
            IFS=',' read -ra TVALS <<< "$VARIANT_TYPES"
            TSUB=()
            for v in "${TVALS[@]}"; do
                bt=$(echo "$v" | tr '[:upper:]' '[:lower:]')
                [ "$bt" = "snv" ] && bt="snp"
                TSUB+=("TYPE=\"$bt\"")
            done
            EXPR_PARTS+=("($(join_by '||' "${TSUB[@]}"))")
        fi

        # push every site-level filter down to bcftools (compiled C, single streaming pass)
        # so pysam below only has to walk the pre-shrunk set to do the per-sample GT check
        if [ ${#EXPR_PARTS[@]} -gt 0 ]; then
            bcftools view -i "$(join_by '&&' "${EXPR_PARTS[@]}")" -Ob -o filtered.bcf ~{vcf}
        else
            bcftools view -Ob -o filtered.bcf ~{vcf}
        fi

        python3 <<'PYCODE'
import random
import pysam

MIN_AF = "~{default='' min_af}"
MAX_AF = "~{default='' max_af}"
MIN_AC = "~{default='' min_ac}"
MAX_AC = "~{default='' max_ac}"
SINGLETON = "~{true='1' false='0' singleton}" == "1"
INCLUDE_SAMPLES = [v for v in "~{sep=',' include_samples}".split(",") if v]
EXCLUDE_SAMPLES = set(v for v in "~{sep=',' exclude_samples}".split(",") if v)
COUNT = ~{count}
SEED = ~{random_seed}

min_af = float(MIN_AF) if MIN_AF else None
max_af = float(MAX_AF) if MAX_AF else None
min_ac = int(MIN_AC) if MIN_AC else None
max_ac = int(MAX_AC) if MAX_AC else None
if SINGLETON:
    min_ac = max_ac = 1

def allele_passes(idx, ac_field, af_field):
    # AC/AF are Number=A (one value per ALT, e.g. multiallelic VAMOS/TR sites) on some
    # callsets, plain scalars on others. The bcftools prefilter above only checked whether
    # ANY element of the vector satisfied the bound, which is correct for keeping the
    # record but wrong for attributing a hit to a specific sample: a sample carrying a
    # common allele (e.g. AC=29) at a site that also has an unrelated true singleton ALT
    # must not be reported as a singleton carrier. Index into the vector by the sample's
    # actual ALT index (1-based in GT, so idx-1 into the tuple); a scalar applies as-is.
    ac = ac_field[idx - 1] if isinstance(ac_field, tuple) else ac_field
    af = af_field[idx - 1] if isinstance(af_field, tuple) else af_field
    if min_ac is not None and (ac is None or ac < min_ac):
        return False
    if max_ac is not None and (ac is None or ac > max_ac):
        return False
    if min_af is not None and (af is None or af < min_af):
        return False
    if max_af is not None and (af is None or af > max_af):
        return False
    return True

rng = random.Random(SEED)
vcf_in = pysam.VariantFile("filtered.bcf")
candidate_samples = INCLUDE_SAMPLES or list(vcf_in.header.samples)
candidate_samples = [s for s in candidate_samples if s not in EXCLUDE_SAMPLES]

reservoir = []
n_seen = 0
for record in vcf_in:
    chrom = record.chrom
    start = record.pos + 1
    end = record.stop
    vid = record.id if record.id else "."
    allele_type = record.info.get("allele_type", ".")
    ac_field = record.info.get("AC")
    af_field = record.info.get("AF")
    need_allele_check = min_ac is not None or max_ac is not None or min_af is not None or max_af is not None

    for sample in candidate_samples:
        gt = record.samples[sample]["GT"]
        if gt is None:
            continue
        alt_indices = {a for a in gt if a is not None and a > 0}
        if not alt_indices:
            continue
        if need_allele_check and not any(allele_passes(idx, ac_field, af_field) for idx in alt_indices):
            continue
        n_seen += 1
        pair = (chrom, start, end, vid, allele_type, sample)
        if len(reservoir) < COUNT:
            reservoir.append(pair)
        else:
            j = rng.randint(0, n_seen - 1)
            if j < COUNT:
                reservoir[j] = pair

reservoir.sort(key=lambda p: (p[0], p[1]))
with open("~{prefix}.pairs.tsv", "w") as out:
    out.write("#chrom\tstart\tend\tID\tallele_type\tsamples\n")
    for chrom, start, end, vid, allele_type, sample in reservoir:
        out.write(f"{chrom}\t{start}\t{end}\t{vid}\t{allele_type}\t{sample}\n")

with open("~{prefix}.candidate_count.txt", "w") as out:
    out.write(str(n_seen) + "\n")
PYCODE
    >>>

    output {
        File pairs_tsv = "~{prefix}.pairs.tsv"
        Int candidate_count = read_int("~{prefix}.candidate_count.txt")
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 4,
        disk_gb: 3 * ceil(size([vcf, vcf_idx], "GB")) + 10,
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

task MergeRandomCalls {
    input {
        Array[File] pair_tsvs
        Array[Int] candidate_counts
        Int count
        Int random_seed
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        python3 <<'PYCODE'
import random

PAIR_FILES = "~{sep=',' pair_tsvs}".split(",")
COUNTS = [int(v) for v in "~{sep=',' candidate_counts}".split(",") if v]
COUNT = ~{count}
SEED = ~{random_seed}

candidates = []
for path in PAIR_FILES:
    with open(path) as handle:
        next(handle)
        candidates.extend(line.rstrip("\n") for line in handle)

rng = random.Random(SEED)
rng.shuffle(candidates)
selected = sorted(candidates[:COUNT], key=lambda line: (line.split("\t")[0], int(line.split("\t")[1])))

with open("~{prefix}.variant_sample_pairs.tsv", "w") as out:
    out.write("#chrom\tstart\tend\tID\tallele_type\tsamples\n")
    for line in selected:
        out.write(line + "\n")

with open("~{prefix}.candidate_summary.txt", "w") as out:
    out.write(f"requested\t{COUNT}\n")
    out.write(f"total_candidates_found\t{sum(COUNTS)}\n")
    out.write(f"pairs_written\t{len(selected)}\n")
PYCODE
    >>>

    output {
        File pairs_tsv = "~{prefix}.variant_sample_pairs.tsv"
        File summary_tsv = "~{prefix}.candidate_summary.txt"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 2,
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
