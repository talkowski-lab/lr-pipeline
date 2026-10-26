version 1.0

import "../annotation/AnnotateRegion.wdl"
import "../utils/ExactMatch.wdl"
import "../utils/Structs.wdl"
import "../utils/TruvariMatch.wdl"

workflow EvaluateCallsetAccuracy {
    meta {
        description: [
            "This utility benchmarks one contig of a cohort VCF against a truth set combined from an SNV & indel truth VCF and an SV truth VCF, reporting precision and recall at the site and call level. Each truth VCF is cut to the callset's samples and contig, split into biallelic records and left-aligned, and the SNV & indel truth is restricted to `truth_snv_indel_bed` when given. SV truth records identical in CHROM, POS, REF and ALT to an SNV & indel truth record are dropped before the two are combined, so each variant is counted once. The BED also restricts callset records shorter than `sv_min_length`, matching the truth records it bounds; longer callset records are kept genome-wide like the SV truth. Only records with at least one non-reference call are kept on either side.",
            "Matching reuses the first round of `AnnotateCallsetOverlap`, an exact match on CHROM, POS, REF and ALT, followed by a single Truvari pass at 0.5 sequence similarity, size similarity and overlap on the records left unmatched on both sides, admitting callset records with `abs(INFO/allele_length)` of at least `min_sv_length_truvari_vcf` and truth records of at least `min_sv_length_truvari_truth_vcf`. Every record is therefore paired at most once. Truvari skips records longer than 50 kb, so those count only through the exact match.",
            "A site is one record and a call is one sample's non-reference genotype at it, with phasing and zygosity ignored. A callset site is a true positive when it matched a truth record and a callset call when that sample is also non-reference in the matched truth record; a truth site or call is recovered the same way from the other side. Precision counts callset sites and calls, stratified by `INFO/allele_type` and `abs(INFO/allele_length)`, and recall counts truth sites and calls, stratified by the type and length of the truth record - `SVTYPE` and `SVLEN` when present, otherwise REF and ALT. A matched pair can therefore fall in different strata on the two sides. Truth records that are not an SNV, insertion or deletion, such as multi-nucleotide variants and inversions, are left out of recall.",
            "Each record is also assigned the genomic region of `AnnotateRegion` - simple repeat, segmental duplication, RepeatMasker or unique sequence - and rows are stratified into `US/RM` and `SD/SR`, by the callset record's region for precision and the truth record's for recall.",
            "Every row is reported at both matching stages, cumulatively, so the `truvari` rows give the totals. All count columns sum across contigs."
        ]
    }

    parameter_meta {
        vcf: "Callset VCF to benchmark, carrying `INFO/allele_type` and `INFO/allele_length`."
        vcf_idx: "Index for `vcf`."
        truth_snv_indel_vcf: "Truth VCF of SNVs and indels, holding every callset sample."
        truth_snv_indel_vcf_idx: "Index for `truth_snv_indel_vcf`."
        truth_sv_vcf: "Truth VCF of SVs, holding every callset sample."
        truth_sv_vcf_idx: "Index for `truth_sv_vcf`."
        truth_snv_indel_bed: "Confident regions of `truth_snv_indel_vcf`, restricting it and the callset records shorter than `sv_min_length`."
        contig: "Contig to evaluate."
        args_string_vcf: "`bcftools view` include expression selecting the callset records to benchmark."
        sv_min_length: "Callset `abs(INFO/allele_length)` from which records are kept outside `truth_snv_indel_bed`."
        min_sv_length_truvari_vcf: "Minimum `abs(INFO/allele_length)` for an unmatched callset record to enter the Truvari pass."
        min_sv_length_truvari_truth_vcf: "Minimum absolute length for an unmatched truth record to enter the Truvari pass."
        length_bins: "Size-bin edges stratifying records by absolute length."
        shard_bin_size_exact_match: "Width in base pairs of the regions the exact-match round is sharded into."
        records_per_shard_region: "Number of records per shard when assigning genomic regions."
        ref_fa: "From references."
        ref_fai: "From references."
        simple_repeats_bed: "From references."
        seg_dup_bed: "From references."
        repeat_masked_bed: "From references."
        accuracy_tsv: "One row per level (`site` or `call`), sample, allele type, length bin, region (`US/RM` or `SD/SR`) and cumulative matching stage (`exact` or `truvari`), with `tp_eval`, `fp`, `tp_truth`, `fn`, precision, recall and F1; strata with no records on either side are left out. Site rows use sample `all`; type, length bin and region each have an `all` marginal."
    }

    input {
        File vcf
        File vcf_idx
        File truth_snv_indel_vcf
        File truth_snv_indel_vcf_idx
        File truth_sv_vcf
        File truth_sv_vcf_idx
        File? truth_snv_indel_bed
        String contig
        String prefix

        String args_string_vcf = 'INFO/allele_type!="trv"'
        Int sv_min_length = 50
        Int min_sv_length_truvari_vcf = 5
        Int min_sv_length_truvari_truth_vcf = 10
        Array[Int] length_bins = [0, 1, 50, 500]

        Int shard_bin_size_exact_match = 5000000
        Int records_per_shard_region = 10000

        File ref_fa
        File ref_fai
        File simple_repeats_bed
        File seg_dup_bed
        File repeat_masked_bed

        String utils_docker

        RuntimeAttr? runtime_attr_prepare_eval
        RuntimeAttr? runtime_attr_prepare_truth
        RuntimeAttr? runtime_attr_drop_matched_truth
        RuntimeAttr? runtime_attr_run_truvari
        RuntimeAttr? runtime_attr_count_accuracy
        RuntimeAttr? runtime_attr_create_exact_shards
        RuntimeAttr? runtime_attr_subset_exact_vcf
        RuntimeAttr? runtime_attr_subset_exact_truth
        RuntimeAttr? runtime_attr_exact_match
        RuntimeAttr? runtime_attr_append_exact_annotations
        RuntimeAttr? runtime_attr_truvari_subset_vcf
        RuntimeAttr? runtime_attr_truvari_subset_truth
        RuntimeAttr? runtime_attr_concat_exact_annotations
        RuntimeAttr? runtime_attr_concat_exact_unmatched
        RuntimeAttr? runtime_attr_concat_exact_truth
        RuntimeAttr? runtime_attr_region_shard
        RuntimeAttr? runtime_attr_annotate_region
        RuntimeAttr? runtime_attr_concat_region_shards
    }

    call PrepareAccuracyEvalVcf {
        input:
            vcf = vcf,
            vcf_idx = vcf_idx,
            contig = contig,
            include_args = args_string_vcf,
            sv_min_length = sv_min_length,
            restrict_bed = truth_snv_indel_bed,
            prefix = "~{prefix}.~{contig}.eval",
            docker = utils_docker,
            runtime_attr_override = runtime_attr_prepare_eval
    }

    call PrepareAccuracyTruthVcf {
        input:
            snv_indel_vcf = truth_snv_indel_vcf,
            snv_indel_vcf_idx = truth_snv_indel_vcf_idx,
            sv_vcf = truth_sv_vcf,
            sv_vcf_idx = truth_sv_vcf_idx,
            contig = contig,
            samples = PrepareAccuracyEvalVcf.samples,
            restrict_bed = truth_snv_indel_bed,
            ref_fa = ref_fa,
            ref_fai = ref_fai,
            prefix = "~{prefix}.~{contig}.truth",
            docker = utils_docker,
            runtime_attr_override = runtime_attr_prepare_truth
    }

    call ExactMatch.ExactMatch {
        input:
            vcf = PrepareAccuracyEvalVcf.eval_vcf,
            vcf_idx = PrepareAccuracyEvalVcf.eval_vcf_idx,
            truth_snv_indel_vcf = PrepareAccuracyTruthVcf.truth_vcf,
            truth_snv_indel_vcf_idx = PrepareAccuracyTruthVcf.truth_vcf_idx,
            contig = contig,
            prefix = "~{prefix}.~{contig}",
            run_exact_matching = true,
            run_truvari_matching = true,
            shard_bin_size_exact_match = shard_bin_size_exact_match,
            min_sv_length_truvari_vcf = min_sv_length_truvari_vcf,
            min_sv_length_truvari_truth_snv_indel_vcf = min_sv_length_truvari_truth_vcf,
            length_field_vcf = "allele_length",
            length_field_truth_snv_indel_vcf = "allele_length",
            source_tag_truth_snv_indel_vcf = "truth",
            utils_docker = utils_docker,
            runtime_attr_create_exact_shards = runtime_attr_create_exact_shards,
            runtime_attr_subset_exact_vcf = runtime_attr_subset_exact_vcf,
            runtime_attr_subset_exact_truth = runtime_attr_subset_exact_truth,
            runtime_attr_exact_match = runtime_attr_exact_match,
            runtime_attr_append_exact_annotations = runtime_attr_append_exact_annotations,
            runtime_attr_truvari_subset_vcf = runtime_attr_truvari_subset_vcf,
            runtime_attr_truvari_subset_truth = runtime_attr_truvari_subset_truth,
            runtime_attr_concat_exact_annotations = runtime_attr_concat_exact_annotations,
            runtime_attr_concat_exact_unmatched = runtime_attr_concat_exact_unmatched,
            runtime_attr_concat_exact_truth = runtime_attr_concat_exact_truth
    }

    # ExactMatch hands Truvari every long enough truth record, so remove those already matched to keep pairs one-to-one
    call DropMatchedTruthRecords {
        input:
            vcf = select_first([ExactMatch.truvari_truth_vcf]),
            vcf_idx = select_first([ExactMatch.truvari_truth_vcf_idx]),
            match_tsv = select_first([ExactMatch.annotated_tsv]),
            prefix = "~{prefix}.~{contig}.truvari_truth",
            docker = utils_docker,
            runtime_attr_override = runtime_attr_drop_matched_truth
    }

    call TruvariMatch.RunTruvari {
        input:
            vcf = ExactMatch.truvari_eval_vcf,
            vcf_idx = ExactMatch.truvari_eval_vcf_idx,
            truth_snv_indel_vcf = DropMatchedTruthRecords.unmatched_vcf,
            truth_snv_indel_vcf_idx = DropMatchedTruthRecords.unmatched_vcf_idx,
            pctseq = 0.5,
            pctsize = 0.5,
            pctovl = 0.5,
            sizemin = 0,
            sizefilt = 0,
            tag_value = "TRUVARI",
            source_tag = "truth",
            ref_fa = ref_fa,
            ref_fai = ref_fai,
            prefix = "~{prefix}.~{contig}.truvari",
            docker = utils_docker,
            runtime_attr_override = runtime_attr_run_truvari
    }

    call AnnotateRegion.AnnotateRegion as AnnotateEvalRegion {
        input:
            vcf = PrepareAccuracyEvalVcf.eval_vcf,
            vcf_idx = PrepareAccuracyEvalVcf.eval_vcf_idx,
            contigs = [contig],
            prefix = "~{prefix}.~{contig}.eval",
            records_per_shard = records_per_shard_region,
            simple_repeats_bed = simple_repeats_bed,
            seg_dup_bed = seg_dup_bed,
            repeat_masked_bed = repeat_masked_bed,
            utils_docker = utils_docker,
            runtime_attr_shard = runtime_attr_region_shard,
            runtime_attr_annotate_region = runtime_attr_annotate_region,
            runtime_attr_concat_shards = runtime_attr_concat_region_shards
    }

    call AnnotateRegion.AnnotateRegion as AnnotateTruthRegion {
        input:
            vcf = PrepareAccuracyTruthVcf.truth_vcf,
            vcf_idx = PrepareAccuracyTruthVcf.truth_vcf_idx,
            contigs = [contig],
            prefix = "~{prefix}.~{contig}.truth",
            records_per_shard = records_per_shard_region,
            simple_repeats_bed = simple_repeats_bed,
            seg_dup_bed = seg_dup_bed,
            repeat_masked_bed = repeat_masked_bed,
            utils_docker = utils_docker,
            runtime_attr_shard = runtime_attr_region_shard,
            runtime_attr_annotate_region = runtime_attr_annotate_region,
            runtime_attr_concat_shards = runtime_attr_concat_region_shards
    }

    call CountCallsetAccuracy {
        input:
            eval_vcf = PrepareAccuracyEvalVcf.eval_vcf,
            eval_vcf_idx = PrepareAccuracyEvalVcf.eval_vcf_idx,
            truth_vcf = PrepareAccuracyTruthVcf.truth_vcf,
            truth_vcf_idx = PrepareAccuracyTruthVcf.truth_vcf_idx,
            eval_region_tsv = AnnotateEvalRegion.annotations_tsv_region,
            truth_region_tsv = AnnotateTruthRegion.annotations_tsv_region,
            match_tsvs = select_all([ExactMatch.annotated_tsv, RunTruvari.annotation_tsv]),
            length_bins = length_bins,
            prefix = "~{prefix}.~{contig}.callset_accuracy",
            docker = utils_docker,
            runtime_attr_override = runtime_attr_count_accuracy
    }

    output {
        File accuracy_tsv = CountCallsetAccuracy.accuracy_tsv
    }
}

task PrepareAccuracyEvalVcf {
    input {
        File vcf
        File vcf_idx
        String contig
        String include_args
        Int sv_min_length
        File? restrict_bed
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        bcftools query -l ~{vcf} > samples.txt

        # Short records are bounded by the BED like the SNV & indel truth, long ones are kept genome-wide like the SV truth
        bcftools view \
            -r ~{contig} \
            ~{if defined(restrict_bed) then "-T " + restrict_bed else ""} \
            -i '(~{include_args}) && abs(INFO/allele_length)<~{sv_min_length}' \
            -Oz -o short.vcf.gz \
            ~{vcf}
        bcftools view \
            -r ~{contig} \
            -i '(~{include_args}) && abs(INFO/allele_length)>=~{sv_min_length}' \
            -Oz -o long.vcf.gz \
            ~{vcf}
        tabix -p vcf short.vcf.gz
        tabix -p vcf long.vcf.gz

        # Drop the cohort-wide AC/AN so --min-ac counts the genotypes themselves
        bcftools concat -a short.vcf.gz long.vcf.gz -Ou \
            | bcftools annotate -x INFO/AC,INFO/AN,INFO/AF -Ou \
            | bcftools view --min-ac 1 -Oz -o ~{prefix}.vcf.gz

        tabix -p vcf -f ~{prefix}.vcf.gz
    >>>

    output {
        File eval_vcf = "~{prefix}.vcf.gz"
        File eval_vcf_idx = "~{prefix}.vcf.gz.tbi"
        Array[String] samples = read_lines("samples.txt")
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 4,
        disk_gb: 3 * ceil(size(vcf, "GB")) + 10,
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

task PrepareAccuracyTruthVcf {
    input {
        File snv_indel_vcf
        File snv_indel_vcf_idx
        File sv_vcf
        File sv_vcf_idx
        String contig
        Array[String] samples
        File? restrict_bed
        File ref_fa
        File ref_fai
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        # Split and left-align so alleles match the callset exactly, warning on placeholder N REFs; -s puts samples in the callset's order for concat, and sort fixes records left-aligned past their neighbours
        bcftools view -r ~{contig} -s ~{sep=',' samples} -Ou ~{snv_indel_vcf} \
            | bcftools annotate -x INFO/AC,INFO/AN,INFO/AF -Ou \
            | bcftools norm -m- -f ~{ref_fa} -c w -Ou \
            | bcftools view --min-ac 1 ~{if defined(restrict_bed) then "-T " + restrict_bed else ""} -Ou \
            | bcftools sort -T sort_tmp -Oz -o snv_indel.vcf.gz
        bcftools view -r ~{contig} -s ~{sep=',' samples} -Ou ~{sv_vcf} \
            | bcftools annotate -x INFO/AC,INFO/AN,INFO/AF -Ou \
            | bcftools norm -m- -f ~{ref_fa} -c w -Ou \
            | bcftools view --min-ac 1 -Ou \
            | bcftools sort -T sort_tmp -Oz -o sv.vcf.gz
        tabix -p vcf snv_indel.vcf.gz
        tabix -p vcf sv.vcf.gz

        # Keep one copy of a variant present in both truths, dropping the SV truth's
        bcftools isec -C -c none -w1 sv.vcf.gz snv_indel.vcf.gz -Oz -o sv.unique.vcf.gz
        tabix -p vcf sv.unique.vcf.gz
        echo "Dropped $(( $(bcftools index -n sv.vcf.gz) - $(bcftools index -n sv.unique.vcf.gz) )) SV truth records also in the SNV & indel truth"

        bcftools concat -a snv_indel.vcf.gz sv.unique.vcf.gz -Oz -o combined.vcf.gz

        python3 <<'CODE'
import pysam

# Give every record a unique ID and the allele type and signed length the callset carries, from SVTYPE and SVLEN when present
with pysam.VariantFile('combined.vcf.gz') as vcf_in:
    header = vcf_in.header.copy()
    header.add_line('##INFO=<ID=allele_type,Number=1,Type=String,Description="Allele type">')
    header.add_line('##INFO=<ID=allele_length,Number=1,Type=Integer,Description="Allele length">')
    with pysam.VariantFile("~{prefix}.vcf", 'w', header=header) as vcf_out:
        for n, record in enumerate(vcf_in, 1):
            record.translate(header)
            ref, alt = record.ref, record.alts[0]
            if 'SVTYPE' in record.info:
                allele_type = record.info['SVTYPE'].lower()
                svlen = record.info.get('SVLEN')
                svlen = svlen[0] if isinstance(svlen, tuple) else svlen
                length = abs(svlen) if svlen is not None else abs(len(alt) - len(ref))
                length = -length if allele_type == 'del' else length
            elif len(alt) == len(ref):
                allele_type, length = ('snv' if len(ref) == 1 else 'mnv'), 0
            else:
                allele_type, length = ('ins' if len(alt) > len(ref) else 'del'), len(alt) - len(ref)
            record.id = f"truth_~{contig}_{n}"
            record.info['allele_type'] = allele_type
            record.info['allele_length'] = length
            vcf_out.write(record)

pysam.tabix_index("~{prefix}.vcf", preset='vcf', force=True)
CODE
    >>>

    output {
        File truth_vcf = "~{prefix}.vcf.gz"
        File truth_vcf_idx = "~{prefix}.vcf.gz.tbi"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 4,
        disk_gb: ceil(3 * (size(snv_indel_vcf, "GB") + size(sv_vcf, "GB")) + size(ref_fa, "GB")) + 10,
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

task DropMatchedTruthRecords {
    input {
        File vcf
        File vcf_idx
        File match_tsv
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        # The truth ID is the seventh column of the exact-match annotations
        grep -v '^#' ~{match_tsv} | cut -f7 > matched_ids.txt || true

        bcftools view -e 'ID=@matched_ids.txt' -Oz -o ~{prefix}.vcf.gz ~{vcf}
        tabix -p vcf -f ~{prefix}.vcf.gz
    >>>

    output {
        File unmatched_vcf = "~{prefix}.vcf.gz"
        File unmatched_vcf_idx = "~{prefix}.vcf.gz.tbi"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 4,
        disk_gb: 2 * ceil(size(vcf, "GB") + size(match_tsv, "GB")) + 10,
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

task CountCallsetAccuracy {
    input {
        File eval_vcf
        File eval_vcf_idx
        File truth_vcf
        File truth_vcf_idx
        File eval_region_tsv
        File truth_region_tsv
        Array[File] match_tsvs
        Array[Int] length_bins
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        python3 <<'CODE'
import csv
import sys
from collections import defaultdict

import pysam

LENGTH_BINS = [~{sep=", " length_bins}]
SIZE_LABELS = [f'{start}-{end - 1}' if end - 1 > start else f'{start}' for start, end in zip(LENGTH_BINS, LENGTH_BINS[1:])]
SIZE_LABELS.append(f'{LENGTH_BINS[-1]}+')
# Matching stages in the order they run, as tagged by ExactMatch and RunTruvari
STAGES = ['EXACT', 'TRUVARI']
TYPES = ['snv', 'ins', 'del']
REGION_GROUPS = {'US': 'US/RM', 'RM': 'US/RM', 'SD': 'SD/SR', 'SR': 'SD/SR'}
REGIONS = ['US/RM', 'SD/SR']
COUNT_COLUMNS = ['tp_eval', 'fp', 'tp_truth', 'fn']

def size_bin(length):
    for index, label in enumerate(SIZE_LABELS):
        if index + 1 == len(LENGTH_BINS) or length < LENGTH_BINS[index + 1]:
            return label

def first(value):
    return value[0] if isinstance(value, tuple) else value

def carriers(record):
    return frozenset(name for name, sample in record.samples.items() if any(a not in (None, 0) for a in sample['GT'] or ()))

def load(path, region_path):
    with open(region_path) as handle:
        regions = {fields[4]: REGION_GROUPS[fields[5]] for fields in (line.rstrip('\n').split('\t') for line in handle)}
    with pysam.VariantFile(path) as vcf:
        records = {}
        for record in vcf:
            if record.id not in regions:
                sys.exit(f'Record {record.id} has no genomic region')
            allele_type = str(first(record.info['allele_type'])).lower()
            records[record.id] = (allele_type, size_bin(abs(int(first(record.info['allele_length'])))), regions[record.id], carriers(record))
        return list(vcf.header.samples), records

samples, evals = load("~{eval_vcf}", "~{eval_region_tsv}")
_, truths = load("~{truth_vcf}", "~{truth_region_tsv}")
skipped = defaultdict(int)
for allele_type, _, _, _ in truths.values():
    if allele_type not in TYPES:
        skipped[allele_type] += 1

# Each callset record is matched at most once, by the first stage that pairs it
matches = {}
for path in "~{sep=',' match_tsvs}".split(','):
    with open(path) as handle:
        for line in handle:
            if line.startswith('#'):
                continue
            fields = line.rstrip('\n').split('\t')
            eval_id, stage, truth_id = fields[4], fields[5], fields[6]
            if eval_id in matches:
                sys.exit(f'Callset record {eval_id} is matched more than once')
            matches[eval_id] = (STAGES.index(stage), truth_id)

counts = defaultdict(lambda: [0] * len(COUNT_COLUMNS))

def add(stage, level, sample, allele_type, size, region, column):
    for type_key in (allele_type, 'all'):
        for size_key in (size, 'all'):
            for region_key in (region, 'all'):
                counts[(stage, level, sample, type_key, size_key, region_key)][COUNT_COLUMNS.index(column)] += 1

for stage_index, stage in enumerate(STAGES):
    recovered = defaultdict(set)
    for eval_id, (allele_type, size, region, eval_carriers) in evals.items():
        match = matches.get(eval_id)
        truth_id = match[1] if match is not None and match[0] <= stage_index else None
        matched_carriers = truths[truth_id][3] if truth_id is not None else frozenset()
        if truth_id is not None:
            recovered[truth_id] |= eval_carriers
        add(stage, 'site', 'all', allele_type, size, region, 'tp_eval' if truth_id is not None else 'fp')
        for sample in eval_carriers:
            column = 'tp_eval' if sample in matched_carriers else 'fp'
            add(stage, 'call', sample, allele_type, size, region, column)
            add(stage, 'call', 'all', allele_type, size, region, column)
    for truth_id, (allele_type, size, region, truth_carriers) in truths.items():
        if allele_type not in TYPES:
            continue
        found = recovered.get(truth_id)
        add(stage, 'site', 'all', allele_type, size, region, 'tp_truth' if found is not None else 'fn')
        for sample in truth_carriers:
            column = 'tp_truth' if found is not None and sample in found else 'fn'
            add(stage, 'call', sample, allele_type, size, region, column)
            add(stage, 'call', 'all', allele_type, size, region, column)

def ratio(numerator, denominator):
    return numerator / denominator if denominator else None

def fmt(value):
    return '.' if value is None else f'{value:.4f}'

types = TYPES + sorted({allele_type for allele_type, _, _, _ in evals.values()} - set(TYPES)) + ['all']
with open("~{prefix}.tsv", 'w', newline='') as handle:
    writer = csv.writer(handle, delimiter='\t', lineterminator='\n')
    writer.writerow(['level', 'sample', 'allele_type', 'length_bin', 'region', 'match_stage'] + COUNT_COLUMNS + ['precision', 'recall', 'f1'])
    for stage in STAGES:
        for level, level_samples in (('site', ['all']), ('call', samples + ['all'])):
            for sample in level_samples:
                for allele_type in types:
                    for size in SIZE_LABELS + ['all']:
                        for region in REGIONS + ['all']:
                            tp_eval, fp, tp_truth, fn = counts[(stage, level, sample, allele_type, size, region)]
                            # A stratum empty at one stage is empty at every stage, such as SNVs above 0 bp
                            if tp_eval + fp + tp_truth + fn == 0:
                                continue
                            precision = ratio(tp_eval, tp_eval + fp)
                            recall = ratio(tp_truth, tp_truth + fn)
                            f1 = ratio(2 * precision * recall, precision + recall) if precision is not None and recall is not None else None
                            writer.writerow([level, sample, allele_type, size, region, stage.lower(),
                                             tp_eval, fp, tp_truth, fn, fmt(precision), fmt(recall), fmt(f1)])

print(f'{len(evals)} callset and {len(truths)} truth records, {len(matches)} matched; truth left out of recall: {dict(skipped)}')
CODE
    >>>

    output {
        File accuracy_tsv = "~{prefix}.tsv"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 8,
        disk_gb: ceil(size(eval_vcf, "GB") + size(truth_vcf, "GB") + size([eval_region_tsv, truth_region_tsv], "GB") + size(match_tsvs, "GB")) + 10,
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
