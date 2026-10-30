version 1.0

import "../annotation/AnnotateRegion.wdl"
import "../utils/ExactMatch.wdl"
import "../utils/Structs.wdl"

workflow EvaluateCallsetAccuracy {
    meta {
        description: [
            "This utility benchmarks a callset against a truth VCF, reporting precision and recall at the site and call level across a set of contigs. Each contig's callset VCF is evaluated on its own, in parallel, and the counts are summed across contigs into one table and a set of plots. The truth VCF holds every contig and is read once per contig; it may be paired with a second truth VCF of SVs, whose records identical in CHROM, POS, REF and ALT to a `truth_vcf` record are dropped before the two are combined, so each variant is counted once.",
            "Both sides are cut to the contig and the callset's samples, split into biallelic records and left-aligned, and filtered by their own `bcftools view` include expressions; only records with at least one non-reference call are kept. Every record is given an allele type and signed length, from `INFO/allele_type` and `INFO/allele_length` when present, then from `SVTYPE` and `SVLEN`, then from REF and ALT, with a duplication written as an insertion (a sequence ALT longer than its REF) counted and matched as an insertion, so any VCF can be benchmarked. `truth_bed` restricts `truth_vcf` to its confident regions and bounds the callset records it covers: every callset record, or only those shorter than `sv_min_length` when `truth_sv_vcf` supplies the longer variants genome-wide.",
            "Matching reuses the first round of `AnnotateCallsetOverlap`, an exact match on CHROM, POS, REF and ALT, followed by a single Truvari pass between the callset records left unmatched with an absolute length of at least `min_sv_length_truvari_vcf` and every truth record of at least `min_sv_length_truvari_truth_vcf`, up to the `max_sv_length_truvari_*` lengths, at the `truvari_*` thresholds. The Truvari pass keeps every pair that passes these thresholds rather than one partner per record, so a record can match several records on the other side, and each record still counts once in its own side's totals. Callset records above `max_sv_length_truvari_vcf` count only through the exact match.",
            "A site is one record and a call is one sample's non-reference genotype at it, with phasing and zygosity ignored. A callset site is a true positive when it matched any truth record and a callset call when that sample is also non-reference in any matched truth record; a truth site or call is recovered the same way from the other side. Precision counts callset sites and calls, stratified by allele type and absolute length, and recall counts truth sites and calls, stratified by the type and length of the truth record, so a matched pair can fall in different strata on the two sides. Truth records that are not an SNV, insertion or deletion, such as multi-nucleotide variants, inversions and symbolic duplications, are left out of recall.",
            "Each record is also assigned the genomic region of `AnnotateRegion` - simple repeat, segmental duplication, RepeatMasker or unique sequence - and rows are stratified into `US/RM` and `SD/SR`, by the callset record's region for precision and the truth record's for recall. Every row is reported at both matching stages, cumulatively, so the `truvari` rows give the totals."
        ]
    }

    parameter_meta {
        vcfs: "Callset VCFs to benchmark, one per contig."
        vcf_idxs: "Indexes for `vcfs`."
        contigs: "Contig of each of `vcfs`, in the same order."
        truth_vcf: "Truth VCF covering every contig in `contigs` and holding every callset sample."
        truth_vcf_idx: "Index for `truth_vcf`."
        truth_sv_vcf: "Optional second truth VCF of SVs, holding every callset sample and kept genome-wide."
        truth_sv_vcf_idx: "Index for `truth_sv_vcf`."
        truth_bed: "Confident regions of `truth_vcf`, restricting it and the callset records it covers."
        args_string_vcf: "`bcftools view` include expression selecting the callset records to benchmark, applied after typing, so it may use `INFO/allele_type` and `INFO/allele_length`."
        args_string_truth_vcf: "`bcftools view` include expression selecting the `truth_vcf` records, applied to them as read."
        args_string_truth_sv_vcf: "`bcftools view` include expression selecting the `truth_sv_vcf` records, applied to them as read."
        sv_min_length: "Absolute length from which callset records are kept outside `truth_bed` when `truth_sv_vcf` is given."
        min_sv_length_truvari_vcf: "Minimum absolute length for an unmatched callset record to enter the Truvari pass."
        min_sv_length_truvari_truth_vcf: "Minimum absolute length for a truth record to enter the Truvari pass."
        max_sv_length_truvari_vcf: "Maximum absolute length for an unmatched callset record to enter the Truvari pass."
        max_sv_length_truvari_truth_vcf: "Maximum absolute length for a truth record to enter the Truvari pass."
        truvari_refdist: "Maximum distance in base pairs between a Truvari pair's breakpoints."
        truvari_pctseq: "Minimum Truvari sequence similarity of a pair; 0 compares by position and size only."
        truvari_pctsize: "Minimum Truvari size similarity of a pair."
        truvari_pctovl: "Minimum Truvari reciprocal overlap of a pair."
        truvari_typeignore: "Let Truvari pair records of different SV types."
        length_bins: "Size-bin edges stratifying records by absolute length."
        shard_bin_size_exact_match: "Width in base pairs of the regions each contig's exact-match round is sharded into."
        records_per_shard_region: "Number of records per shard when assigning genomic regions."
        ref_fa: "From references."
        ref_fai: "From references."
        simple_repeats_bed: "From references."
        seg_dup_bed: "From references."
        repeat_masked_bed: "From references."
        exact_match_tsvs: "Exact-match pairs of each contig, in the order of `contigs`, one row per callset record with the callset record ID in the fifth column and the truth record ID in the seventh."
        truvari_match_tsvs: "Truvari pairs of each contig in the same layout, one row per passing callset and truth record pair."
        accuracy_tsv: "Counts summed across contigs, one row per level (`site` or `call`), sample, allele type, length bin, region (`US/RM` or `SD/SR`) and cumulative matching stage (`exact` or `truvari`), with `tp_eval`, `fp`, `tp_truth`, `fn`, precision, recall and F1; strata with no records on either side are left out. Site rows use sample `all`; type, length bin and region each have an `all` marginal."
        accuracy_per_contig_tsv: "The rows of `accuracy_tsv` for each contig, under a leading `contig` column."
        accuracy_plots_tar: "Tarball of precision and recall plots by allele type, length bin and region, by sample and by matching stage, with `summary.tsv` holding the site and call metrics side by side."
    }

    input {
        Array[File] vcfs
        Array[File] vcf_idxs
        Array[String] contigs
        File truth_vcf
        File truth_vcf_idx
        File? truth_sv_vcf
        File? truth_sv_vcf_idx
        File? truth_bed
        String prefix

        String? args_string_vcf
        String? args_string_truth_vcf
        String? args_string_truth_sv_vcf
        Int sv_min_length = 50
        Int min_sv_length_truvari_vcf = 5
        Int min_sv_length_truvari_truth_vcf = 10
        Int max_sv_length_truvari_vcf = 50000
        Int max_sv_length_truvari_truth_vcf = 50000
        Int truvari_refdist = 500
        Float truvari_pctseq = 0.5
        Float truvari_pctsize = 0.5
        Float truvari_pctovl = 0.5
        Boolean truvari_typeignore = false
        Array[Int] length_bins = [0, 1, 50, 2000]

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
        RuntimeAttr? runtime_attr_truvari_match
        RuntimeAttr? runtime_attr_count_accuracy
        RuntimeAttr? runtime_attr_summarize_accuracy
        RuntimeAttr? runtime_attr_create_exact_shards
        RuntimeAttr? runtime_attr_subset_exact_vcf
        RuntimeAttr? runtime_attr_subset_exact_truth
        RuntimeAttr? runtime_attr_exact_match
        RuntimeAttr? runtime_attr_append_exact_annotations
        RuntimeAttr? runtime_attr_truvari_subset_vcf
        RuntimeAttr? runtime_attr_truvari_subset_truth
        RuntimeAttr? runtime_attr_truvari_subset_oversized
        RuntimeAttr? runtime_attr_concat_exact_annotations
        RuntimeAttr? runtime_attr_concat_exact_unmatched
        RuntimeAttr? runtime_attr_concat_exact_truth
        RuntimeAttr? runtime_attr_concat_exact_oversized
        RuntimeAttr? runtime_attr_region_shard
        RuntimeAttr? runtime_attr_annotate_region
        RuntimeAttr? runtime_attr_concat_region_shards
    }

    # Without a separate SV truth, the BED's truth holds every variant length, so it bounds every callset record
    Int bounded_below_length = if defined(truth_sv_vcf) then sv_min_length else 2147483647

    scatter (i in range(length(vcfs))) {
        String contig = contigs[i]

        call PrepareAccuracyEvalVcf {
            input:
                vcf = vcfs[i],
                vcf_idx = vcf_idxs[i],
                contig = contig,
                include_args = args_string_vcf,
                bounded_below_length = bounded_below_length,
                restrict_bed = truth_bed,
                ref_fa = ref_fa,
                ref_fai = ref_fai,
                prefix = "~{prefix}.~{contig}.eval",
                docker = utils_docker,
                runtime_attr_override = runtime_attr_prepare_eval
        }

        call PrepareAccuracyTruthVcf {
            input:
                truth_vcf = truth_vcf,
                truth_vcf_idx = truth_vcf_idx,
                sv_vcf = truth_sv_vcf,
                sv_vcf_idx = truth_sv_vcf_idx,
                include_args = args_string_truth_vcf,
                sv_include_args = args_string_truth_sv_vcf,
                contig = contig,
                samples = PrepareAccuracyEvalVcf.samples,
                restrict_bed = truth_bed,
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
                truth_snv_indel_vcf = PrepareAccuracyTruthVcf.prepared_vcf,
                truth_snv_indel_vcf_idx = PrepareAccuracyTruthVcf.prepared_vcf_idx,
                contig = contig,
                prefix = "~{prefix}.~{contig}",
                run_exact_matching = true,
                run_truvari_matching = true,
                shard_bin_size_exact_match = shard_bin_size_exact_match,
                min_sv_length_truvari_vcf = min_sv_length_truvari_vcf,
                min_sv_length_truvari_truth_snv_indel_vcf = min_sv_length_truvari_truth_vcf,
                max_sv_length_truvari_vcf = max_sv_length_truvari_vcf,
                max_sv_length_truvari_truth_snv_indel_vcf = max_sv_length_truvari_truth_vcf,
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
                runtime_attr_truvari_subset_oversized = runtime_attr_truvari_subset_oversized,
                runtime_attr_concat_exact_annotations = runtime_attr_concat_exact_annotations,
                runtime_attr_concat_exact_unmatched = runtime_attr_concat_exact_unmatched,
                runtime_attr_concat_exact_truth = runtime_attr_concat_exact_truth,
                runtime_attr_concat_exact_oversized = runtime_attr_concat_exact_oversized
        }

        # ExactMatch hands Truvari every long enough truth record, so a truth record matched exactly can still pair with other callset records
        call TruvariMatchAllPairs {
            input:
                vcf = ExactMatch.truvari_eval_vcf,
                vcf_idx = ExactMatch.truvari_eval_vcf_idx,
                truth_vcf = select_first([ExactMatch.truvari_truth_vcf]),
                truth_vcf_idx = select_first([ExactMatch.truvari_truth_vcf_idx]),
                refdist = truvari_refdist,
                pctseq = truvari_pctseq,
                pctsize = truvari_pctsize,
                pctovl = truvari_pctovl,
                sizemax = if max_sv_length_truvari_vcf > max_sv_length_truvari_truth_vcf then max_sv_length_truvari_vcf else max_sv_length_truvari_truth_vcf,
                typeignore = truvari_typeignore,
                ref_fa = ref_fa,
                ref_fai = ref_fai,
                prefix = "~{prefix}.~{contig}.truvari",
                docker = utils_docker,
                runtime_attr_override = runtime_attr_truvari_match
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
                vcf = PrepareAccuracyTruthVcf.prepared_vcf,
                vcf_idx = PrepareAccuracyTruthVcf.prepared_vcf_idx,
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
                truth_vcf = PrepareAccuracyTruthVcf.prepared_vcf,
                truth_vcf_idx = PrepareAccuracyTruthVcf.prepared_vcf_idx,
                eval_region_tsv = AnnotateEvalRegion.annotations_tsv_region,
                truth_region_tsv = AnnotateTruthRegion.annotations_tsv_region,
                match_tsvs = select_all([ExactMatch.annotated_tsv, TruvariMatchAllPairs.annotation_tsv]),
                length_bins = length_bins,
                prefix = "~{prefix}.~{contig}.callset_accuracy",
                docker = utils_docker,
                runtime_attr_override = runtime_attr_count_accuracy
        }

        File exact_match_tsv = select_first([ExactMatch.annotated_tsv])
    }

    call SummarizeCallsetAccuracy {
        input:
            accuracy_tsvs = CountCallsetAccuracy.accuracy_tsv,
            contigs = contigs,
            prefix = "~{prefix}.callset_accuracy",
            docker = utils_docker,
            runtime_attr_override = runtime_attr_summarize_accuracy
    }

    output {
        Array[File] exact_match_tsvs = exact_match_tsv
        Array[File] truvari_match_tsvs = TruvariMatchAllPairs.annotation_tsv
        File accuracy_tsv = SummarizeCallsetAccuracy.combined_tsv
        File accuracy_per_contig_tsv = SummarizeCallsetAccuracy.per_contig_tsv
        File accuracy_plots_tar = SummarizeCallsetAccuracy.plots_tar
    }
}

task PrepareAccuracyEvalVcf {
    input {
        File vcf
        File vcf_idx
        String contig
        String? include_args
        Int bounded_below_length
        File? restrict_bed
        File ref_fa
        File ref_fai
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        bcftools query -l ~{vcf} > samples.txt

        # Split and left-align as the truth is, dropping the cohort-wide AC/AN so --min-ac counts the genotypes themselves
        bcftools view -r ~{contig} -Ou ~{vcf} \
            | bcftools annotate -x INFO/AC,INFO/AN,INFO/AF -Ou \
            | bcftools norm -m- -f ~{ref_fa} -c w -Ou \
            | bcftools sort -T sort_tmp -Oz -o split.vcf.gz
        python3 /opt/scripts/qc/annotate_accuracy_alleles.py --vcf split.vcf.gz --output typed.vcf.gz --id_prefix "eval_~{contig}"

        # Records the BED's truth covers are bounded by it; longer ones, covered by a separate SV truth, are kept genome-wide
        bcftools view \
            ~{if defined(restrict_bed) then "-T " + restrict_bed else ""} \
            -i '~{"(" + include_args + ") && "}abs(INFO/allele_length)<~{bounded_below_length}' \
            --min-ac 1 \
            -Oz -o short.vcf.gz \
            typed.vcf.gz
        bcftools view \
            -i '~{"(" + include_args + ") && "}abs(INFO/allele_length)>=~{bounded_below_length}' \
            --min-ac 1 \
            -Oz -o long.vcf.gz \
            typed.vcf.gz
        tabix -p vcf short.vcf.gz
        tabix -p vcf long.vcf.gz
        bcftools concat -a short.vcf.gz long.vcf.gz -Oz -o ~{prefix}.vcf.gz
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
        disk_gb: ceil(4 * size(vcf, "GB") + size(ref_fa, "GB")) + 10,
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
        File truth_vcf
        File truth_vcf_idx
        File? sv_vcf
        File? sv_vcf_idx
        String? include_args
        String? sv_include_args
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
        bcftools view -r ~{contig} -s ~{sep=',' samples} ~{"-i '" + include_args + "'"} -Ou ~{truth_vcf} \
            | bcftools annotate -x INFO/AC,INFO/AN,INFO/AF -Ou \
            | bcftools norm -m- -f ~{ref_fa} -c w -Ou \
            | bcftools view --min-ac 1 ~{if defined(restrict_bed) then "-T " + restrict_bed else ""} -Ou \
            | bcftools sort -T sort_tmp -Oz -o truth.vcf.gz
        tabix -p vcf truth.vcf.gz

        if [[ -n "~{sv_vcf}" ]]; then
            bcftools view -r ~{contig} -s ~{sep=',' samples} ~{"-i '" + sv_include_args + "'"} -Ou ~{sv_vcf} \
                | bcftools annotate -x INFO/AC,INFO/AN,INFO/AF -Ou \
                | bcftools norm -m- -f ~{ref_fa} -c w -Ou \
                | bcftools view --min-ac 1 -Ou \
                | bcftools sort -T sort_tmp -Oz -o sv.vcf.gz
            tabix -p vcf sv.vcf.gz

            # Keep one copy of a variant present in both truths, dropping the SV truth's
            bcftools isec -C -c none -w1 sv.vcf.gz truth.vcf.gz -Oz -o sv.unique.vcf.gz
            tabix -p vcf sv.unique.vcf.gz
            echo "Dropped $(( $(bcftools index -n sv.vcf.gz) - $(bcftools index -n sv.unique.vcf.gz) )) SV truth records also in the truth VCF"
            bcftools concat -a truth.vcf.gz sv.unique.vcf.gz -Oz -o combined.vcf.gz
        else
            mv truth.vcf.gz combined.vcf.gz
        fi

        python3 /opt/scripts/qc/annotate_accuracy_alleles.py --vcf combined.vcf.gz --output ~{prefix}.vcf.gz --id_prefix "truth_~{contig}" --replace_ids
    >>>

    output {
        File prepared_vcf = "~{prefix}.vcf.gz"
        File prepared_vcf_idx = "~{prefix}.vcf.gz.tbi"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 4,
        disk_gb: ceil(3 * (size(truth_vcf, "GB") + size(sv_vcf, "GB")) + size(ref_fa, "GB")) + 10,
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

task TruvariMatchAllPairs {
    input {
        File vcf
        File vcf_idx
        File truth_vcf
        File truth_vcf_idx
        Int refdist
        Float pctseq
        Float pctsize
        Float pctovl
        Int sizemax
        Boolean typeignore
        File ref_fa
        File ref_fai
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        # Truvari looks for each index beside its VCF, but Cromwell localizes separately globbed shard indexes elsewhere
        if [[ "~{vcf_idx}" != "~{vcf}.tbi" ]]; then
            ln -sf "~{vcf_idx}" "~{vcf}.tbi"
        fi
        if [[ "~{truth_vcf_idx}" != "~{truth_vcf}.tbi" ]]; then
            ln -sf "~{truth_vcf_idx}" "~{truth_vcf}.tbi"
        fi

        python3 <<'CODE'
import truvari
from truvari.bench import check_inputs, parse_args

# Build the parameters as `truvari bench` does, so its filters and defaults apply unchanged
args = parse_args(['-b', "~{truth_vcf}", '-c', "~{vcf}", '-o', 'unused', '--reference', "~{ref_fa}",
                   '--refdist', "~{refdist}", '--pctseq', "~{pctseq}", '--pctsize', "~{pctsize}", '--pctovl', "~{pctovl}",
                   '--sizemin', '0', '--sizefilt', '0', '--sizemax', "~{sizemax}"~{if typeignore then ", '--typeignore'" else ""}])
check_inputs(args)
params = truvari.VariantParams(args, skip_gt=args.bSample is None or args.cSample is None,
                               short_circuit=args.short, decompose=args.no_decompose)
bench = truvari.Bench(params)
base = truvari.VariantFile("~{truth_vcf}", params=params)
comp = truvari.VariantFile("~{vcf}", params=params)
tree = truvari.build_region_tree(base, comp)
truvari.merge_region_tree_overlaps(tree)

# Keep every pair passing the thresholds instead of letting a picker choose one partner per record
pairs = set()
for chunk, chunk_id in truvari.chunker(params, ('base', base.fetch_regions(tree)), ('comp', comp.fetch_regions(tree))):
    if not chunk['base'] or not chunk['comp']:
        continue
    for row in bench.build_matrix(chunk['base'], chunk['comp'], chunk_id):
        for match in row:
            if match.state:
                record = match.comp
                pairs.add((record.chrom, record.pos + 1, record.ref, record.alts[0], record.id, match.base.id))

with open("~{prefix}.annotation.tsv", 'w') as handle:
    for chrom, pos, ref, alt, eval_id, truth_id in sorted(pairs):
        handle.write(f'{chrom}\t{pos}\t{ref}\t{alt}\t{eval_id}\tTRUVARI\t{truth_id}\ttruth\t.\n')
print(f'{len(pairs)} pairs between {len({p[4] for p in pairs})} callset and {len({p[5] for p in pairs})} truth records')
CODE
    >>>

    output {
        File annotation_tsv = "~{prefix}.annotation.tsv"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 4,
        disk_gb: 2 * ceil(size(vcf, "GB") + size(truth_vcf, "GB") + size(ref_fa, "GB")) + 10,
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

# A callset record can pair with several truth records, each pair kept at the first stage that finds it
matches = defaultdict(dict)
for path in "~{sep=',' match_tsvs}".split(','):
    with open(path) as handle:
        for line in handle:
            if line.startswith('#'):
                continue
            fields = line.rstrip('\n').split('\t')
            eval_id, stage, truth_id = fields[4], fields[5], fields[6]
            matches[eval_id].setdefault(truth_id, STAGES.index(stage))

counts = defaultdict(lambda: [0] * len(COUNT_COLUMNS))

def add(stage, level, sample, allele_type, size, region, column):
    for type_key in (allele_type, 'all'):
        for size_key in (size, 'all'):
            for region_key in (region, 'all'):
                counts[(stage, level, sample, type_key, size_key, region_key)][COUNT_COLUMNS.index(column)] += 1

for stage_index, stage in enumerate(STAGES):
    recovered = defaultdict(set)
    for eval_id, (allele_type, size, region, eval_carriers) in evals.items():
        truth_ids = [truth_id for truth_id, index in matches.get(eval_id, {}).items() if index <= stage_index]
        matched_carriers = frozenset().union(*(truths[truth_id][3] for truth_id in truth_ids))
        for truth_id in truth_ids:
            recovered[truth_id] |= eval_carriers
        add(stage, 'site', 'all', allele_type, size, region, 'tp_eval' if truth_ids else 'fp')
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

pairs = sum(len(matched) for matched in matches.values())
print(f'{len(evals)} callset and {len(truths)} truth records, {len(matches)} callset records in {pairs} pairs; truth left out of recall: {dict(skipped)}')
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

task SummarizeCallsetAccuracy {
    input {
        Array[File] accuracy_tsvs
        Array[String] contigs
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        python3 /opt/scripts/qc/plot_callset_accuracy.py \
            --tsvs ~{sep=' ' accuracy_tsvs} \
            --contigs ~{sep=' ' contigs} \
            --combined_tsv ~{prefix}.tsv \
            --per_contig_tsv ~{prefix}.per_contig.tsv \
            --output_dir ~{prefix}_plots

        tar -czf ~{prefix}_plots.tar.gz ~{prefix}_plots
    >>>

    output {
        File combined_tsv = "~{prefix}.tsv"
        File per_contig_tsv = "~{prefix}.per_contig.tsv"
        File plots_tar = "~{prefix}_plots.tar.gz"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 4,
        disk_gb: 2 * ceil(size(accuracy_tsvs, "GB")) + 10,
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
