version 1.0

import "../annotation/AnnotateInSilicoPredictors.wdl"
import "../annotation/AnnotateGQMetrics.wdl"
import "../annotation/AnnotateRegion.wdl"
import "../annotation/AnnotateSQMetrics.wdl"
import "../annotation/AnnotateVRS.wdl"
import "../tools/MergeHiPhaseCallsets.wdl"
import "../utils/Helpers.wdl"
import "../utils/Structs.wdl"
import "AnnotateVcf.wdl"

workflow PostProcessTRLociHPRCHGSVC {
    input {
        File vcf
        File vcf_idx
        String contig
        Array[File] trgt_vcfs
        Array[File] trgt_vcf_idxs
        Array[String] sample_ids
        Array[File] base_vcfs
        Array[File] base_vcf_idxs
        String prefix

        Boolean run_flag_homopolymer_trvs
        Boolean run_normalize_ploidy
        Boolean replace_gnomad_str
        File? ped
        File? swap_samples_base
        Int max_phase_edit_distance = 10
        Float max_phase_edit_distance_pct = 10.0

        File gnomad_tr_json
        File ref_fa
        File ref_fai
        File seqrepo_tar
        File simple_repeats_bed
        File seg_dup_bed
        File repeat_masked_bed
        String cadd_ht
        String pangolin_ht
        String phylop_ht
        String revel_ht
        String spliceai_ht
        String annotate_in_silico_predictors_script = "https://raw.githubusercontent.com/talkowski-lab/lr-annotation/main/scripts/annotation/annotate_insilico_predictors.py"
        String genome_build = "GRCh38"

        String utils_docker
        String trgt_docker
        String vrs_docker
        String hail_docker

        RuntimeAttr? runtime_attr_discover
        RuntimeAttr? runtime_attr_subset_trgt
        RuntimeAttr? runtime_attr_add_trgt_end
        RuntimeAttr? runtime_attr_fix_AL_header
        RuntimeAttr? runtime_attr_merge_trgt
        RuntimeAttr? runtime_attr_filter_merged
        RuntimeAttr? runtime_attr_prepare
        RuntimeAttr? runtime_attr_phase_replacements
        RuntimeAttr? runtime_attr_replacement_sq_metrics
        RuntimeAttr? runtime_attr_replacement_sd_metrics
        RuntimeAttr? runtime_attr_replacement_ab_metrics
        RuntimeAttr? runtime_attr_vrs_annotate
        RuntimeAttr? runtime_attr_vrs_extract
        RuntimeAttr? runtime_attr_region_annotate
        RuntimeAttr? runtime_attr_insilico_annotate
        RuntimeAttr? runtime_attr_attach_annotations
        RuntimeAttr? runtime_attr_apply
    }

    Boolean do_normalize_ploidy = run_normalize_ploidy && defined(ped)

    call DiscoverTRLoci {
        input:
            vcf = vcf,
            vcf_idx = vcf_idx,
            trgt_vcfs = trgt_vcfs,
            trgt_vcf_idxs = trgt_vcf_idxs,
            sample_ids = sample_ids,
            gnomad_tr_json = gnomad_tr_json,
            contig = contig,
            prefix = "~{prefix}.~{contig}",
            docker = utils_docker,
            runtime_attr_override = runtime_attr_discover
    }

    if (DiscoverTRLoci.trgt_candidate_locus_count > 0) {
        scatter (i in range(length(trgt_vcfs))) {
            call SubsetTRGTForCatalogLoci {
                input:
                    vcf = trgt_vcfs[i],
                    vcf_idx = trgt_vcf_idxs[i],
                    match_keys = DiscoverTRLoci.trgt_match_keys,
                    prefix = "~{prefix}.~{contig}.trgt.~{i}",
                    docker = utils_docker,
                    runtime_attr_override = runtime_attr_subset_trgt
            }

            call Helpers.AddTREndTag as AddTRGTEndTag {
                input:
                    vcf = SubsetTRGTForCatalogLoci.subset_vcf,
                    vcf_idx = SubsetTRGTForCatalogLoci.subset_vcf_idx,
                    prefix = "~{prefix}.~{contig}.trgt.~{i}.with_end",
                    docker = utils_docker,
                    runtime_attr_override = runtime_attr_add_trgt_end
            }

            call MergeHiPhaseCallsets.FixALHeader {
                input:
                    vcf = AddTRGTEndTag.vcf_with_end,
                    vcf_idx = AddTRGTEndTag.vcf_with_end_idx,
                    prefix = "~{prefix}.~{contig}.trgt.~{i}.merge_ready",
                    docker = utils_docker,
                    runtime_attr_override = runtime_attr_fix_AL_header
            }
        }

        call MergeHiPhaseCallsets.TRGTMergeContig {
            input:
                vcfs = FixALHeader.fixed_vcf,
                vcf_idxs = FixALHeader.fixed_vcf_idx,
                contig = contig,
                ref_fa = ref_fa,
                ref_fai = ref_fai,
                prefix = "~{prefix}.~{contig}.missing_loci",
                docker = trgt_docker,
                runtime_attr_override = runtime_attr_merge_trgt
        }

        call KeepMergedTRGTWithAC {
            input:
                vcf = TRGTMergeContig.merged_vcf,
                vcf_idx = TRGTMergeContig.merged_vcf_idx,
                prefix = "~{prefix}.~{contig}.missing_loci.ac",
                docker = utils_docker,
                runtime_attr_override = runtime_attr_filter_merged
        }

        if (KeepMergedTRGTWithAC.retained_count > 0) {
            call PrepareReplacementLoci {
                input:
                    vcf = vcf,
                    vcf_idx = vcf_idx,
                    replacement_vcf = KeepMergedTRGTWithAC.retained_vcf,
                    replacement_vcf_idx = KeepMergedTRGTWithAC.retained_vcf_idx,
                    sample_ids = sample_ids,
                    run_flag_homopolymer_trvs = run_flag_homopolymer_trvs,
                    normalize_ploidy = do_normalize_ploidy,
                    ped = ped,
                    prefix = "~{prefix}.~{contig}.replacement_seed",
                    docker = utils_docker,
                    runtime_attr_override = runtime_attr_prepare
            }

            if (PrepareReplacementLoci.retained_count > 0) {
            call PhaseReplacementLoci {
                input:
                    vcf = PrepareReplacementLoci.prepared_vcf,
                    vcf_idx = PrepareReplacementLoci.prepared_vcf_idx,
                    original_vcf = vcf,
                    original_vcf_idx = vcf_idx,
                    replacement_map_tsv = PrepareReplacementLoci.replacement_map_tsv,
                    base_vcfs = base_vcfs,
                    base_vcf_idxs = base_vcf_idxs,
                    swap_samples_base = swap_samples_base,
                    max_phase_edit_distance = max_phase_edit_distance,
                    max_phase_edit_distance_pct = max_phase_edit_distance_pct,
                    prefix = "~{prefix}.~{contig}.replacement_seed.trv_phasing",
                    docker = utils_docker,
                    runtime_attr_override = runtime_attr_phase_replacements
            }

            call AnnotateSQMetrics.CalculateSiteMetrics as CalculateReplacementSQMetrics {
                input:
                    vcf = PhaseReplacementLoci.phased_vcf,
                    vcf_idx = PhaseReplacementLoci.phased_vcf_idx,
                    prefix = "~{prefix}.~{contig}.replacement.sq_metrics",
                    docker = utils_docker,
                    runtime_attr_override = runtime_attr_replacement_sq_metrics
            }

            call AnnotateGQMetrics.GenerateGQAnnotationTsv as CalculateReplacementSDMetrics {
                input:
                    vcf = PhaseReplacementLoci.phased_vcf,
                    vcf_idx = PhaseReplacementLoci.phased_vcf_idx,
                    gq_field = "SD",
                    gq_bins = [0, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 85, 90, 95, 100],
                    gq_variant_filter = ".",
                    gq_larger_field = false,
                    prefix = "~{prefix}.~{contig}.replacement.sd_metrics",
                    docker = utils_docker,
                    runtime_attr_override = runtime_attr_replacement_sd_metrics
            }

            call AnnotateGQMetrics.GenerateABAnnotationTsv as CalculateReplacementABMetrics {
                input:
                    vcf = PhaseReplacementLoci.phased_vcf,
                    vcf_idx = PhaseReplacementLoci.phased_vcf_idx,
                    ab_bins = [0.00, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45, 0.50, 0.55, 0.60, 0.65, 0.70, 0.75, 0.80, 0.85, 0.90, 0.95, 1.00],
                    prefix = "~{prefix}.~{contig}.replacement.ab_metrics",
                    docker = utils_docker,
                    runtime_attr_override = runtime_attr_replacement_ab_metrics
            }

            call AnnotateVRS.AnnotateVcfWithVRS as AnnotateReplacementVRS {
                input:
                    vcf = PhaseReplacementLoci.phased_vcf,
                    vcf_idx = PhaseReplacementLoci.phased_vcf_idx,
                    prefix = "~{prefix}.~{contig}.replacement.vrs",
                    seqrepo_tar = seqrepo_tar,
                    docker = vrs_docker,
                    runtime_attr_override = runtime_attr_vrs_annotate
            }

            call AnnotateVRS.ExtractVRSAnnotations as ExtractReplacementVRS {
                input:
                    vcf = AnnotateReplacementVRS.annotated_vcf,
                    vcf_idx = AnnotateReplacementVRS.annotated_vcf_idx,
                    prefix = "~{prefix}.~{contig}.replacement.vrs",
                    docker = utils_docker,
                    runtime_attr_override = runtime_attr_vrs_extract
            }

            call AnnotateRegion.AnnotateGenomicContext as AnnotateReplacementRegion {
                input:
                    vcf = PhaseReplacementLoci.phased_vcf,
                    vcf_idx = PhaseReplacementLoci.phased_vcf_idx,
                    simple_repeats_bed = simple_repeats_bed,
                    seg_dup_bed = seg_dup_bed,
                    repeat_masked_bed = repeat_masked_bed,
                    prefix = "~{prefix}.~{contig}.replacement.region",
                    docker = utils_docker,
                    runtime_attr_override = runtime_attr_region_annotate
            }

            call AnnotateInSilicoPredictors.AnnotateInSilicoPredictorsTask as AnnotateReplacementInSilico {
                input:
                    vcf = PhaseReplacementLoci.phased_vcf,
                    vcf_idx = PhaseReplacementLoci.phased_vcf_idx,
                    cadd_ht = cadd_ht,
                    pangolin_ht = pangolin_ht,
                    phylop_ht = phylop_ht,
                    revel_ht = revel_ht,
                    spliceai_ht = spliceai_ht,
                    annotate_in_silico_predictors_script = annotate_in_silico_predictors_script,
                    genome_build = genome_build,
                    prefix = "~{prefix}.~{contig}.replacement.in_silico",
                    docker = hail_docker,
                    runtime_attr_override = runtime_attr_insilico_annotate
            }

            call AnnotateVcf.AnnotateSequentially as AttachReplacementAnnotations {
                input:
                    vcf = PhaseReplacementLoci.phased_vcf,
                    vcf_idx = PhaseReplacementLoci.phased_vcf_idx,
                    annotations_tsvs = [AnnotateReplacementRegion.annotations_tsv, AnnotateReplacementInSilico.annotations_tsv, ExtractReplacementVRS.annotations_tsv, CalculateReplacementSQMetrics.annotations_tsv, CalculateReplacementSDMetrics.annotation_tsv, CalculateReplacementABMetrics.annotation_tsv],
                    prefix = "~{prefix}.~{contig}.replacement_annotated",
                    info_names = [["REGION"], ["cadd_raw_score", "cadd_phred", "pangolin_largest", "revel_max", "phylop", "spliceai_ds_max"], ["VRS_Allele_IDs", "VRS_Error", "VRS_Starts", "VRS_Ends", "VRS_States", "VRS_Lengths", "VRS_RepeatSubunitLengths"], ["inbreeding_coeff", "AS_pab_max", "AS_QUALapprox", "AS_QD", "AS_VarDP", "HWE"], ["sd_hist_all_bin_freq", "sd_hist_alt_bin_freq"], ["ab_hist_alt_bin_freq"]],
                    info_descriptions = [["Genomic context of variant"], ["CADD raw score", "CADD PHRED score", "Largest Pangolin delta score", "Maximum REVEL score", "PhyloP score", "Maximum SpliceAI delta score"], ["VRS allele identifiers", "VRS annotation error", "VRS start positions", "VRS end positions", "VRS allele states", "VRS allele lengths", "VRS repeat subunit lengths"], ["Inbreeding coefficient, the excess heterozygosity at a variant site, computed as 1 - (the number of heterozygous genotypes)/(the number of heterozygous genotypes expected under Hardy-Weinberg equilibrium).", "Allele-specific maximum p-value over callset for binomial test of observed allele balance for a heterozygous genotype, given expectation of AB=0.5.", "Allele-specific sum of PL[0] values; used to approximate the QUAL score.", "Allele-specific variant call confidence normalized by depth of sample reads supporting a variant.", "Allele-specific depth over variant genotypes (does not include depth of reference samples).", "Hardy-Weinberg equilibrium p-value."], ["Histogram for SD calculated on high quality genotypes; bin edges are: 0|5|10|15|20|25|30|35|40|45|50|55|60|65|70|75|80|85|90|95|100.", "Histogram for SD in heterozygous individuals calculated on high quality genotypes; bin edges are: 0|5|10|15|20|25|30|35|40|45|50|55|60|65|70|75|80|85|90|95|100."], ["Histogram for AB in heterozygous individuals calculated on high quality genotypes; bin edges are: 0.00|0.05|0.10|0.15|0.20|0.25|0.30|0.35|0.40|0.45|0.50|0.55|0.60|0.65|0.70|0.75|0.80|0.85|0.90|0.95|1.00."]],
                    info_types = [["String"], ["Float", "Float", "Float", "Float", "Float", "Float"], ["String", "String", "Integer", "Integer", "String", "String", "String"], ["Float", "Float", "Integer", "Float", "Integer", "Float"], ["String", "String"], ["String"]],
                    info_numbers = [["1"], ["1", "1", "1", "1", "1", "1"], ["R", ".", "R", "R", ".", ".", "."], ["A", "A", "A", "A", "A", "1"], ["1", "1"], ["1"]],
                    subset_vcf_strings = [],
                    awk_tsv_conditions = [],
                    subset_tsv_columns = [[6], [6, 7, 8, 9, 10, 11], [6, 7, 8, 9, 10, 11, 12], [6, 7, 8, 9, 10, 11], [6, 7], [6]],
                    strip_info_fields_per_tsv = [false, false, false, false, false, false],
                    docker = utils_docker,
                    runtime_attr_override = runtime_attr_attach_annotations
            }
            }
        }
    }

    call ApplyTRLocusUpdates {
        input:
            vcf = vcf,
            vcf_idx = vcf_idx,
            replacement_vcf = AttachReplacementAnnotations.annotated_vcf,
            replacement_vcf_idx = AttachReplacementAnnotations.annotated_vcf_idx,
            replacement_map_tsv = PrepareReplacementLoci.replacement_map_tsv,
            input_trv_phasing_summary_tsv = PhaseReplacementLoci.trv_phasing_summary_tsv,
            input_trv_catalog_match_tsv = DiscoverTRLoci.trv_catalog_match_tsv,
            replace_gnomad_str = replace_gnomad_str,
            prefix = "~{prefix}.~{contig}",
            docker = utils_docker,
            runtime_attr_override = runtime_attr_apply
    }

    output {
        File trv_postprocessed_vcf = ApplyTRLocusUpdates.trv_postprocessed_vcf
        File trv_postprocessed_vcf_idx = ApplyTRLocusUpdates.trv_postprocessed_vcf_idx
        File trv_subsetted_vcf = ApplyTRLocusUpdates.trv_subsetted_vcf
        File trv_subsetted_vcf_idx = ApplyTRLocusUpdates.trv_subsetted_vcf_idx
        File trv_catalog_match_tsv = ApplyTRLocusUpdates.trv_catalog_match_tsv
        File trv_phasing_summary_tsv = ApplyTRLocusUpdates.trv_phasing_summary_tsv
    }
}

# Build one contig-level catalog report, then select eligible unmatched TRGT loci for replacement
task DiscoverTRLoci {
    input {
        File vcf
        File vcf_idx
        Array[File] trgt_vcfs
        Array[File] trgt_vcf_idxs
        Array[String] sample_ids
        File gnomad_tr_json
        String contig
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        # Validate sample/file alignment before using indexed windows for strict TRID matching
        python3 <<'PY'
import json
import os
import pysam

def normalize_contig(value):  # noqa: E302
    return value[3:] if value.lower().startswith('chr') else value

def parse_explorer(value):  # noqa: E302
    fields = value.rsplit('-', 3)
    if len(fields) != 4:
        return None
    try:
        return fields[0], int(fields[1]), int(fields[2])
    except ValueError:
        return None

def values(value):  # noqa: E302
    if value is None:
        return []
    return [str(item) for item in value] if isinstance(value, tuple) else [str(value)]

trgt_paths = [line.rstrip('\n') for line in open('~{write_lines(trgt_vcfs)}') if line.strip()]  # noqa: E305
trgt_indexes = [line.rstrip('\n') for line in open('~{write_lines(trgt_vcf_idxs)}') if line.strip()]
sample_ids = [line.rstrip('\n') for line in open('~{write_lines(sample_ids)}') if line.strip()]
if not trgt_paths or len(trgt_paths) != len(trgt_indexes) or len(trgt_paths) != len(sample_ids):
    raise RuntimeError('trgt_vcfs, trgt_vcf_idxs, and sample_ids must be non-empty parallel arrays')
if len(set(sample_ids)) != len(sample_ids):
    raise RuntimeError('sample_ids must be unique')
if not os.path.exists('~{vcf_idx}') or not all(os.path.exists(path) for path in trgt_indexes):
    raise RuntimeError('All VCF indexes must be localized')

def record_end(rec):  # noqa: E302
    return rec.stop if rec.stop is not None else rec.pos + len(rec.ref) - 1

def overlaps(rec, start, stop):  # noqa: E302
    return rec.pos <= stop and record_end(rec) >= start

def trid_text(rec):  # noqa: E302
    # Preserve VCF comma separation within one INFO/TRID field
    return ','.join(values(rec.info.get('TRID')))

def nonref_allele_count(rec):  # noqa: E302
    return sum(
        1
        for sample in rec.samples.values()
        for allele in (sample.get('GT') or [])
        if allele is not None and allele > 0
    )

with open('~{gnomad_tr_json}') as handle:
    catalog = json.load(handle)

# Emit one report row per contig-relevant catalog entry, showing its TRExplorerV1 values together but matching each
loci = []
for entry in catalog:
    if not entry or not entry.get('LocusId'):
        continue
    values_from_catalog = entry.get('TRExplorerV1')
    values_from_catalog = values_from_catalog if isinstance(values_from_catalog, list) else [values_from_catalog]
    explorers = []
    for value in values_from_catalog:
        parsed = parse_explorer(str(value)) if value else None
        if parsed and normalize_contig(parsed[0]) == normalize_contig('~{contig}'):
            explorers.append((str(value), parsed[1], parsed[2]))
    if explorers:
        diseases = entry.get('Diseases')
        loci.append({
            'locus_id': str(entry['LocusId']),
            'explorers': explorers,
            # Only catalog loci with one or more associated diseases are eligible
            'disease_eligible': isinstance(diseases, list) and len(diseases) > 0,
        })

with pysam.VariantFile('~{vcf}') as main:
    if list(main.header.samples) != sample_ids:
        raise RuntimeError('sample_ids must exactly match main VCF sample order')
    for locus in loci:
        locus['main_trids'] = []
        locus['main_strict'] = False
        if not locus['disease_eligible']:
            continue
        seen = set()
        for explorer, start, stop in locus['explorers']:
            for rec in main.fetch('~{contig}', max(0, start - 1), stop + 1):
                if rec.info.get('allele_type') != 'trv' or not overlaps(rec, start, stop):
                    continue
                trid = trid_text(rec)
                if trid and trid not in seen:
                    locus['main_trids'].append(trid)
                    seen.add(trid)
                if explorer in trid:
                    locus['main_strict'] = True

for locus in loci:
    locus['trgt_trids'] = []
    locus['trgt_strict'] = False
    locus['trgt_ac'] = 0
    locus['trgt_seen'] = set()
    # Count a sample's non-reference alleles once even when overlapping TRExplorerV1 windows return the same record
    locus['trgt_strict_seen'] = set()

# Open each per-sample VCF once, then query every eligible catalog interval
for i, path in enumerate(trgt_paths):
    with pysam.VariantFile(path) as trgt:
        if list(trgt.header.samples) != [sample_ids[i]]:
            raise RuntimeError(f'TRGT VCF {path} must contain only sample {sample_ids[i]}')
        for locus in loci:
            # A strict input-VCF match takes precedence, so blank the TRGT report cells and skip the raw-TRGT lookup
            if not locus['disease_eligible'] or locus['main_strict']:
                continue
            for explorer, start, stop in locus['explorers']:
                for rec in trgt.fetch('~{contig}', max(0, start - 1), stop + 1):
                    if not overlaps(rec, start, stop):
                        continue
                    trid = trid_text(rec)
                    if trid and trid not in locus['trgt_seen']:
                        locus['trgt_trids'].append(trid)
                        locus['trgt_seen'].add(trid)
                    if explorer in trid:
                        locus['trgt_strict'] = True
                        record_identity = (
                            i, rec.contig, rec.pos, record_end(rec), rec.ref,
                            tuple(rec.alts or ()), trid,
                        )
                        if record_identity not in locus['trgt_strict_seen']:
                            locus['trgt_ac'] += nonref_allele_count(rec)
                            locus['trgt_strict_seen'].add(record_identity)

with open('~{prefix}.trv_catalog_match.tsv', 'w') as out:
    out.write(
        'locus_id\tTRExplorerV1\thas_diseases\tinput_overlapping_trids'
        '\tinput_has_TRExplorerV1_substring\ttrgt_overlapping_trids'
        '\ttrgt_has_TRExplorerV1_substring\ttrgt_matching_allele_count\n'
    )
    for locus in loci:
        explorers = ','.join(value for value, _, _ in locus['explorers'])
        if not locus['disease_eligible']:
            out.write(f'{locus["locus_id"]}\t{explorers}\tfalse\t.\t.\t.\t.\t.\n')
            continue
        if locus['main_strict']:
            out.write(
                f'{locus["locus_id"]}\t{explorers}\ttrue\t'
                f'{"|".join(locus["main_trids"]) or "."}\ttrue\t\t\t\n'
            )
            continue
        out.write(
            f'{locus["locus_id"]}\t{explorers}\ttrue\t'
            f'{"|".join(locus["main_trids"]) or "."}\t{str(locus["main_strict"]).lower()}\t'
            f'{"|".join(locus["trgt_trids"]) or "."}\t{str(locus["trgt_strict"]).lower()}\t'
            f'{locus["trgt_ac"]}\n'
        )

# Only disease-eligible, input-unmatched, strict TRGT matches with AC >= 1 feed cohort merging and replacement
with open('~{prefix}.trgt_match_keys.txt', 'w') as out:
    for locus in loci:
        if (locus['disease_eligible'] and not locus['main_strict']
                and locus['trgt_strict'] and locus['trgt_ac'] >= 1):
            for explorer, _, _ in locus['explorers']:
                out.write(explorer + '\n')
PY
        wc -l < ~{prefix}.trgt_match_keys.txt
    >>>

    output {
        File trv_catalog_match_tsv = "~{prefix}.trv_catalog_match.tsv"
        File trgt_match_keys = "~{prefix}.trgt_match_keys.txt"
        Int trgt_candidate_locus_count = read_int(stdout())
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 8,
        disk_gb: 3 * ceil(size(vcf, "GB") + size(trgt_vcfs, "GB")) + 10,
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

# Retain only disease-associated fallback loci selected during discovery from one sample VCF
task SubsetTRGTForCatalogLoci {
    input {
        File vcf
        File vcf_idx
        File match_keys
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        # Fetch narrow catalog intervals and confirm strict TRID matches before writing records
        python3 <<'PY'
import pysam

keys = [line.strip() for line in open('~{match_keys}') if line.strip()]

def parse_explorer(value):  # noqa: E302
    fields = value.rsplit('-', 3)
    return fields[0], int(fields[1]), int(fields[2])

def vals(value):  # noqa: E302
    return [str(item) for item in value] if isinstance(value, tuple) else ([] if value is None else [str(value)])

src = pysam.VariantFile('~{vcf}')  # noqa: E305
out = pysam.VariantFile('~{prefix}.vcf.gz', 'wz', header=src.header)
records = {}
for match_key in keys:
    key_contig, start, stop = parse_explorer(match_key)
    fetch_contig = key_contig if key_contig in src.header.contigs else 'chr' + key_contig
    for rec in src.fetch(fetch_contig, max(0, start - 1), stop + 1):
        if match_key in '|'.join(vals(rec.info.get('TRID'))):
            records[(rec.contig, rec.start, rec.stop, rec.ref, rec.alts)] = rec.copy()
contig_order = {name: i for i, name in enumerate(src.header.contigs)}
for rec in sorted(records.values(), key=lambda item: (contig_order[item.contig], item.start, item.stop)):
    out.write(rec)
src.close()
out.close()
PY
        tabix -f -p vcf ~{prefix}.vcf.gz
    >>>
    output {
        File subset_vcf = "~{prefix}.vcf.gz"
        File subset_vcf_idx = "~{prefix}.vcf.gz.tbi"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 4,
        disk_gb: 3 * ceil(size(vcf, "GB")) + 5,
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

# Recompute cohort allele counts and retain AC-positive merged loci
task KeepMergedTRGTWithAC {
    input {
        File vcf
        File vcf_idx
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        # Calculate AC directly from merged genotypes rather than trusting upstream INFO/AC
        python3 <<'PY'
import pysam

src = pysam.VariantFile('~{vcf}')  # noqa: E305
if 'AC' not in src.header.info:
    src.header.add_meta(
        'INFO',
        items=[('ID', 'AC'), ('Number', 'A'), ('Type', 'Integer'), ('Description', 'Number of alleles observed')],
    )
out = pysam.VariantFile('~{prefix}.vcf.gz', 'wz', header=src.header)
kept = 0
for rec in src:
    counts = [0] * len(rec.alts or [])
    for sample in rec.samples.values():
        gt = sample.get('GT')
        if gt:
            for allele in gt:
                if allele is not None and 0 < allele <= len(counts):
                    counts[allele - 1] += 1
    ac = sum(counts)
    if ac >= 1:
        rec.info['AC'] = tuple(counts)
        out.write(rec)
        kept += 1
src.close()
out.close()
PY
        tabix -f -p vcf ~{prefix}.vcf.gz
        bcftools view -H ~{prefix}.vcf.gz | wc -l
    >>>
    output {
        File retained_vcf = "~{prefix}.vcf.gz"
        File retained_vcf_idx = "~{prefix}.vcf.gz.tbi"
        Int retained_count = read_int(stdout())
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 8,
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

# Reconcile replacement headers and select old TRVs by positive overlap
task PrepareReplacementLoci {
    input {
        File vcf
        File vcf_idx
        File replacement_vcf
        File replacement_vcf_idx
        Array[String] sample_ids
        Boolean run_flag_homopolymer_trvs
        Boolean normalize_ploidy
        File? ped
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        # Match FORMAT/AL to the main VCF before records move between pysam headers
        python3 <<'PY'
import subprocess

main_header = subprocess.check_output(['bcftools', 'view', '-h', '~{vcf}'], text=True).splitlines()
incoming_header = subprocess.check_output(['bcftools', 'view', '-h', '~{replacement_vcf}'], text=True).splitlines()
main_al = next((line for line in main_header if line.startswith('##FORMAT=<ID=AL,')), None)
if main_al:
    incoming_header = [main_al if line.startswith('##FORMAT=<ID=AL,') else line for line in incoming_header]
with open('replacement.header', 'w') as handle:
    handle.write('\n'.join(incoming_header) + '\n')
PY

        bcftools reheader -h replacement.header ~{replacement_vcf} \
            | bcftools view -Ov -o replacement.normalized.vcf

        # Require unique positive-overlap replacements; phase after all replacements are prepared
        python3 <<'PY'
import pysam

def record_key(rec):  # noqa: E302
    return rec.id if rec.id and rec.id != '.' else f'{rec.chrom}:{rec.pos}:{rec.ref}:{",".join(rec.alts or [])}'

def record_end(rec):  # noqa: E302
    return rec.stop if rec.stop is not None else rec.pos + len(rec.ref) - 1

def overlap(left, right):  # noqa: E302
    return max(0, min(record_end(left), record_end(right)) - max(left.pos, right.pos) + 1)

def shortest_motif_length(rec):  # noqa: E302
    """Match PostprocessCallset: a shortest MOTIFS element of length one is homopolymer."""
    motifs = rec.info.get('MOTIFS')
    if motifs is None:
        return None
    raw_values = motifs if isinstance(motifs, (list, tuple)) else [motifs]
    motif_values = [
        part
        for value in raw_values if value is not None
        for part in str(value).split(',') if part and part != '.'
    ]
    return min((len(motif) for motif in motif_values), default=None)

def parse_ped(path):  # noqa: E302
    """Read PED sex codes using the same mapping as PostprocessCallset."""
    sex_by_sample = {}
    with open(path) as handle:
        for line in handle:
            fields = line.strip().split()
            if not fields:
                continue
            sample_id = fields[1]
            sex_code = fields[4]
            if sex_code == '1':
                sex_by_sample[sample_id] = 'M'
            elif sex_code == '2':
                sex_by_sample[sample_id] = 'F'
            else:
                sex_by_sample[sample_id] = None
    return sex_by_sample

def clear_format_fields(sample_data):  # noqa: E302
    sample_data['GT'] = (None, None)
    sample_data.phased = False

def right_align_unphased(gt):  # noqa: E302
    if gt is None:
        return gt
    return tuple(sorted(gt, key=lambda allele: (allele is not None, allele if allele is not None else -1)))

def make_male_hemizygous(gt, phased):  # noqa: E302
    if gt is None:
        return gt
    alleles = list(gt)
    called_positions = [index for index, allele in enumerate(alleles) if allele is not None]
    if len(called_positions) <= 1:
        return tuple(alleles)
    alt_positions = [index for index, allele in enumerate(alleles) if allele is not None and allele > 0]
    if phased:
        if len(alt_positions) == 1:
            keep_index = alt_positions[0]
        elif alt_positions:
            keep_index = alt_positions[-1]
        else:
            keep_index = called_positions[-1]
        new_gt = [None] * len(alleles)
        new_gt[keep_index] = alleles[keep_index]
        return tuple(new_gt)
    keep_allele = alleles[alt_positions[-1]] if alt_positions else alleles[called_positions[-1]]
    new_gt = [None] * len(alleles)
    new_gt[-1] = keep_allele
    return right_align_unphased(tuple(new_gt))

def normalize_replacement_genotypes(rec):  # noqa: E302
    """Apply replacement-only ploidy normalization before recalculating AC."""
    if not normalize_ploidy:
        return
    for sample in samples:
        sample_data = rec.samples[sample]
        sample_sex = sex_by_sample.get(sample)
        if rec.chrom == 'chrY' and sample_sex == 'F':
            clear_format_fields(sample_data)
            continue
        if rec.chrom in {'chrX', 'chrY'} and sample_sex == 'M':
            sample_data['GT'] = make_male_hemizygous(sample_data.get('GT'), sample_data.phased)
        current_gt = sample_data.get('GT')
        if current_gt is None:
            sample_data['GT'] = (None, None)
        elif len(current_gt) == 1:
            sample_data['GT'] = (None, current_gt[0])
        if not sample_data.phased:
            sample_data['GT'] = right_align_unphased(sample_data.get('GT'))

def recompute_ac(rec):  # noqa: E302
    """Derive allele-specific AC from final replacement GTs, never stale INFO."""
    counts = [0] * len(rec.alts or [])
    for sample_data in rec.samples.values():
        for allele in sample_data.get('GT') or ():
            if allele is not None and 0 < allele <= len(counts):
                counts[allele - 1] += 1
    rec.info['AC'] = tuple(counts)
    return sum(counts)

base = pysam.VariantFile('~{vcf}')  # noqa: E305
incoming = pysam.VariantFile('replacement.normalized.vcf')
if 'allele_type' not in incoming.header.info:
    incoming.header.add_meta(
        'INFO',
        items=[('ID', 'allele_type'), ('Number', '1'), ('Type', 'String'), ('Description', 'Allele type')],
    )
if 'SOURCE' not in incoming.header.info:
    incoming.header.add_meta(
        'INFO',
        items=[('ID', 'SOURCE'), ('Number', '1'), ('Type', 'String'), ('Description', 'Source of variant call')],
    )
if 'allele_length' not in incoming.header.info:
    # Declare allele_length because AnnotateRegion requires it before deriving TRV length from REF
    incoming.header.add_meta(
        'INFO',
        items=[('ID', 'allele_length'), ('Number', '1'), ('Type', 'Integer'), ('Description', 'Allele length')],
    )
if 'AC' not in incoming.header.info:
    incoming.header.add_meta(
        'INFO',
        items=[('ID', 'AC'), ('Number', 'A'), ('Type', 'Integer'), ('Description', 'Number of alleles observed')],
    )
if 'PS' not in incoming.header.formats:
    incoming.header.add_meta(
        'FORMAT',
        items=[('ID', 'PS'), ('Number', '1'), ('Type', 'Integer'), ('Description', 'Phase set')],
    )
run_flag_homopolymer_trvs = ~{true="True" false="False" run_flag_homopolymer_trvs}
normalize_ploidy = ~{true="True" false="False" normalize_ploidy}
sex_by_sample = parse_ped('~{default="NONE" ped}') if normalize_ploidy else {}
if run_flag_homopolymer_trvs and 'HOMOPOLYMER_TRV' not in incoming.header.info:
    incoming.header.info.add(
        'HOMOPOLYMER_TRV',
        0,
        'Flag',
        'Tandem repeat call where the shortest motif has length 1.',
    )
header = incoming.header.copy()
out = pysam.VariantFile('~{prefix}.vcf.gz', 'wz', header=header)
mapping = open('~{prefix}.map.tsv', 'w')
mapping.write('new_record_id\told_record_id\toverlap_bp\tstatus\tphase_summary\told_contig\told_pos\told_end\n')

samples = [line.rstrip('\n') for line in open('~{write_lines(sample_ids)}') if line.strip()]
if list(base.header.samples) != samples or list(header.samples) != samples:
    raise RuntimeError('Main, merged TRGT, and sample_ids sample order must match exactly')
used_old_records = set()
# Match IntegrateTRs.SetTrVariantIds exactly by counting canonical IDs first, then suffixing duplicates in input order
id_counts = {}
id_input = pysam.VariantFile('replacement.normalized.vcf')
for record in id_input:
    normalize_replacement_genotypes(record)
    if recompute_ac(record) < 1:
        continue
    new_id = f'{record.chrom}-{record.pos}-TRV-{len(record.ref) - 1}'
    id_counts[new_id] = id_counts.get(new_id, 0) + 1
id_input.close()
id_seen = {}
for rec in incoming:
    # Recompute AC after replacement-only FORMAT and ploidy changes so a zero-AC record cannot displace an input TRV
    normalize_replacement_genotypes(rec)
    if recompute_ac(rec) < 1:
        continue
    old_candidates = [
        old.copy()
        for old in base.fetch(rec.contig, max(0, rec.start - 1), record_end(rec) + 1)
        if old.info.get('allele_type') == 'trv' and overlap(rec, old) > 0
    ]
    best = max(old_candidates, key=lambda old: overlap(rec, old), default=None)
    best_overlap = overlap(rec, best) if best else 0
    if not best or best_overlap == 0:
        raise RuntimeError(f'No overlapping main VCF TRV found for {record_key(rec)}')
    old_key = record_key(best)
    if old_key in used_old_records:
        raise RuntimeError(f'Multiple replacement records selected main VCF TRV {old_key}')
    used_old_records.add(old_key)
    # Match IntegrateTRs.SetTrVariantIds naming so envelope links use canonical IDs
    new_id = f'{rec.chrom}-{rec.pos}-TRV-{len(rec.ref) - 1}'
    if id_counts[new_id] > 1:
        id_seen[new_id] = id_seen.get(new_id, 0) + 1
        rec.id = f'{new_id}_{id_seen[new_id]}'
    else:
        rec.id = new_id
    rec.info['allele_type'] = 'trv'
    rec.info['SOURCE'] = 'TRExplorer'
    # Flag only replacement records; original main-VCF records never enter this task
    if run_flag_homopolymer_trvs and shortest_motif_length(rec) == 1:
        rec.info['HOMOPOLYMER_TRV'] = True
    # Normalize only replacements to PASS because TRGT leaves FILTER unset, preserving any named filters
    if not tuple(rec.filter.keys()):
        rec.filter.add('PASS')
    mapping.write(
        f'{record_key(rec)}\t{old_key}\t{best_overlap}\treplace\tsee_trv_phasing_summary_tsv'
        f'\t{best.contig}\t{best.pos}\t{record_end(best)}\n'
    )
    out.write(rec)
base.close()
incoming.close()
out.close()
mapping.close()
PY
        rm -f replacement.header replacement.normalized.vcf
        tabix -f -p vcf ~{prefix}.vcf.gz
        bcftools view -H ~{prefix}.vcf.gz | wc -l
    >>>
    output {
        File prepared_vcf = "~{prefix}.vcf.gz"
        File prepared_vcf_idx = "~{prefix}.vcf.gz.tbi"
        File replacement_map_tsv = "~{prefix}.map.tsv"
        Int retained_count = read_int(stdout())
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 2,
        mem_gb: 8,
        disk_gb: 3 * ceil(size(vcf, "GB") + size(replacement_vcf, "GB")) + 10,
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

# Phase eligible replacement calls by sequence agreement with matching base-VCF haplotypes
task PhaseReplacementLoci {
    input {
        File vcf
        File vcf_idx
        File original_vcf
        File original_vcf_idx
        File replacement_map_tsv
        Array[File] base_vcfs
        Array[File] base_vcf_idxs
        File? swap_samples_base
        Int max_phase_edit_distance
        Float max_phase_edit_distance_pct
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        # Fetch each replacement interval from indexed base VCFs; swap map uses BackbonePhase raw-to-canonical form
        python3 <<'PY'
import csv
import os

import edlib
import pysam


def record_key(rec):  # noqa: E302
    return rec.id if rec.id and rec.id != '.' else f'{rec.chrom}:{rec.pos}:{rec.ref}:{",".join(rec.alts or [])}'


def record_end(rec):  # noqa: E302
    return rec.stop if rec.stop is not None else rec.pos + len(rec.ref) - 1


def gt_string(gt, phased=False):  # noqa: E302
    if not gt:
        return '.'
    return ('|' if phased else '/').join('.' if allele is None else str(allele) for allele in gt)


def trid_text(rec):  # noqa: E302
    value = rec.info.get('TRID')
    if value is None:
        return '.'
    return ','.join(str(item) for item in value) if isinstance(value, tuple) else str(value)


def contig_for(handle, contig):  # noqa: E302
    if contig in handle.header.contigs:
        return contig
    alternate = contig[3:] if contig.startswith('chr') else f'chr{contig}'
    return alternate if alternate in handle.header.contigs else None


def normalized_base_gt(call):  # noqa: E302
    """Keep only base calls whose two haplotypes can be unambiguously reconstructed."""
    gt = call.get('GT')
    if not gt or len(gt) != 2:
        return None, 'base_missing_gt'
    if gt == (None, None):
        # Sparse truth records encode no call at this site; reconstruct both haplotypes as REF
        return (0, 0), None
    if any(allele is not None and allele < 0 for allele in gt):
        return None, 'base_invalid_gt'
    if None in gt:
        # Do not fabricate the missing haplotype as reference, since a phased partial call describes only one haplotype
        return None, 'base_partial_gt'
    if not call.phased and gt[0] != gt[1]:
        # Unphased heterozygotes have unknown haplotype orientation
        return None, 'base_ambiguous_unphased_gt'
    return tuple(gt), None


def reconstruct_haplotypes(base_handle, contig, sample, rec):  # noqa: E302
    """Apply complete, non-overlapping base variants to replacement REF for both haplotypes."""
    locus_start = rec.pos
    locus_end = record_end(rec)
    variants = []
    for base_rec in base_handle.fetch(contig, rec.start, rec.stop):
        call = base_rec.samples[sample]
        gt, status = normalized_base_gt(call)
        if status:
            return None, None, status
        # Cohort VCFs contain many overlapping records carried as reference for this sample
        if gt == (0, 0):
            continue
        base_end = record_end(base_rec)
        if base_rec.pos < locus_start or base_end > locus_end:
            return None, None, 'base_boundary_overlapping_variant'
        if any(allele > len(base_rec.alts or []) for allele in gt):
            return None, None, 'base_invalid_allele_index'
        selected_alts = [base_rec.alts[allele - 1] for allele in gt if allele > 0]
        if any(not alt or alt == '*' or alt.startswith('<') or '[' in alt or ']' in alt for alt in selected_alts):
            return None, None, 'base_symbolic_allele'
        offset = base_rec.pos - locus_start
        if rec.ref[offset:offset + len(base_rec.ref)].upper() != base_rec.ref.upper():
            return None, None, 'base_reference_mismatch'
        variants.append((offset, base_rec, gt))
    variants.sort(key=lambda item: item[0])
    output = [[], []]
    for haplotype in range(2):
        cursor = 0
        for offset, base_rec, gt in variants:
            allele = gt[haplotype]
            if allele == 0:
                continue
            if offset < cursor:
                return None, None, 'base_overlapping_variants'
            output[haplotype].append(rec.ref[cursor:offset])
            output[haplotype].append(base_rec.alts[allele - 1])
            cursor = offset + len(base_rec.ref)
        output[haplotype].append(rec.ref[cursor:])
    return ''.join(output[0]), ''.join(output[1]), 'ok'


def reconstruct_available_haploid(base_handle, contig, sample, rec):  # noqa: E302
    """Reconstruct one fully observed biological base haplotype for a haploid replacement call."""
    locus_start = rec.pos
    locus_end = record_end(rec)
    variants = []
    viable_haplotypes = {0, 1}
    for base_rec in base_handle.fetch(contig, rec.start, rec.stop):
        call = base_rec.samples[sample]
        gt = call.get('GT')
        if not gt or len(gt) not in (1, 2):
            return None, None, 'base_missing_gt'
        if all(allele is None for allele in gt):
            # Sparse cohort truth records use fully missing calls for noncarriers
            continue
        if any(allele is not None and allele < 0 for allele in gt):
            return None, None, 'base_invalid_gt'
        if len(gt) == 2 and not call.phased and gt[0] != gt[1]:
            return None, None, 'base_ambiguous_unphased_gt'
        known = {0: gt[0]}
        if len(gt) == 2:
            known[1] = gt[1]
        viable_haplotypes &= {index for index, allele in known.items() if allele is not None}
        if not viable_haplotypes:
            return None, None, 'base_partial_gt_no_complete_haplotype'
        offset = base_rec.pos - locus_start
        variants.append((offset, base_rec, known))
    haplotype = min(viable_haplotypes)
    sequence = []
    cursor = 0
    for offset, base_rec, known in sorted(variants, key=lambda item: item[0]):
        allele = known.get(haplotype)
        if allele is None:
            return None, None, 'base_partial_gt_no_complete_haplotype'
        if allele == 0:
            continue
        base_end = record_end(base_rec)
        if base_rec.pos < locus_start or base_end > locus_end:
            return None, None, 'base_boundary_overlapping_variant'
        if allele > len(base_rec.alts or []):
            return None, None, 'base_invalid_allele_index'
        selected_alt = base_rec.alts[allele - 1]
        if (not selected_alt or selected_alt == '*' or selected_alt.startswith('<')
                or '[' in selected_alt or ']' in selected_alt):
            return None, None, 'base_symbolic_allele'
        offset = base_rec.pos - locus_start
        if rec.ref[offset:offset + len(base_rec.ref)].upper() != base_rec.ref.upper():
            return None, None, 'base_reference_mismatch'
        if offset < cursor:
            return None, None, 'base_overlapping_variants'
        sequence.append(rec.ref[cursor:offset])
        sequence.append(selected_alt)
        cursor = offset + len(base_rec.ref)
    sequence.append(rec.ref[cursor:])
    return ''.join(sequence), haplotype, 'ok'


def distance(left, right):  # noqa: E302
    return edlib.align(left.upper(), right.upper(), task='distance')['editDistance']


def format_distance(total, first, second):  # noqa: E302
    return f'{total} ({first}, {second})'


def combined_edit_distance_pct(first_distance, second_distance, first_left, first_right, second_left, second_right):  # noqa: E302
    """Return the length-weighted percentage across both haplotype pairs."""
    denominator_1 = max(len(first_left), len(first_right), 1)
    denominator_2 = max(len(second_left), len(second_right), 1)
    return 100.0 * (first_distance + second_distance) / (denominator_1 + denominator_2)


def format_distance_pct(first_distance, second_distance, first_left, first_right, second_left, second_right):  # noqa: E302
    """Report weighted total percent, followed by per-haplotype percentages."""
    denominator_1 = max(len(first_left), len(first_right), 1)
    denominator_2 = max(len(second_left), len(second_right), 1)
    first_pct = 100.0 * first_distance / denominator_1
    second_pct = 100.0 * second_distance / denominator_2
    total = combined_edit_distance_pct(
        first_distance, second_distance, first_left, first_right, second_left, second_right)
    return f'{total:.2f}% ({first_pct:.2f}%, {second_pct:.2f}%)'


def audit_status(value):  # noqa: E302
    """Keep reconstruction failures readable in the compact per-call audit."""
    return {
        'base_ambiguous_unphased_gt': 'base_unphased_gt',
        'base_partial_gt_no_complete_haplotype': 'base_partial_gt',
        'base_boundary_overlapping_variant': 'base_boundary_overlap',
        'base_invalid_allele_index': 'base_invalid_allele',
        'base_symbolic_allele': 'base_symbolic',
        'base_reference_mismatch': 'base_ref_mismatch',
        'base_overlapping_variants': 'base_overlap',
    }.get(value, value)


base_paths = [line.strip() for line in open('~{write_lines(base_vcfs)}') if line.strip()]
base_idx_paths = [line.strip() for line in open('~{write_lines(base_vcf_idxs)}') if line.strip()]
if len(base_paths) != len(base_idx_paths):
    raise RuntimeError('base_vcfs and base_vcf_idxs must have equal lengths')
if ~{max_phase_edit_distance} < 0:
    raise RuntimeError('max_phase_edit_distance must be non-negative')
if not 0.0 <= ~{max_phase_edit_distance_pct} <= 100.0:
    raise RuntimeError('max_phase_edit_distance_pct must be between 0 and 100')

sample_swaps = {}
swap_path = '~{swap_samples_base}'
if swap_path and os.path.exists(swap_path):
    with open(swap_path) as source:
        for line in source:
            fields = line.split()
            if not fields:
                continue
            if len(fields) < 2:
                raise RuntimeError('swap_samples_base must contain raw and canonical sample IDs in columns 1 and 2')
            sample_swaps[fields[0]] = fields[1]

# First matching base VCF wins, mirroring BackbonePhase sample assignment
sample_to_base = {}
base_handles = []
for index, path in enumerate(base_paths):
    # Use explicit localized index path; Cromwell need not preserve sibling filenames
    handle = pysam.VariantFile(path, index_filename=base_idx_paths[index])
    base_handles.append(handle)
    canonical_seen = set()
    for raw_sample in handle.header.samples:
        canonical_sample = sample_swaps.get(raw_sample, raw_sample)
        if canonical_sample in canonical_seen:
            raise RuntimeError(f'Duplicate canonical sample {canonical_sample} in base VCF index {index}')
        canonical_seen.add(canonical_sample)
        sample_to_base.setdefault(canonical_sample, (index, raw_sample))

source = pysam.VariantFile('~{vcf}')
original = pysam.VariantFile('~{original_vcf}', index_filename='~{original_vcf_idx}')
if list(original.header.samples) != list(source.header.samples):
    raise RuntimeError('Original main VCF and replacement VCF sample order must match exactly')

# The preparation map links each normalized replacement to its displaced main-VCF TRV
replacement_to_old = {}
with open('~{replacement_map_tsv}') as mapping:
    next(mapping, None)
    for line in mapping:
        fields = line.rstrip('\n').split('\t')
        if len(fields) >= 8 and fields[3] == 'replace':
            replacement_to_old[fields[0]] = (fields[1], fields[5], int(fields[6]), int(fields[7]))

# Indexed lookups avoid scanning a multi-million-record cohort VCF for a handful of displaced TRVs
old_records = {}
for old_key, old_contig, old_pos, old_end in set(replacement_to_old.values()):
    for old_rec in original.fetch(old_contig, max(0, old_pos - 1), old_end):
        if record_key(old_rec) == old_key:
            old_records[old_key] = old_rec.copy()
            break
    if old_key not in old_records:
        raise RuntimeError(f'Replacement map target not found in original VCF: {old_key}')
original.close()

header = source.header.copy()
if 'POSTHOC_BACKBONE_PHASED' not in header.info:
    header.add_meta(
        'INFO',
        items=[('ID', 'POSTHOC_BACKBONE_PHASED'), ('Number', '0'), ('Type', 'Flag'),
               ('Description', 'At least one heterozygous call phased by base-VCF sequence agreement')],
    )
if 'PS' not in header.formats:
    header.add_meta(
        'FORMAT',
        items=[('ID', 'PS'), ('Number', '1'), ('Type', 'Integer'), ('Description', 'Phase set')],
    )

audit_fields = [
    'base_trid', 'replace_trid', 'sample_id', 'base_gt', 'replace_gt',
    'base_hap1_seq', 'base_hap2_seq', 'replace_hap1_seq', 'replace_hap2_seq',
    'edit_dist_aligned', 'edit_dist_unaligned', 'edit_dist_pct',
    'max_edit_dist', 'max_edit_dist_pct', 'final_gt', 'status',
]
output = pysam.VariantFile('~{prefix}.vcf.gz', 'wz', header=header)
audit = open('~{prefix}.trv_phasing_summary.tsv', 'w')
with output, audit:
    writer = csv.DictWriter(audit, fieldnames=audit_fields, delimiter='\t', lineterminator='\n')
    writer.writeheader()
    for rec in source:
        rec.translate(header)
        any_phased = False
        old_mapping = replacement_to_old.get(record_key(rec))
        old_rec = old_records.get(old_mapping[0]) if old_mapping else None
        replacement_trid = trid_text(rec)
        original_trid = trid_text(old_rec) if old_rec is not None else '.'
        for sample in header.samples:
            call = rec.samples[sample]
            gt = call.get('GT')
            # Replacement calls must be unphased unless sequence evidence below selects an orientation
            call.phased = False
            if call.get('PS') is not None:
                call['PS'] = None
            row = dict.fromkeys(audit_fields, '.')
            row.update({
                'base_trid': original_trid,
                'replace_trid': replacement_trid,
                'sample_id': sample,
                'base_gt': gt_string(old_rec.samples[sample].get('GT'), old_rec.samples[sample].phased)
                    if old_rec is not None else '.',
                'replace_gt': gt_string(gt, False),
                'max_edit_dist': ~{max_phase_edit_distance},
                'max_edit_dist_pct': f'{~{max_phase_edit_distance_pct}:.2f}%',
                'final_gt': gt_string(gt, False),
                'status': 'not_replacement' if old_rec is None else 'pending',
            })
            # Haploid calls stay unphased, but retain their one observed sequence comparison
            if not gt:
                row['status'] = 'missing_gt'
                writer.writerow(row)
                continue
            if len(gt) == 1:
                if gt[0] is None:
                    row['status'] = 'missing_gt'
                    writer.writerow(row)
                    continue
                if gt[0] < 0 or gt[0] > len(rec.alts or []):
                    row['status'] = 'invalid_gt'
                    writer.writerow(row)
                    continue
                trgt_hap = rec.ref if gt[0] == 0 else rec.alts[gt[0] - 1]
                row['replace_hap1_seq'] = trgt_hap
                assignment = sample_to_base.get(sample)
                if assignment is None:
                    row['status'] = 'base_sample_missing'
                    writer.writerow(row)
                    continue
                base_index, base_sample = assignment
                base_handle = base_handles[base_index]
                base_contig = contig_for(base_handle, rec.contig)
                if base_contig is None:
                    row['status'] = 'base_contig_missing'
                    writer.writerow(row)
                    continue
                base_hap, base_haplotype, status = reconstruct_available_haploid(
                    base_handle, base_contig, base_sample, rec)
                if status != 'ok':
                    row['status'] = audit_status(status)
                    writer.writerow(row)
                    continue
                if base_haplotype == 1:
                    row['replace_hap1_seq'] = '.'
                row[f'replace_hap{base_haplotype + 1}_seq'] = trgt_hap
                row[f'base_hap{base_haplotype + 1}_seq'] = base_hap
                direct = distance(trgt_hap, base_hap)
                components = ['.', '.']
                components[base_haplotype] = str(direct)
                row['edit_dist_aligned'] = f'{direct} ({components[0]}, {components[1]})'
                row['status'] = 'haploid'
                writer.writerow(row)
                continue
            # Construct both sequences for complete diploid calls; only non-ref heterozygotes get an orientation
            if len(gt) != 2 or any(allele is None for allele in gt):
                row['status'] = 'non_diploid_gt'
                writer.writerow(row)
                continue
            if any(allele < 0 or allele > len(rec.alts or []) for allele in gt):
                row['status'] = 'invalid_gt'
                writer.writerow(row)
                continue
            phase_eligible = gt[0] != gt[1] and any(allele > 0 for allele in gt)
            trgt_haps = [rec.ref if allele == 0 else rec.alts[allele - 1] for allele in gt]
            row['replace_hap1_seq'], row['replace_hap2_seq'] = trgt_haps
            assignment = sample_to_base.get(sample)
            if assignment is None:
                row['status'] = 'base_sample_missing'
                writer.writerow(row)
                continue
            base_index, base_sample = assignment
            base_handle = base_handles[base_index]
            base_contig = contig_for(base_handle, rec.contig)
            if base_contig is None:
                row['status'] = 'base_contig_missing'
                writer.writerow(row)
                continue
            base_hap_1, base_hap_2, status = reconstruct_haplotypes(
                base_handle, base_contig, base_sample, rec)
            if status != 'ok':
                row['status'] = audit_status(status)
                writer.writerow(row)
                continue
            row['base_hap1_seq'] = base_hap_1
            row['base_hap2_seq'] = base_hap_2
            # Aligned compares same haplotype indices; unaligned compares crossed indices
            aligned_1, aligned_2 = distance(base_hap_1, trgt_haps[0]), distance(base_hap_2, trgt_haps[1])
            unaligned_1, unaligned_2 = distance(base_hap_1, trgt_haps[1]), distance(base_hap_2, trgt_haps[0])
            aligned, unaligned = aligned_1 + aligned_2, unaligned_1 + unaligned_2
            row.update({
                'edit_dist_aligned': format_distance(aligned, aligned_1, aligned_2),
                'edit_dist_unaligned': format_distance(unaligned, unaligned_1, unaligned_2),
            })
            if not phase_eligible:
                row['status'] = 'not_het'
            elif aligned == unaligned:
                row['status'] = 'tie'
            elif aligned < unaligned:
                row['edit_dist_pct'] = format_distance_pct(
                    aligned_1, aligned_2, base_hap_1, trgt_haps[0], base_hap_2, trgt_haps[1])
                aligned_pct = combined_edit_distance_pct(
                    aligned_1, aligned_2, base_hap_1, trgt_haps[0], base_hap_2, trgt_haps[1])
                # Gate only summed distance and weighted percentage of winning orientation
                if aligned > ~{max_phase_edit_distance}:
                    row['status'] = 'edit_dist_too_high'
                elif aligned_pct > ~{max_phase_edit_distance_pct}:
                    row['status'] = 'edit_dist_pct_too_high'
                else:
                    call['GT'] = gt
                    call.phased = True
                    call['PS'] = rec.pos
                    row['final_gt'] = gt_string(gt, True)
                    row['status'] = 'phased_aligned'
                    any_phased = True
            else:
                row['edit_dist_pct'] = format_distance_pct(
                    unaligned_1, unaligned_2, base_hap_1, trgt_haps[1], base_hap_2, trgt_haps[0])
                unaligned_pct = combined_edit_distance_pct(
                    unaligned_1, unaligned_2, base_hap_1, trgt_haps[1], base_hap_2, trgt_haps[0])
                # Gate only summed distance and weighted percentage of winning orientation
                if unaligned > ~{max_phase_edit_distance}:
                    row['status'] = 'edit_dist_too_high'
                elif unaligned_pct > ~{max_phase_edit_distance_pct}:
                    row['status'] = 'edit_dist_pct_too_high'
                else:
                    call['GT'] = (gt[1], gt[0])
                    call.phased = True
                    call['PS'] = rec.pos
                    row['final_gt'] = gt_string((gt[1], gt[0]), True)
                    row['status'] = 'phased_unaligned'
                    any_phased = True
            writer.writerow(row)
        if any_phased:
            rec.info['POSTHOC_BACKBONE_PHASED'] = True
        output.write(rec)

source.close()
for handle in base_handles:
    handle.close()
PY
        tabix -f -p vcf ~{prefix}.vcf.gz
    >>>

    output {
        File phased_vcf = "~{prefix}.vcf.gz"
        File phased_vcf_idx = "~{prefix}.vcf.gz.tbi"
        File trv_phasing_summary_tsv = "~{prefix}.trv_phasing_summary.tsv"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 2,
        mem_gb: 12,
        disk_gb: 2 * ceil(size(vcf, "GB") + size(original_vcf, "GB") + size(base_vcfs, "GB")) + 20,
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

# Replace selected TRVs, rebuild envelope relationships, and produce final VCFs and audit
task ApplyTRLocusUpdates {
    input {
        File vcf
        File vcf_idx
        File? replacement_vcf
        File? replacement_vcf_idx
        File? replacement_map_tsv
        File? input_trv_phasing_summary_tsv
        File input_trv_catalog_match_tsv
        Boolean replace_gnomad_str
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        # Stream the contig twice to rebuild TR envelope annotations: collect final TRV intervals, then merge and write
        # When requested, derive gnomAD_STR assignments solely from the catalog report
        python3 <<'PY'
import os
import pysam

def record_key(rec):  # noqa: E302
    return rec.id if rec.id and rec.id != '.' else f'{rec.chrom}:{rec.pos}:{rec.ref}:{",".join(rec.alts or [])}'

def vals(value):  # noqa: E302
    return [str(item) for item in value] if isinstance(value, tuple) else ([] if value is None else [str(value)])

def record_end(rec):  # noqa: E302
    return rec.stop if rec.stop is not None else rec.pos + len(rec.ref) - 1

def normalize_contig(value):  # noqa: E302
    return value[3:] if value.lower().startswith('chr') else value

def parse_explorer(value):  # noqa: E302
    fields = value.rsplit('-', 3)
    if len(fields) != 4:
        return None
    try:
        return fields[0], int(fields[1]), int(fields[2])
    except ValueError:
        return None

def overlaps(rec, start, stop):  # noqa: E302
    return rec.pos <= stop and record_end(rec) >= start

def inclusive_overlap(rec, start, stop):  # noqa: E302
    return max(0, min(record_end(rec), stop) - max(rec.pos, start) + 1)

def trid_text(rec):  # noqa: E302
    return ','.join(vals(rec.info.get('TRID')))

def truth(value):  # noqa: E302
    return value.strip().lower() == 'true'

def signature(rec):  # noqa: E302
    return (rec.contig, rec.pos, record_end(rec), rec.ref,
            ','.join(rec.alts or []), rec.id or '.')

def matching_candidates(records, explorers):  # noqa: E302
    candidates = []
    for rec in records:
        if rec.info.get('allele_type') != 'trv':
            continue
        trid = trid_text(rec)
        for explorer, start, stop in explorers:
            if explorer in trid and overlaps(rec, start, stop):
                candidates.append((inclusive_overlap(rec, start, stop), rec))
                break
    return sorted(candidates, key=lambda item: (-item[0], signature(item[1]), trid_text(item[1])))

# Drive gnomAD_STR selection from the catalog report to avoid reinterpreting the JSON during output assembly
catalog_rows = []
with open('~{input_trv_catalog_match_tsv}') as handle:
    header = next(handle, '').rstrip('\n').split('\t')
    expected = [
        'locus_id', 'TRExplorerV1', 'has_diseases', 'input_overlapping_trids',
        'input_has_TRExplorerV1_substring', 'trgt_overlapping_trids',
        'trgt_has_TRExplorerV1_substring', 'trgt_matching_allele_count',
    ]
    # Accept the canonical columns as a prefix so cohort variants may append extra trailing columns
    if header[:len(expected)] != expected:
        raise RuntimeError('Unexpected trv_catalog_match_tsv header')
    for line_number, line in enumerate(handle, start=2):
        fields = line.rstrip('\n').split('\t')
        if len(fields) < len(expected):
            raise RuntimeError(f'Malformed catalog report row {line_number}')
        if not truth(fields[2]):
            continue
        explorers = []
        for value in fields[1].split(','):
            parsed = parse_explorer(value)
            if parsed:
                explorers.append((value, parsed[1], parsed[2]))
        if not explorers:
            raise RuntimeError(f'Eligible catalog row {line_number} has no parseable TRExplorerV1 value')
        catalog_rows.append({
            'line_number': line_number,
            'locus_id': fields[0],
            'explorers': explorers,
            'main_strict': truth(fields[4]),
            'trgt_strict': truth(fields[6]),
            # Input matches leave TRGT cells blank; numeric AC gates recovery for rows needing a TRGT lookup
            'trgt_ac': int(fields[7]) if fields[7] not in ('', '.') else 0,
        })

replacement_path = '~{replacement_vcf}'
map_path = '~{replacement_map_tsv}'
replace_gnomad_str = '~{replace_gnomad_str}'.lower() == 'true'
replacements = []

def recompute_ac(rec):  # noqa: E302
    """Recalculate allele-specific AC after all replacement transformations."""
    counts = [0] * len(rec.alts or [])
    for sample_data in rec.samples.values():
        for allele in sample_data.get('GT') or ():
            if allele is not None and 0 < allele <= len(counts):
                counts[allele - 1] += 1
    rec.info['AC'] = tuple(counts)
    return sum(counts)

if replacement_path and os.path.exists(replacement_path):
    with pysam.VariantFile(replacement_path) as handle:
        if 'AC' not in handle.header.info:
            handle.header.add_meta(
                'INFO',
                items=[('ID', 'AC'), ('Number', 'A'), ('Type', 'Integer'), ('Description', 'Number of alleles observed')],
            )
        # Final safety gate: keep INFO/AC synchronized with GT so an AC=0 call cannot remove its mapped input TRV
        for rec in handle:
            rec = rec.copy()
            if recompute_ac(rec) >= 1:
                replacements.append(rec)

retained_replacement_keys = {record_key(rec) for rec in replacements}
map_rows_by_new = {}
if map_path and os.path.exists(map_path):
    with open(map_path) as handle:
        next(handle, None)
        for line in handle:
            fields = line.rstrip('\n').split('\t')
            if len(fields) >= 4 and fields[3] in ('replace', 'add') and fields[0] in retained_replacement_keys:
                map_rows_by_new.setdefault(fields[0], []).append(fields)

# Every final replacement maps to exactly one map row: 'replace' displaces an original TRV, 'add' inserts a novel locus
# Map rows for AC=0 records are deliberately ignored, preserving the original call
for new_key in retained_replacement_keys:
    rows = map_rows_by_new.get(new_key, [])
    if len(rows) != 1:
        raise RuntimeError(f'Retained replacement {new_key} requires exactly one map row')
    status_field, old_field = rows[0][3], rows[0][1]
    if status_field == 'replace' and old_field == '.':
        raise RuntimeError(f'Replacement {new_key} status=replace requires a displaced TRV')
    if status_field == 'add' and old_field != '.':
        raise RuntimeError(f'Insertion {new_key} status=add must not displace a TRV')
replace_old = {rows[0][1] for rows in map_rows_by_new.values() if rows[0][3] == 'replace'}

base = pysam.VariantFile('~{vcf}', index_filename='~{vcf_idx}')
header = base.header.copy()

def vcf_contig_for(header, catalog_contig):  # noqa: E302
    for name in header.contigs:
        if normalize_contig(name) == normalize_contig(catalog_contig):
            return name
    raise RuntimeError(f'Catalog contig {catalog_contig} absent from input VCF')

# Select exactly one output target per assignable catalog row, with main-VCF rows taking precedence over raw-TRGT rows
gnomad_assignments = {}
if replace_gnomad_str:
    for row in catalog_rows:
        if row['main_strict']:
            main_records = []
            for _, start, stop in row['explorers']:
                contig = vcf_contig_for(base.header, row['explorers'][0][0].rsplit('-', 3)[0])
                main_records.extend(
                    rec.copy() for rec in base.fetch(contig, max(0, start - 1), stop + 1)
                    if record_key(rec) not in replace_old
                )
            candidates = matching_candidates(main_records, row['explorers'])
            target_kind = 'main VCF'
        elif row['trgt_strict'] and row['trgt_ac']:
            candidates = matching_candidates(replacements, row['explorers'])
            target_kind = 'replacement VCF'
        else:
            continue
        # A raw TRGT match dropped by final AC filtering emits no replacement, leaving its original VCF record in place
        if not candidates and target_kind == 'replacement VCF':
            continue
        if not candidates:
            raise RuntimeError(
                f'Catalog row {row["line_number"]} ({row["locus_id"]}) has no matching {target_kind} record'
            )
        target = signature(candidates[0][1])
        gnomad_assignments.setdefault(target, set()).add(row['locus_id'])
base.close()

# Pass A: the final TRV interval set is the non-displaced input TRVs plus the retained replacements
# Collect it up front because a non-TRV can be enveloped by a TRV at the same POS that streams in later
intervals = []
with pysam.VariantFile('~{vcf}') as base_scan:
    for rec in base_scan:
        if rec.info.get('allele_type') == 'trv' and record_key(rec) not in replace_old:
            intervals.append((rec.contig, rec.pos, record_end(rec), record_key(rec)))
for rec in replacements:
    if rec.info.get('allele_type') == 'trv':
        intervals.append((rec.contig, rec.pos, record_end(rec), record_key(rec)))
intervals.sort()
by_contig = {}
for interval in intervals:
    by_contig.setdefault(interval[0], []).append(interval)

if replacements:
    # Preserve all annotations generated on the small replacement VCF
    with pysam.VariantFile(replacement_path) as repl_header_source:
        header.merge(repl_header_source.header)
for name, number, type_, description in [
        ('TR_ENVELOPED', '0', 'Flag', 'Variant enveloped by tandem repeat'),
        ('TRID', '1', 'String', 'ID of enveloping tandem repeat'),
        ('allele_type', '1', 'String', 'Allele type'),
        ('SOURCE', '1', 'String', 'Source of variant call'),
        ('allele_length', '1', 'Integer', 'Allele length'),
        ('POSTHOC_BACKBONE_PHASED', '0', 'Flag', 'At least one heterozygous call phased by base-VCF sequence agreement'),
]:
    if name not in header.info:
        header.add_meta(
            'INFO',
            items=[('ID', name), ('Number', number), ('Type', type_), ('Description', description)],
        )
if replace_gnomad_str:
    # Multiple disease-associated catalog loci can intentionally select one TRV
    # Rebuild rather than mutate in place: removing a definition shifts htslib field IDs and corrupts translated fields
    rebuilt_header = pysam.VariantHeader()
    for line in str(header).splitlines():
        if not line.startswith('##'):
            continue
        if (line.startswith('##fileformat=') or line.startswith('##FILTER=<ID=PASS,')
                or line.startswith('##INFO=<ID=gnomAD_STR,')):
            continue
        rebuilt_header.add_line(line)
    for sample in header.samples:
        rebuilt_header.add_sample(sample)
    header = rebuilt_header
    header.add_meta(
        'INFO',
        items=[('ID', 'gnomAD_STR'), ('Number', '.'), ('Type', 'String'),
               ('Description', 'Matched gnomAD tandem-repeat locus ID')],
    )
if 'PS' not in header.formats:
    header.add_meta(
        'FORMAT',
        items=[('ID', 'PS'), ('Number', '1'), ('Type', 'Integer'), ('Description', 'Phase set')],
    )

contig_order = {name: i for i, name in enumerate(header.contigs)}
replacements.sort(key=lambda rec: (contig_order[rec.contig], rec.pos, record_end(rec)))


def merged_records():  # noqa: E302
    """Yield final records in output order: input records minus displaced TRVs,
    interleaved with the position-sorted replacements."""
    replacement_index = 0
    with pysam.VariantFile('~{vcf}') as base_stream:
        for rec in base_stream:
            while replacement_index < len(replacements) and (
                contig_order[replacements[replacement_index].contig] < contig_order[rec.contig]
                or (
                    replacements[replacement_index].contig == rec.contig
                    and replacements[replacement_index].pos <= rec.pos
                )
            ):
                yield replacements[replacement_index].copy()
                replacement_index += 1
            if record_key(rec) not in replace_old:
                yield rec
    while replacement_index < len(replacements):
        yield replacements[replacement_index].copy()
        replacement_index += 1


out = pysam.VariantFile('~{prefix}.trv_postprocessed.vcf.gz', 'wz', header=header)
emitted_gnomad_targets = set()
with out:
    current_contig = None
    contig_intervals = []
    interval_index = 0
    active = []
    for rec in merged_records():
        rec.translate(header)
        if rec.contig != current_contig:
            current_contig = rec.contig
            contig_intervals = by_contig.get(current_contig, [])
            interval_index = 0
            active = []
        while interval_index < len(contig_intervals) and contig_intervals[interval_index][1] <= rec.pos:
            active.append(contig_intervals[interval_index])
            interval_index += 1
        active = [interval for interval in active if interval[2] >= rec.pos]
        if 'TR_ENVELOPED' in rec.info:
            del rec.info['TR_ENVELOPED']
            if 'TRID' in rec.info:
                del rec.info['TRID']
        if replace_gnomad_str and 'gnomAD_STR' in rec.info:
            del rec.info['gnomAD_STR']
        if rec.info.get('allele_type') != 'trv':
            for _, start, stop, trid in active:
                if rec.pos >= start and record_end(rec) <= stop:
                    rec.info['TR_ENVELOPED'] = True
                    rec.info['TRID'] = trid
                    break
        target = signature(rec)
        if replace_gnomad_str and target in gnomad_assignments:
            if target in emitted_gnomad_targets:
                raise RuntimeError(f'gnomAD_STR target was emitted more than once: {target}')
            rec.info['gnomAD_STR'] = tuple(sorted(gnomad_assignments[target]))
            emitted_gnomad_targets.add(target)
        out.write(rec)

if replace_gnomad_str and emitted_gnomad_targets != set(gnomad_assignments):
    missing = sorted(set(gnomad_assignments) - emitted_gnomad_targets)
    raise RuntimeError(f'gnomAD_STR targets were not emitted: {missing}')

with open('~{prefix}.trv_catalog_match.tsv', 'w') as out, open('~{input_trv_catalog_match_tsv}') as source:
    out.write(source.read())

# Always materialize phase audit so no-replacement contigs have stable workflow output
phase_audit_path = '~{input_trv_phasing_summary_tsv}'
phase_header = (
    'base_trid\treplace_trid\tsample_id\tbase_gt\treplace_gt'
    '\tbase_hap1_seq\tbase_hap2_seq\treplace_hap1_seq\treplace_hap2_seq'
    '\tedit_dist_aligned\tedit_dist_unaligned\tedit_dist_pct'
    '\tmax_edit_dist\tmax_edit_dist_pct\tfinal_gt\tstatus\n'
)
with open('~{prefix}.trv_phasing_summary.tsv', 'w') as out:
    if phase_audit_path and os.path.exists(phase_audit_path):
        with open(phase_audit_path) as source:
            out.write(source.read())
    else:
        out.write(phase_header)
PY
        tabix -f -p vcf ~{prefix}.trv_postprocessed.vcf.gz
        bcftools view -i 'INFO/allele_type="trv"' -Oz -o ~{prefix}.trv_subsetted.vcf.gz ~{prefix}.trv_postprocessed.vcf.gz
        tabix -f -p vcf ~{prefix}.trv_subsetted.vcf.gz
    >>>
    output {
        File trv_postprocessed_vcf = "~{prefix}.trv_postprocessed.vcf.gz"
        File trv_postprocessed_vcf_idx = "~{prefix}.trv_postprocessed.vcf.gz.tbi"
        File trv_subsetted_vcf = "~{prefix}.trv_subsetted.vcf.gz"
        File trv_subsetted_vcf_idx = "~{prefix}.trv_subsetted.vcf.gz.tbi"
        File trv_catalog_match_tsv = "~{prefix}.trv_catalog_match.tsv"
        File trv_phasing_summary_tsv = "~{prefix}.trv_phasing_summary.tsv"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 2,
        mem_gb: 12,
        disk_gb: 3 * ceil(size(vcf, "GB")) + 20,
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
