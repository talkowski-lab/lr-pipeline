version 1.0

import "../annotation/AnnotateInSilicoPredictors.wdl"
import "../annotation/AnnotateRegion.wdl"
import "../annotation/AnnotateVRS.wdl"
import "../tools/MergeHiPhaseCallsets.wdl"
import "../utils/Structs.wdl"
import "AnnotateVcf.wdl"

workflow PostprocessTRLoci {
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

        File? swap_samples_base
        Int min_phase_edit_distance_delta = 10
        Float min_phase_similarity = 0.90
        Boolean run_decrement_trv_ids

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
        RuntimeAttr? runtime_attr_vrs_annotate
        RuntimeAttr? runtime_attr_vrs_extract
        RuntimeAttr? runtime_attr_region_annotate
        RuntimeAttr? runtime_attr_insilico_annotate
        RuntimeAttr? runtime_attr_attach_annotations
        RuntimeAttr? runtime_attr_apply
    }

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

            call MergeHiPhaseCallsets.AddTRGTEndTag {
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
                match_keys = DiscoverTRLoci.trgt_match_keys,
                gnomad_tr_json = gnomad_tr_json,
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
                    run_decrement_trv_ids = run_decrement_trv_ids,
                    prefix = "~{prefix}.~{contig}.replacement_seed",
                    docker = utils_docker,
                    runtime_attr_override = runtime_attr_prepare
            }

            call PhaseReplacementLoci {
                input:
                    vcf = PrepareReplacementLoci.prepared_vcf,
                    vcf_idx = PrepareReplacementLoci.prepared_vcf_idx,
                    base_vcfs = base_vcfs,
                    base_vcf_idxs = base_vcf_idxs,
                    swap_samples_base = swap_samples_base,
                    min_phase_edit_distance_delta = min_phase_edit_distance_delta,
                    min_phase_similarity = min_phase_similarity,
                    prefix = "~{prefix}.~{contig}.replacement_seed.trv_phasing",
                    docker = utils_docker,
                    runtime_attr_override = runtime_attr_phase_replacements
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
                    annotations_tsvs = [ExtractReplacementVRS.annotations_tsv, AnnotateReplacementRegion.annotations_tsv, AnnotateReplacementInSilico.annotations_tsv],
                    prefix = "~{prefix}.~{contig}.replacement_annotated",
                    info_names = [["VRS_Allele_IDs", "VRS_Error", "VRS_Starts", "VRS_Ends", "VRS_States", "VRS_Lengths", "VRS_RepeatSubunitLengths"], ["REGION"], ["cadd_raw_score", "cadd_phred", "pangolin_largest", "revel_max", "phylop", "spliceai_ds_max"]],
                    info_descriptions = [["VRS allele identifiers", "VRS annotation error", "VRS start positions", "VRS end positions", "VRS allele states", "VRS allele lengths", "VRS repeat subunit lengths"], ["Genomic context of variant"], ["CADD raw score", "CADD PHRED score", "Largest Pangolin delta score", "Maximum REVEL score", "PhyloP score", "Maximum SpliceAI delta score"]],
                    info_types = [["String", "String", "Integer", "Integer", "String", "String", "String"], ["String"], ["Float", "Float", "Float", "Float", "Float", "Float"]],
                    info_numbers = [["R", ".", "R", "R", ".", ".", "."], ["1"], ["1", "1", "1", "1", "1", "1"]],
                    subset_vcf_strings = [],
                    awk_tsv_conditions = [],
                    subset_tsv_columns = [[6, 7, 8, 9, 10, 11, 12], [6], [6, 7, 8, 9, 10, 11]],
                    strip_info_fields_per_tsv = [false, false, false],
                    docker = utils_docker,
                    runtime_attr_override = runtime_attr_attach_annotations
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
            merged_status_tsv = KeepMergedTRGTWithAC.status_tsv,
            gnomad_tr_json = gnomad_tr_json,
            prefix = "~{prefix}.~{contig}",
            docker = utils_docker,
            runtime_attr_override = runtime_attr_apply
    }

    output {
        File trv_postprocessed_vcf = ApplyTRLocusUpdates.trv_postprocessed_vcf
        File trv_postprocessed_vcf_idx = ApplyTRLocusUpdates.trv_postprocessed_vcf_idx
        File trv_updated_vcf = ApplyTRLocusUpdates.trv_updated_vcf
        File trv_updated_vcf_idx = ApplyTRLocusUpdates.trv_updated_vcf_idx
        File trv_catalog_match_tsv = ApplyTRLocusUpdates.trv_catalog_match_tsv
        File trv_phasing_summary_tsv = ApplyTRLocusUpdates.trv_phasing_summary_tsv
    }
}

# Find disease-associated catalog loci in the main VCF, then in per-sample TRGT VCFs.
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

        # Validate sample/file alignment before using indexed windows for strict TRID matching.
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

def record_key(rec):  # noqa: E302
    return rec.id if rec.id and rec.id != '.' else f'{rec.chrom}:{rec.pos}:{rec.ref}:{",".join(rec.alts or [])}'

trgt_paths = [line.rstrip('\n') for line in open('~{write_lines(trgt_vcfs)}') if line.strip()]  # noqa: E305
trgt_indexes = [line.rstrip('\n') for line in open('~{write_lines(trgt_vcf_idxs)}') if line.strip()]
sample_ids = [line.rstrip('\n') for line in open('~{write_lines(sample_ids)}') if line.strip()]
if not trgt_paths or len(trgt_paths) != len(trgt_indexes) or len(trgt_paths) != len(sample_ids):
    raise RuntimeError('trgt_vcfs, trgt_vcf_idxs, and sample_ids must be non-empty parallel arrays')
if len(set(sample_ids)) != len(sample_ids):
    raise RuntimeError('sample_ids must be unique')
if not os.path.exists('~{vcf_idx}') or not all(os.path.exists(path) for path in trgt_indexes):
    raise RuntimeError('All VCF indexes must be localized')

with open('~{gnomad_tr_json}') as handle:
    catalog = json.load(handle)
loci = []
for entry in catalog:
    # Catalog key is capitalized. Missing, non-array, and empty Diseases values are ignored.
    diseases = entry.get('Diseases') if entry else None
    if not isinstance(diseases, list) or not diseases:
        continue
    if not entry or not entry.get('LocusId'):
        continue
    explorers = entry.get('TRExplorerV1')
    explorers = explorers if isinstance(explorers, list) else [explorers]
    for explorer in explorers:
        parsed = parse_explorer(str(explorer)) if explorer else None
        if parsed and normalize_contig(parsed[0]) == normalize_contig('~{contig}'):
            loci.append((str(entry['LocusId']), str(explorer), parsed[1], parsed[2]))

with pysam.VariantFile('~{vcf}') as main:
    if list(main.header.samples) != sample_ids:
        raise RuntimeError('sample_ids must exactly match main VCF sample order')
    main_matches = {explorer: [] for _, explorer, _, _ in loci}
    for _, explorer, start, stop in loci:
        for rec in main.fetch('~{contig}', max(0, start - 1), stop + 1):
            if rec.info.get('allele_type') != 'trv':
                continue
            if explorer in '|'.join(values(rec.info.get('TRID'))):
                main_matches[explorer].append(record_key(rec))

for i, path in enumerate(trgt_paths):
    with pysam.VariantFile(path) as trgt:
        if list(trgt.header.samples) != [sample_ids[i]]:
            raise RuntimeError(f'TRGT VCF {path} must contain only sample {sample_ids[i]}')

trgt_matches = {explorer: [] for _, explorer, _, _ in loci}
for path in trgt_paths:
    with pysam.VariantFile(path) as trgt:
        for _, explorer, start, stop in loci:
            if main_matches[explorer]:
                continue
            for rec in trgt.fetch('~{contig}', max(0, start - 1), stop + 1):
                if explorer in '|'.join(values(rec.info.get('TRID'))):
                    trgt_matches[explorer].append(record_key(rec))

with open('~{prefix}.trv_catalog_match.tsv', 'w') as out:
    out.write('source\tlocus_id\tTRExplorerV1\tmatching_vcf_ids\tstatus\n')
    for locus_id, explorer, _, _ in loci:
        if main_matches[explorer]:
            out.write(f'main_vcf\t{locus_id}\t{explorer}\t{",".join(sorted(set(main_matches[explorer])))}\tmain_match\n')
        elif trgt_matches[explorer]:
            out.write(f'trgt_vcf\t{locus_id}\t{explorer}\t{",".join(sorted(set(trgt_matches[explorer])))}\ttrgt_candidate\n')
        else:
            out.write(f'none\t{locus_id}\t{explorer}\t.\tno_match\n')
with open('~{prefix}.trgt_match_keys.txt', 'w') as out:
    for _, explorer, _, _ in loci:
        if not main_matches[explorer] and trgt_matches[explorer]:
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

# Retain only disease-associated fallback loci selected during discovery from one sample VCF.
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

        # Fetch narrow catalog intervals and confirm strict TRID matches before writing records.
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

# Recompute cohort allele counts, retain AC-positive merged loci, and audit AC-zero drops.
task KeepMergedTRGTWithAC {
    input {
        File vcf
        File vcf_idx
        File match_keys
        File gnomad_tr_json
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        # Rebuild selected disease-locus metadata and calculate AC directly from merged genotypes.
        python3 <<'PY'
import json
import pysam

def record_key(rec):  # noqa: E302
    return rec.id if rec.id and rec.id != '.' else f'{rec.chrom}:{rec.pos}:{rec.ref}:{",".join(rec.alts or [])}'

def vals(value):  # noqa: E302
    return [str(item) for item in value] if isinstance(value, tuple) else ([] if value is None else [str(value)])

def normalize_contig(value):  # noqa: E302
    return value[3:] if value.lower().startswith('chr') else value

def parse_explorer(value):  # noqa: E302
    fields = value.rsplit('-', 3)
    return fields[0], int(fields[1]), int(fields[2])

with open('~{gnomad_tr_json}') as handle:  # noqa: E305
    catalog = json.load(handle)
selected_keys = {line.strip() for line in open('~{match_keys}') if line.strip()}
locus_ids = {}
for entry in catalog:
    # Apply the same non-empty Diseases-array rule used during initial discovery.
    diseases = entry.get('Diseases') if entry else None
    if not isinstance(diseases, list) or not diseases:
        continue
    explorer = entry.get('TRExplorerV1') if entry else None
    explorers = explorer if isinstance(explorer, list) else [explorer]
    for value in explorers:
        if entry and entry.get('LocusId') and value and str(value) in selected_keys:
            locus_ids.setdefault(str(value), set()).add(str(entry['LocusId']))
loci = [(locus_id, explorer, *parse_explorer(explorer)) for explorer, ids in locus_ids.items() for locus_id in ids]

src = pysam.VariantFile('~{vcf}')
if 'AC' not in src.header.info:
    src.header.add_meta(
        'INFO',
        items=[('ID', 'AC'), ('Number', 'A'), ('Type', 'Integer'), ('Description', 'Number of alleles observed')],
    )
out = pysam.VariantFile('~{prefix}.vcf.gz', 'wz', header=src.header)
status = open('~{prefix}.status.tsv', 'w')
status.write('source\tlocus_id\tTRExplorerV1\tmatching_vcf_ids\tstatus\n')
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
    searchable = '|'.join(vals(rec.info.get('TRID')))
    found = [(locus, explorer) for locus, explorer, _, _, _ in loci if explorer in searchable]
    if not found:
        found = [
            (locus, explorer)
            for locus, explorer, locus_contig, start, stop in loci
            if normalize_contig(locus_contig) == normalize_contig(rec.contig)
            and max(rec.pos, start) < min(rec.stop, stop)
        ]
    if not found:
        raise RuntimeError(f'Merged TRGT record {record_key(rec)} could not be attributed to a selected catalog locus')
    if ac >= 1:
        rec.info['AC'] = tuple(counts)
        out.write(rec)
        kept += 1
        for locus, explorer in found:
            status.write(f'trgt_merged\t{locus}\t{explorer}\t{record_key(rec)}\ttrgt_merged_retained_AC={ac}\n')
    else:
        for locus, explorer in found:
            status.write(f'trgt_merged\t{locus}\t{explorer}\t{record_key(rec)}\tdropped_AC0\n')
src.close()
out.close()
status.close()
PY
        tabix -f -p vcf ~{prefix}.vcf.gz
        bcftools view -H ~{prefix}.vcf.gz | wc -l
    >>>
    output {
        File retained_vcf = "~{prefix}.vcf.gz"
        File retained_vcf_idx = "~{prefix}.vcf.gz.tbi"
        File status_tsv = "~{prefix}.status.tsv"
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

# Reconcile replacement headers and select old TRVs by positive overlap.
task PrepareReplacementLoci {
    input {
        File vcf
        File vcf_idx
        File replacement_vcf
        File replacement_vcf_idx
        Array[String] sample_ids
        Boolean run_decrement_trv_ids
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        # Match FORMAT/AL to the main VCF before records move between pysam headers.
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

        # Require unique positive-overlap replacements; phase after all replacements are prepared.
        python3 <<'PY'
import pysam

def record_key(rec):  # noqa: E302
    return rec.id if rec.id and rec.id != '.' else f'{rec.chrom}:{rec.pos}:{rec.ref}:{",".join(rec.alts or [])}'

def record_end(rec):  # noqa: E302
    return rec.stop if rec.stop is not None else rec.pos + len(rec.ref) - 1

def overlap(left, right):  # noqa: E302
    return max(0, min(record_end(left), record_end(right)) - max(left.pos, right.pos) + 1)

def decrement_trv_id(trv_id):  # noqa: E302
    """Match PostprocessCallset decrementing for only newly introduced TRVs."""
    head, _, tail = trv_id.rpartition('-')
    if not head.endswith('-TRV'):
        return trv_id
    try:
        return f'{head}-{int(tail) - 1}'
    except ValueError:
        return trv_id

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
    # AnnotateRegion requires declaration, then derives TRV length from REF.
    incoming.header.add_meta(
        'INFO',
        items=[('ID', 'allele_length'), ('Number', '1'), ('Type', 'Integer'), ('Description', 'Allele length')],
    )
if 'PS' not in incoming.header.formats:
    incoming.header.add_meta(
        'FORMAT',
        items=[('ID', 'PS'), ('Number', '1'), ('Type', 'Integer'), ('Description', 'Phase set')],
    )
header = incoming.header.copy()
out = pysam.VariantFile('~{prefix}.vcf.gz', 'wz', header=header)
mapping = open('~{prefix}.map.tsv', 'w')
mapping.write('new_record_id\told_record_id\toverlap_bp\tstatus\tphase_summary\n')

samples = [line.rstrip('\n') for line in open('~{write_lines(sample_ids)}') if line.strip()]
if list(base.header.samples) != samples or list(header.samples) != samples:
    raise RuntimeError('Main, merged TRGT, and sample_ids sample order must match exactly')
used_old_records = set()
id_seen = {}
run_decrement_trv_ids = ~{true="True" false="False" run_decrement_trv_ids}
for rec in incoming:
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
    # Match IntegrateTRs naming so enveloped records can point to a stable replacement ID.
    # Decrement new replacement IDs before suffixing; existing main-VCF IDs remain untouched.
    new_id = f'{rec.chrom}-{rec.pos}-TRV-{len(rec.ref) - 1}'
    if run_decrement_trv_ids:
        new_id = decrement_trv_id(new_id)
    id_seen[new_id] = id_seen.get(new_id, 0) + 1
    rec.id = new_id if id_seen[new_id] == 1 else f'{new_id}_{id_seen[new_id]}'
    rec.info['allele_type'] = 'trv'
    rec.info['SOURCE'] = 'TRExplorer'
    mapping.write(f'{record_key(rec)}\t{old_key}\t{best_overlap}\treplace\tsee_trv_phasing_summary_tsv\n')
    out.write(rec)
base.close()
incoming.close()
out.close()
mapping.close()
PY
        rm -f replacement.header replacement.normalized.vcf
        tabix -f -p vcf ~{prefix}.vcf.gz
    >>>
    output {
        File prepared_vcf = "~{prefix}.vcf.gz"
        File prepared_vcf_idx = "~{prefix}.vcf.gz.tbi"
        File replacement_map_tsv = "~{prefix}.map.tsv"
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

# Phase eligible replacement calls by sequence agreement with matching base-VCF haplotypes.
task PhaseReplacementLoci {
    input {
        File vcf
        File vcf_idx
        Array[File] base_vcfs
        Array[File] base_vcf_idxs
        File? swap_samples_base
        Int min_phase_edit_distance_delta
        Float min_phase_similarity
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        # Fetch only each replacement interval from indexed base VCFs; swap map uses BackbonePhase raw-to-canonical format.
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


def contig_for(handle, contig):  # noqa: E302
    if contig in handle.header.contigs:
        return contig
    alternate = contig[3:] if contig.startswith('chr') else f'chr{contig}'
    return alternate if alternate in handle.header.contigs else None


def normalized_base_gt(call):  # noqa: E302
    """Keep only base calls whose two haplotypes can be unambiguously reconstructed."""
    gt = call.get('GT')
    if not gt or len(gt) != 2:
        return None, 'base_missing_gt', 0
    if gt == (None, None):
        # Sparse truth records encode no call at this site; reconstruct both haplotypes as REF.
        return (0, 0), None, 2
    if any(allele is not None and allele < 0 for allele in gt):
        return None, 'base_invalid_gt', 0
    if not call.phased and (None in gt or gt[0] != gt[1]):
        # 1/. and unphased heterozygotes have unknown haplotype orientation.
        return None, 'base_ambiguous_unphased_gt', 0
    return tuple(0 if allele is None else allele for allele in gt), None, gt.count(None)


def reconstruct_haplotypes(base_handle, contig, sample, rec):  # noqa: E302
    """Apply complete, non-overlapping base variants to replacement REF for both haplotypes."""
    locus_start = rec.pos
    locus_end = record_end(rec)
    variants = []
    missing_as_reference_count = 0
    for base_rec in base_handle.fetch(contig, rec.start, rec.stop):
        call = base_rec.samples[sample]
        gt, status, missing_as_reference = normalized_base_gt(call)
        if status:
            return None, None, status, [], 0
        missing_as_reference_count += missing_as_reference
        # Cohort VCFs contain many overlapping records carried as reference for this sample.
        if gt == (0, 0):
            continue
        base_end = record_end(base_rec)
        if base_rec.pos < locus_start or base_end > locus_end:
            return None, None, 'base_boundary_overlapping_variant', [], 0
        if any(allele > len(base_rec.alts or []) for allele in gt):
            return None, None, 'base_invalid_allele_index', [], 0
        selected_alts = [base_rec.alts[allele - 1] for allele in gt if allele > 0]
        if any(not alt or alt == '*' or alt.startswith('<') or '[' in alt or ']' in alt for alt in selected_alts):
            return None, None, 'base_symbolic_allele', [], 0
        offset = base_rec.pos - locus_start
        if rec.ref[offset:offset + len(base_rec.ref)].upper() != base_rec.ref.upper():
            return None, None, 'base_reference_mismatch', [], 0
        variants.append((offset, base_rec, gt))
    variants.sort(key=lambda item: item[0])
    output = [[], []]
    variant_ids = []
    for haplotype in range(2):
        cursor = 0
        for offset, base_rec, gt in variants:
            allele = gt[haplotype]
            if allele == 0:
                continue
            if offset < cursor:
                return None, None, 'base_overlapping_variants', [], 0
            output[haplotype].append(rec.ref[cursor:offset])
            output[haplotype].append(base_rec.alts[allele - 1])
            cursor = offset + len(base_rec.ref)
        output[haplotype].append(rec.ref[cursor:])
    for _, base_rec, _ in variants:
        variant_ids.append(record_key(base_rec))
    return ''.join(output[0]), ''.join(output[1]), 'ok', variant_ids, missing_as_reference_count


def distance(left, right):  # noqa: E302
    return edlib.align(left.upper(), right.upper(), task='distance')['editDistance']


def similarity(distance_value, left, right):  # noqa: E302
    return 1.0 if not left and not right else 1.0 - distance_value / max(len(left), len(right))


base_paths = [line.strip() for line in open('~{write_lines(base_vcfs)}') if line.strip()]
base_idx_paths = [line.strip() for line in open('~{write_lines(base_vcf_idxs)}') if line.strip()]
if len(base_paths) != len(base_idx_paths):
    raise RuntimeError('base_vcfs and base_vcf_idxs must have equal lengths')
if ~{min_phase_edit_distance_delta} < 0:
    raise RuntimeError('min_phase_edit_distance_delta must be non-negative')
if not 0.0 <= ~{min_phase_similarity} <= 1.0:
    raise RuntimeError('min_phase_similarity must be between 0 and 1')

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

# First matching base VCF wins, mirroring BackbonePhase sample assignment.
sample_to_base = {}
base_handles = []
for index, path in enumerate(base_paths):
    # Use explicit localized index path; Cromwell need not preserve sibling filenames.
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
    'record_id', 'chrom', 'pos', 'end', 'sample_id', 'input_gt',
    'trgt_haplotype_1_sequence', 'trgt_haplotype_2_sequence', 'base_vcf_index', 'base_sample_id',
    'base_variant_ids', 'base_missing_as_reference_count', 'base_haplotype_1_sequence',
    'base_haplotype_2_sequence', 'direct_haplotype_1_edit_distance',
    'direct_haplotype_2_edit_distance', 'direct_edit_distance', 'swapped_haplotype_1_edit_distance',
    'swapped_haplotype_2_edit_distance', 'swapped_edit_distance', 'distance_delta',
    'direct_similarity', 'swapped_similarity', 'winning_orientation', 'winning_phased_gt',
    'winning_edit_distance', 'winning_similarity', 'min_phase_edit_distance_delta',
    'min_phase_similarity', 'passes_edit_distance_delta', 'passes_similarity',
    'phase_set', 'status', 'details',
]
output = pysam.VariantFile('~{prefix}.vcf.gz', 'wz', header=header)
audit = open('~{prefix}.trv_phasing_summary.tsv', 'w')
with output, audit:
    writer = csv.DictWriter(audit, fieldnames=audit_fields, delimiter='\t', lineterminator='\n')
    writer.writeheader()
    for rec in source:
        rec.translate(header)
        any_phased = False
        for sample in header.samples:
            call = rec.samples[sample]
            gt = call.get('GT')
            input_gt = gt_string(gt, call.phased)
            # Replacement calls must be unphased unless sequence evidence below selects an orientation.
            call.phased = False
            if call.get('PS') is not None:
                call['PS'] = None
            # Only complete, diploid, non-reference heterozygotes can receive an orientation.
            if (
                not gt or len(gt) != 2 or any(allele is None for allele in gt)
                or gt[0] == gt[1] or not any(allele > 0 for allele in gt)
            ):
                continue
            row = dict.fromkeys(audit_fields, '.')
            row.update({
                'record_id': record_key(rec), 'chrom': rec.contig, 'pos': rec.pos,
                'end': record_end(rec), 'sample_id': sample, 'input_gt': input_gt,
                'min_phase_edit_distance_delta': ~{min_phase_edit_distance_delta},
                'min_phase_similarity': f'{~{min_phase_similarity}:.6f}',
            })
            if any(allele < 0 or allele > len(rec.alts or []) for allele in gt):
                row['status'] = 'invalid_trgt_allele_index'
                writer.writerow(row)
                continue
            trgt_haps = [rec.ref if allele == 0 else rec.alts[allele - 1] for allele in gt]
            row['trgt_haplotype_1_sequence'], row['trgt_haplotype_2_sequence'] = trgt_haps
            assignment = sample_to_base.get(sample)
            if assignment is None:
                row['status'] = 'no_matching_base_sample'
                writer.writerow(row)
                continue
            base_index, base_sample = assignment
            base_handle = base_handles[base_index]
            base_contig = contig_for(base_handle, rec.contig)
            row.update({'base_vcf_index': base_index, 'base_sample_id': base_sample})
            if base_contig is None:
                row['status'] = 'base_contig_missing'
                writer.writerow(row)
                continue
            base_hap_1, base_hap_2, status, variant_ids, missing_as_reference = reconstruct_haplotypes(
                base_handle, base_contig, base_sample, rec)
            if status != 'ok':
                row['status'] = status
                writer.writerow(row)
                continue
            row['base_haplotype_1_sequence'] = base_hap_1
            row['base_haplotype_2_sequence'] = base_hap_2
            row['base_variant_ids'] = ','.join(variant_ids) if variant_ids else '.'
            row['base_missing_as_reference_count'] = missing_as_reference
            row['details'] = f'base_missing_as_reference_count={missing_as_reference}'
            direct_1, direct_2 = distance(trgt_haps[0], base_hap_1), distance(trgt_haps[1], base_hap_2)
            swapped_1, swapped_2 = distance(trgt_haps[0], base_hap_2), distance(trgt_haps[1], base_hap_1)
            direct, swapped = direct_1 + direct_2, swapped_1 + swapped_2
            direct_similarity = (
                similarity(direct_1, trgt_haps[0], base_hap_1)
                + similarity(direct_2, trgt_haps[1], base_hap_2)
            ) / 2
            swapped_similarity = (
                similarity(swapped_1, trgt_haps[0], base_hap_2)
                + similarity(swapped_2, trgt_haps[1], base_hap_1)
            ) / 2
            row.update({
                'direct_edit_distance': direct,
                'direct_haplotype_1_edit_distance': direct_1,
                'direct_haplotype_2_edit_distance': direct_2,
                'swapped_edit_distance': swapped,
                'swapped_haplotype_1_edit_distance': swapped_1,
                'swapped_haplotype_2_edit_distance': swapped_2,
                'distance_delta': abs(direct - swapped),
                'direct_similarity': f'{direct_similarity:.6f}',
                'swapped_similarity': f'{swapped_similarity:.6f}',
                'min_phase_edit_distance_delta': ~{min_phase_edit_distance_delta},
                'min_phase_similarity': f'{~{min_phase_similarity}:.6f}',
            })
            if direct == swapped:
                row.update({
                    'winning_edit_distance': direct,
                    'winning_similarity': f'{direct_similarity:.6f}',
                    'passes_edit_distance_delta': 'false',
                    'passes_similarity': str(direct_similarity >= ~{min_phase_similarity}).lower(),
                    'winning_orientation': 'tie',
                    'details': (
                        f'base_missing_as_reference_count={missing_as_reference};'
                        f'min_phase_edit_distance_delta=~{min_phase_edit_distance_delta};'
                        f'min_phase_similarity={~{min_phase_similarity}:.6f}'
                    ),
                })
                row['status'] = 'orientation_tie'
            else:
                direct_wins = direct < swapped
                winning_distance = direct if direct_wins else swapped
                winning_similarity = direct_similarity if direct_wins else swapped_similarity
                passes_delta = abs(direct - swapped) >= ~{min_phase_edit_distance_delta}
                passes_similarity = winning_similarity >= ~{min_phase_similarity}
                row.update({
                    'winning_edit_distance': winning_distance,
                    'winning_similarity': f'{winning_similarity:.6f}',
                    'passes_edit_distance_delta': str(passes_delta).lower(),
                    'passes_similarity': str(passes_similarity).lower(),
                    'winning_orientation': 'direct' if direct_wins else 'swapped',
                    'details': (
                        f'base_missing_as_reference_count={missing_as_reference};'
                        f'min_phase_edit_distance_delta=~{min_phase_edit_distance_delta};'
                        f'min_phase_similarity={~{min_phase_similarity}:.6f}'
                    ),
                })
                if not passes_delta and not passes_similarity:
                    row['status'] = 'phase_thresholds_not_met'
                elif not passes_delta:
                    row['status'] = 'edit_distance_delta_below_threshold'
                elif not passes_similarity:
                    row['status'] = 'winning_similarity_below_threshold'
                elif direct_wins:
                    call['GT'] = gt
                    call.phased = True
                    call['PS'] = rec.pos
                    row.update({'winning_phased_gt': gt_string(gt, True), 'phase_set': rec.pos, 'status': 'phased_direct'})
                    any_phased = True
                else:
                    call['GT'] = (gt[1], gt[0])
                    call.phased = True
                    call['PS'] = rec.pos
                    row.update({
                        'winning_phased_gt': gt_string((gt[1], gt[0]), True),
                        'phase_set': rec.pos,
                        'status': 'phased_swapped',
                    })
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
        disk_gb: 2 * ceil(size(vcf, "GB") + size(base_vcfs, "GB")) + 20,
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

# Replace selected TRVs, rebuild envelope relationships, and produce final VCFs and audit.
task ApplyTRLocusUpdates {
    input {
        File vcf
        File vcf_idx
        File? replacement_vcf
        File? replacement_vcf_idx
        File? replacement_map_tsv
        File? input_trv_phasing_summary_tsv
        File input_trv_catalog_match_tsv
        File? merged_status_tsv
        File gnomad_tr_json
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        # Stream the full contig while clearing and rebuilding TR and gnomAD STR annotations.
        python3 <<'PY'
import json
import os
import pysam

def record_key(rec):  # noqa: E302
    return rec.id if rec.id and rec.id != '.' else f'{rec.chrom}:{rec.pos}:{rec.ref}:{",".join(rec.alts or [])}'

def vals(value):  # noqa: E302
    return [str(item) for item in value] if isinstance(value, tuple) else ([] if value is None else [str(value)])

def record_end(rec):  # noqa: E302
    return rec.stop if rec.stop is not None else rec.pos + len(rec.ref) - 1

with open('~{gnomad_tr_json}') as handle:  # noqa: E305
    catalog = json.load(handle)
loci = []
for entry in catalog:
    # Prevent ignored non-disease loci from receiving gnomAD_STR during final annotation.
    diseases = entry.get('Diseases') if entry else None
    if not isinstance(diseases, list) or not diseases:
        continue
    explorer = entry.get('TRExplorerV1') if entry else None
    explorers = explorer if isinstance(explorer, list) else [explorer]
    for value in explorers:
        if entry and entry.get('LocusId') and value:
            loci.append((str(entry['LocusId']), str(value)))
locus_explorer = {locus: explorer for locus, explorer in loci}

# Preserve matches established before a replacement is merged/renamed. Every
# mapping here originates from the same strict substring test in prior tasks.
record_loci = {}
for audit_path in ['~{input_trv_catalog_match_tsv}', '~{merged_status_tsv}']:
    if not audit_path or not os.path.exists(audit_path):
        continue
    with open(audit_path) as audit:
        next(audit, None)
        for line in audit:
            fields = line.rstrip('\n').split('\t')
            if len(fields) < 4 or fields[3] == '.':
                continue
            for record_id in fields[3].split(','):
                record_loci.setdefault(record_id, set()).add(fields[1])

replacement_path = '~{replacement_vcf}'
map_path = '~{replacement_map_tsv}'
replacements = []
replace_old = set()
if replacement_path and os.path.exists(replacement_path):
    with pysam.VariantFile(replacement_path) as handle:
        replacements = [rec.copy() for rec in handle]
    for rec in replacements:
        searchable = '|'.join(vals(rec.info.get('TRID')))
        for locus, explorer in loci:
            if explorer in searchable:
                record_loci.setdefault(record_key(rec), set()).add(locus)
if map_path and os.path.exists(map_path):
    with open(map_path) as handle:
        next(handle, None)
        for line in handle:
            fields = line.rstrip('\n').split('\t')
            if len(fields) >= 4 and fields[3] == 'replace' and fields[1] != '.':
                replace_old.add(fields[1])

base = pysam.VariantFile('~{vcf}')
header = base.header.copy()
if replacements:
    # Preserve all annotations generated on the small replacement VCF.
    with pysam.VariantFile(replacement_path) as repl_header_source:
        header.merge(repl_header_source.header)
for name, number, type_, description in [
        ('TR_ENVELOPED', '0', 'Flag', 'Variant enveloped by tandem repeat'),
        ('TRID', '1', 'String', 'ID of enveloping tandem repeat'),
        ('gnomAD_STR', '1', 'String', 'Matched gnomAD tandem-repeat locus ID'),
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
if 'PS' not in header.formats:
    header.add_meta(
        'FORMAT',
        items=[('ID', 'PS'), ('Number', '1'), ('Type', 'Integer'), ('Description', 'Phase set')],
    )

contig_order = {name: i for i, name in enumerate(header.contigs)}
replacements.sort(key=lambda rec: (contig_order[rec.contig], rec.pos, record_end(rec)))
replacement_index = 0
with pysam.VariantFile('~{prefix}.pre_envelope.vcf', 'w', header=header) as out:
    for rec in base:
        while replacement_index < len(replacements) and (
            contig_order[replacements[replacement_index].contig] < contig_order[rec.contig]
            or (
                replacements[replacement_index].contig == rec.contig
                and replacements[replacement_index].pos <= rec.pos
            )
        ):
            replacement = replacements[replacement_index].copy()
            replacement.translate(header)
            out.write(replacement)
            replacement_index += 1
        if record_key(rec) not in replace_old:
            out.write(rec)
    while replacement_index < len(replacements):
        replacement = replacements[replacement_index].copy()
        replacement.translate(header)
        out.write(replacement)
        replacement_index += 1
base.close()

pre = pysam.VariantFile('~{prefix}.pre_envelope.vcf')
intervals = []
for rec in pre:
    if rec.info.get('allele_type') == 'trv':
        intervals.append((rec.contig, rec.pos, record_end(rec), record_key(rec)))
pre.close()
intervals.sort()
by_contig = {}
for interval in intervals:
    by_contig.setdefault(interval[0], []).append(interval)

src = pysam.VariantFile('~{prefix}.pre_envelope.vcf')
out = pysam.VariantFile('~{prefix}.trv_postprocessed.vcf.gz', 'wz', header=header)
with src, out:
    current_contig = None
    contig_intervals = []
    interval_index = 0
    active = []
    for rec in src:
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
        if 'gnomAD_STR' in rec.info:
            del rec.info['gnomAD_STR']
        if rec.info.get('allele_type') != 'trv':
            for _, start, stop, trid in active:
                if rec.pos >= start and record_end(rec) <= stop:
                    rec.info['TR_ENVELOPED'] = True
                    rec.info['TRID'] = trid
                    break
        searchable = '|'.join(vals(rec.info.get('TRID')))
        matches = sorted(
            set(record_loci.get(record_key(rec), set()))
            | {locus for locus, explorer in loci if rec.info.get('allele_type') == 'trv' and explorer in searchable}
        )
        if matches:
            rec.info['gnomAD_STR'] = ','.join(matches)
        out.write(rec)

with open('~{prefix}.trv_catalog_match.tsv', 'w') as out, open('~{input_trv_catalog_match_tsv}') as source:
    out.write(source.read())
status_path = '~{merged_status_tsv}'
if status_path and os.path.exists(status_path):
    with open(status_path) as source, open('~{prefix}.trv_catalog_match.tsv', 'a') as out:
        next(source, None)
        for line in source:
            out.write(line)
if map_path and os.path.exists(map_path):
    with open(map_path) as source, open('~{prefix}.trv_catalog_match.tsv', 'a') as out:
        next(source, None)
        for line in source:
            new_id, old_id, overlap_bp, status, phase_summary = line.rstrip('\n').split('\t')
            for locus in sorted(record_loci.get(new_id, [])):
                detail = f'{status};old_id={old_id};overlap_bp={overlap_bp};{phase_summary}'
                out.write(f'replacement\t{locus}\t{locus_explorer.get(locus, ".")}\t{new_id}\t{detail}\n')

# Always materialize phase audit so no-replacement contigs have stable workflow output.
phase_audit_path = '~{input_trv_phasing_summary_tsv}'
phase_header = (
    'record_id\tchrom\tpos\tend\tsample_id\tinput_gt'
    '\ttrgt_haplotype_1_sequence\ttrgt_haplotype_2_sequence\tbase_vcf_index\tbase_sample_id'
    '\tbase_variant_ids\tbase_missing_as_reference_count\tbase_haplotype_1_sequence'
    '\tbase_haplotype_2_sequence\tdirect_haplotype_1_edit_distance'
    '\tdirect_haplotype_2_edit_distance\tdirect_edit_distance\tswapped_haplotype_1_edit_distance'
    '\tswapped_haplotype_2_edit_distance\tswapped_edit_distance\tdistance_delta'
    '\tdirect_similarity\tswapped_similarity\twinning_orientation\twinning_phased_gt'
    '\twinning_edit_distance\twinning_similarity\tmin_phase_edit_distance_delta'
    '\tmin_phase_similarity\tpasses_edit_distance_delta\tpasses_similarity'
    '\tphase_set\tstatus\tdetails\n'
)
with open('~{prefix}.trv_phasing_summary.tsv', 'w') as out:
    if phase_audit_path and os.path.exists(phase_audit_path):
        with open(phase_audit_path) as source:
            out.write(source.read())
    else:
        out.write(phase_header)
PY
        rm -f ~{prefix}.pre_envelope.vcf
        tabix -f -p vcf ~{prefix}.trv_postprocessed.vcf.gz
        bcftools view -i 'INFO/allele_type="trv"' -Oz -o ~{prefix}.trv_updated.vcf.gz ~{prefix}.trv_postprocessed.vcf.gz
        tabix -f -p vcf ~{prefix}.trv_updated.vcf.gz
    >>>
    output {
        File trv_postprocessed_vcf = "~{prefix}.trv_postprocessed.vcf.gz"
        File trv_postprocessed_vcf_idx = "~{prefix}.trv_postprocessed.vcf.gz.tbi"
        File trv_updated_vcf = "~{prefix}.trv_updated.vcf.gz"
        File trv_updated_vcf_idx = "~{prefix}.trv_updated.vcf.gz.tbi"
        File trv_catalog_match_tsv = "~{prefix}.trv_catalog_match.tsv"
        File trv_phasing_summary_tsv = "~{prefix}.trv_phasing_summary.tsv"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 2,
        mem_gb: 12,
        disk_gb: 5 * ceil(size(vcf, "GB")) + 20,
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
