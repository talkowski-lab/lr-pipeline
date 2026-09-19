version 1.0

import "../annotation/AnnotateInSilicoPredictors.wdl"
import "../annotation/AnnotateGQMetrics.wdl"
import "../annotation/AnnotateRegion.wdl"
import "../annotation/AnnotateSQMetrics.wdl"
import "../annotation/AnnotateVRS.wdl"
import "../utils/Structs.wdl"
import "AnnotateVcf.wdl"
import "PostProcessTRLociHPRCHGSVC.wdl"

workflow PostProcessTRLociAoU {
    input {
        File vcf
        File vcf_idx
        String contig
        File trgt_vcf
        File trgt_vcf_idx
        String prefix

        File gnomad_tr_json
        File trgt_catalog_bed_gz

        Boolean run_flag_homopolymer_trvs
        Boolean replace_gnomad_str

        File seqrepo_tar
        File simple_repeats_bed
        File seg_dup_bed
        File repeat_masked_bed
        String cadd_ht
        String pangolin_ht
        String phylop_ht
        String revel_ht
        String spliceai_ht
        String annotate_in_silico_predictors_script = "https://raw.githubusercontent.com/talkowski-lab/lr-pipeline/main/scripts/annotation/annotate_insilico_predictors.py"
        String genome_build = "GRCh38"

        String utils_docker
        String vrs_docker
        String hail_docker

        RuntimeAttr? runtime_attr_prepare
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

    call PhaseReplacementLociAoU {
        input:
            vcf = vcf,
            vcf_idx = vcf_idx,
            trgt_vcf = trgt_vcf,
            trgt_vcf_idx = trgt_vcf_idx,
            trgt_catalog_bed_gz = trgt_catalog_bed_gz,
            gnomad_tr_json = gnomad_tr_json,
            contig = contig,
            run_flag_homopolymer_trvs = run_flag_homopolymer_trvs,
            prefix = "~{prefix}.~{contig}.replacement_seed",
            docker = utils_docker,
            runtime_attr_override = runtime_attr_prepare
    }

    if (PhaseReplacementLociAoU.retained_count > 0) {
        call AnnotateSQMetrics.CalculateSiteMetrics as CalculateReplacementSQMetrics {
            input:
                vcf = PhaseReplacementLociAoU.prepared_vcf,
                vcf_idx = PhaseReplacementLociAoU.prepared_vcf_idx,
                prefix = "~{prefix}.~{contig}.replacement.sq_metrics",
                docker = utils_docker,
                runtime_attr_override = runtime_attr_replacement_sq_metrics
        }

        call AnnotateGQMetrics.GenerateGQAnnotationTsv as CalculateReplacementSDMetrics {
            input:
                vcf = PhaseReplacementLociAoU.prepared_vcf,
                vcf_idx = PhaseReplacementLociAoU.prepared_vcf_idx,
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
                vcf = PhaseReplacementLociAoU.prepared_vcf,
                vcf_idx = PhaseReplacementLociAoU.prepared_vcf_idx,
                ab_bins = [0.00, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45, 0.50, 0.55, 0.60, 0.65, 0.70, 0.75, 0.80, 0.85, 0.90, 0.95, 1.00],
                prefix = "~{prefix}.~{contig}.replacement.ab_metrics",
                docker = utils_docker,
                runtime_attr_override = runtime_attr_replacement_ab_metrics
        }

        call AnnotateVRS.AnnotateVcfWithVRS as AnnotateReplacementVRS {
            input:
                vcf = PhaseReplacementLociAoU.prepared_vcf,
                vcf_idx = PhaseReplacementLociAoU.prepared_vcf_idx,
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
                vcf = PhaseReplacementLociAoU.prepared_vcf,
                vcf_idx = PhaseReplacementLociAoU.prepared_vcf_idx,
                simple_repeats_bed = simple_repeats_bed,
                seg_dup_bed = seg_dup_bed,
                repeat_masked_bed = repeat_masked_bed,
                prefix = "~{prefix}.~{contig}.replacement.region",
                docker = utils_docker,
                runtime_attr_override = runtime_attr_region_annotate
        }

        call AnnotateInSilicoPredictors.AnnotateInSilicoPredictorsTask as AnnotateReplacementInSilico {
            input:
                vcf = PhaseReplacementLociAoU.prepared_vcf,
                vcf_idx = PhaseReplacementLociAoU.prepared_vcf_idx,
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
                vcf = PhaseReplacementLociAoU.prepared_vcf,
                vcf_idx = PhaseReplacementLociAoU.prepared_vcf_idx,
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

    call PostProcessTRLociHPRCHGSVC.ApplyTRLocusUpdates {
        input:
            vcf = vcf,
            vcf_idx = vcf_idx,
            replacement_vcf = AttachReplacementAnnotations.annotated_vcf,
            replacement_vcf_idx = AttachReplacementAnnotations.annotated_vcf_idx,
            replacement_map_tsv = PhaseReplacementLociAoU.replacement_map_tsv,
            input_trv_catalog_match_tsv = PhaseReplacementLociAoU.trv_catalog_match_tsv,
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
    }
}

# Recover disease-associated catalog loci missing from the integrated VCF using a joint TRGT VCF and catalog BED
# Genotypes are taken unphased, exactly as TRGT reports them
task PhaseReplacementLociAoU {
    input {
        File vcf
        File vcf_idx
        File trgt_vcf
        File trgt_vcf_idx
        File trgt_catalog_bed_gz
        File gnomad_tr_json
        String contig
        Boolean run_flag_homopolymer_trvs
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        # Force the TRGT VCF onto the main VCF sample set, order, and FORMAT/AL so records translate cleanly downstream
        SAMPLES=$(bcftools query -l ~{vcf} | paste -sd, -)
        bcftools view -s "$SAMPLES" -Oz -o trgt.reordered.vcf.gz ~{trgt_vcf}
        tabix -f -p vcf trgt.reordered.vcf.gz

        python3 <<'PY'
import subprocess

main_header = subprocess.check_output(['bcftools', 'view', '-h', '~{vcf}'], text=True).splitlines()
incoming_header = subprocess.check_output(['bcftools', 'view', '-h', 'trgt.reordered.vcf.gz'], text=True).splitlines()
main_al = next((line for line in main_header if line.startswith('##FORMAT=<ID=AL,')), None)
if main_al:
    incoming_header = [main_al if line.startswith('##FORMAT=<ID=AL,') else line for line in incoming_header]
with open('trgt.header', 'w') as handle:
    handle.write('\n'.join(incoming_header) + '\n')
PY

        bcftools reheader -h trgt.header trgt.reordered.vcf.gz -o trgt.fixed.vcf.gz
        tabix -f -p vcf trgt.fixed.vcf.gz

        # Emit one report row per contig-relevant catalog entry recording how the disease locus was resolved
        python3 <<'PY'
import gzip
import json

import pysam

CONTIG = '~{contig}'
RUN_HOMOPOLYMER = ~{true="True" false="False" run_flag_homopolymer_trvs}


def norm_contig(value):  # noqa: E302
    return value[3:] if value.lower().startswith('chr') else value


def record_end(rec):  # noqa: E302
    return rec.stop if rec.stop is not None else rec.pos + len(rec.ref) - 1


def vals(value):  # noqa: E302
    return [str(item) for item in value] if isinstance(value, (list, tuple)) else ([] if value is None else [str(value)])


def trid_text(rec):  # noqa: E302
    return ','.join(vals(rec.info.get('TRID')))


def record_key(rec):  # noqa: E302
    return rec.id if rec.id and rec.id != '.' else f'{rec.chrom}:{rec.pos}:{rec.ref}:{",".join(rec.alts or [])}'


def overlap(left, right):  # noqa: E302
    return max(0, min(record_end(left), record_end(right)) - max(left.pos, right.pos) + 1)


def parse_explorer(value):  # noqa: E302
    fields = value.rsplit('-', 3)
    if len(fields) != 4:
        return None
    try:
        return fields[0], int(fields[1]), int(fields[2]), fields[3]
    except ValueError:
        return None


def shortest_motif_length(rec):  # noqa: E302
    """Match PostprocessCallset: a shortest MOTIFS element of length one is homopolymer."""
    motifs = rec.info.get('MOTIFS')
    if motifs is None:
        return None
    raw_values = motifs if isinstance(motifs, (list, tuple)) else [motifs]
    parts = [
        part
        for value in raw_values if value is not None
        for part in str(value).split(',') if part and part != '.'
    ]
    return min((len(part) for part in parts), default=None)


def recompute_ac(rec):  # noqa: E302
    counts = [0] * len(rec.alts or [])
    for sample_data in rec.samples.values():
        for allele in sample_data.get('GT') or ():
            if allele is not None and 0 < allele <= len(counts):
                counts[allele - 1] += 1
    rec.info['AC'] = tuple(counts)
    return sum(counts)


# Build a contig-local catalog BED of canonical coordinates keyed by embedded ID
bed_entries = []
with gzip.open('~{trgt_catalog_bed_gz}', 'rt') as handle:
    for line in handle:
        if not line.strip() or line.startswith('#'):
            continue
        cols = line.rstrip('\n').split('\t')
        if len(cols) < 4 or norm_contig(cols[0]) != norm_contig(CONTIG):
            continue
        bed_id = None
        for field in cols[3].split(';'):
            if field.startswith('ID='):
                bed_id = field[3:]
                break
        if bed_id is not None:
            bed_entries.append((bed_id, cols[0], int(cols[1]), int(cols[2])))


def find_bed(explorer):  # noqa: E302
    # Prefer an exact ID match, otherwise treat TRExplorerV1 as a substring of the catalog TRID
    for bed_id, _, start0, end in bed_entries:
        if bed_id == explorer:
            return start0, end
    for bed_id, _, start0, end in bed_entries:
        if explorer in bed_id:
            return start0, end
    return None


def contig_for(handle, contig):  # noqa: E302
    if contig in handle.header.contigs:
        return contig
    alternate = contig[3:] if contig.startswith('chr') else f'chr{contig}'
    return alternate if alternate in handle.header.contigs else None


main = pysam.VariantFile('~{vcf}', index_filename='~{vcf_idx}')
main_samples = list(main.header.samples)
trgt = pysam.VariantFile('trgt.fixed.vcf.gz')
if list(trgt.header.samples) != main_samples:
    raise RuntimeError('trgt_vcf samples must match main VCF sample set and order after subsetting')

for name, number, type_, description in [
    ('allele_type', '1', 'String', 'Allele type'),
    ('SOURCE', '1', 'String', 'Source of variant call'),
    # Declare allele_length because AnnotateRegion requires it before deriving TRV length from REF
    ('allele_length', '1', 'Integer', 'Allele length'),
    ('AC', 'A', 'Integer', 'Number of alleles observed'),
]:
    if name not in trgt.header.info:
        trgt.header.add_meta('INFO', items=[('ID', name), ('Number', number), ('Type', type_), ('Description', description)])
if 'PS' not in trgt.header.formats:
    trgt.header.add_meta('FORMAT', items=[('ID', 'PS'), ('Number', '1'), ('Type', 'Integer'), ('Description', 'Phase set')])
if RUN_HOMOPOLYMER and 'HOMOPOLYMER_TRV' not in trgt.header.info:
    trgt.header.info.add('HOMOPOLYMER_TRV', 0, 'Flag', 'Tandem repeat call where the shortest motif has length 1.')
out_header = trgt.header.copy()

main_contig = contig_for(main, CONTIG)
trgt_contig = contig_for(trgt, CONTIG)

with open('~{gnomad_tr_json}') as handle:
    catalog = json.load(handle)


def explorers_on_contig(entry):  # noqa: E302
    raw = entry.get('TRExplorerV1')
    raw = raw if isinstance(raw, list) else [raw]
    result = []
    for value in raw:
        parsed = parse_explorer(str(value)) if value else None
        if parsed and norm_contig(parsed[0]) == norm_contig(CONTIG):
            result.append(str(value))
    return result


tsv_columns = [
    'locus_id', 'TRExplorerV1', 'has_diseases', 'input_overlapping_trids',
    'input_has_TRExplorerV1_substring', 'trgt_overlapping_trids',
    'trgt_has_TRExplorerV1_substring', 'trgt_matching_allele_count', 'status',
]
tsv_rows = []
selected = []

for entry in catalog:
    if not entry or not entry.get('LocusId'):
        continue
    explorers = explorers_on_contig(entry)
    if not explorers:
        continue
    diseases = entry.get('Diseases')
    eligible = isinstance(diseases, list) and len(diseases) > 0
    row = dict.fromkeys(tsv_columns, '.')
    row['locus_id'] = str(entry['LocusId'])
    row['TRExplorerV1'] = ','.join(explorers)
    row['has_diseases'] = 'true' if eligible else 'false'
    row['status'] = 'not_eligible'
    if not eligible:
        tsv_rows.append(row)
        continue

    bed_match = None
    matched_explorer = None
    for explorer in explorers:
        found = find_bed(explorer)
        if found:
            bed_match = found
            matched_explorer = explorer
            break
    if bed_match is None:
        row['status'] = 'no_catalog_bed_match'
        tsv_rows.append(row)
        continue
    bed_start0, bed_end = bed_match
    # TRGT records carry a left anchor base, so the anchored POS is numerically the BED
    # 0-based start and POS+len(REF)-1 is the BED end. Verified against the catalog BED:
    # 33607/33607 input TRVs matched by TRID satisfy POS==col2 and POS+len(REF)-1==col3.
    locus_pos = bed_start0

    input_overlaps = []
    if main_contig is not None:
        for rec in main.fetch(main_contig, max(0, bed_start0 - 1), bed_end + 1):
            if rec.info.get('allele_type') == 'trv' and record_end(rec) >= locus_pos and rec.pos <= bed_end:
                input_overlaps.append(rec.copy())
    row['input_overlapping_trids'] = '|'.join(sorted({trid_text(r) for r in input_overlaps if trid_text(r)})) or '.'
    coord_matches = [r for r in input_overlaps if r.pos == locus_pos and record_end(r) == bed_end]
    row['input_has_TRExplorerV1_substring'] = 'true' if any(matched_explorer in trid_text(r) for r in input_overlaps) else 'false'
    if coord_matches:
        row['status'] = 'already_in_input_vcf'
        tsv_rows.append(row)
        continue

    trgt_hits = []
    if trgt_contig is not None:
        for rec in trgt.fetch(trgt_contig, max(0, bed_start0 - 1), bed_end + 1):
            if matched_explorer in trid_text(rec):
                trgt_hits.append(rec.copy())
    row['trgt_overlapping_trids'] = '|'.join(sorted({trid_text(r) for r in trgt_hits if trid_text(r)})) or '.'
    row['trgt_has_TRExplorerV1_substring'] = 'true' if trgt_hits else 'false'
    if not trgt_hits:
        row['status'] = 'no_trgt_match'
        tsv_rows.append(row)
        continue

    trgt_hits.sort(key=lambda r: (0 if (r.pos == locus_pos and record_end(r) == bed_end) else 1, r.pos))
    trec = trgt_hits[0]
    ac = recompute_ac(trec)
    row['trgt_matching_allele_count'] = str(ac)
    if ac < 1:
        row['status'] = 'trgt_ac0_skipped'
        tsv_rows.append(row)
        continue

    best = max(input_overlaps, key=lambda r: overlap(trec, r), default=None)
    old = best if (best is not None and overlap(trec, best) > 0) else None
    row['status'] = 'replaced_from_trgt' if old is not None else 'added_from_trgt'
    selected.append({'trec': trec, 'old': old})
    tsv_rows.append(row)

# Assign IDs exactly as IntegrateTRs.SetTrVariantIds does, including the _1/_2 suffixes on duplicates
id_counts = {}
for item in selected:
    item['base_id'] = f"{item['trec'].chrom}-{item['trec'].pos}-TRV-{len(item['trec'].ref) - 1}"
    id_counts[item['base_id']] = id_counts.get(item['base_id'], 0) + 1

id_seen = {}
used_old = set()
records = []
map_rows = []
for item in selected:
    rec = item['trec']
    base_id = item['base_id']
    if id_counts[base_id] > 1:
        id_seen[base_id] = id_seen.get(base_id, 0) + 1
        rec.id = f'{base_id}_{id_seen[base_id]}'
    else:
        rec.id = base_id
    rec.info['allele_type'] = 'trv'
    rec.info['SOURCE'] = 'TRExplorer'
    if RUN_HOMOPOLYMER and shortest_motif_length(rec) == 1:
        rec.info['HOMOPOLYMER_TRV'] = True
    # Normalize only recovered records to PASS because TRGT leaves FILTER unset, preserving named filters
    if not tuple(rec.filter.keys()):
        rec.filter.add('PASS')
    new_key = record_key(rec)
    if item['old'] is not None:
        old_key = record_key(item['old'])
        if old_key in used_old:
            raise RuntimeError(f'Multiple recovered loci selected main VCF TRV {old_key}')
        used_old.add(old_key)
        map_rows.append((new_key, old_key, str(overlap(rec, item['old'])), 'replace', '.',
                         item['old'].contig, str(item['old'].pos), str(record_end(item['old']))))
    else:
        map_rows.append((new_key, '.', '0', 'add', '.', '.', '.', '.'))
    records.append(rec)

contig_order = {name: i for i, name in enumerate(out_header.contigs)}
records.sort(key=lambda rec: (contig_order[rec.contig], rec.pos, record_end(rec)))
with pysam.VariantFile('~{prefix}.vcf.gz', 'wz', header=out_header) as out:
    for rec in records:
        out.write(rec)

with open('~{prefix}.map.tsv', 'w') as handle:
    handle.write('new_record_id\told_record_id\toverlap_bp\tstatus\tphase_summary\told_contig\told_pos\told_end\n')
    for map_row in map_rows:
        handle.write('\t'.join(map_row) + '\n')

with open('~{prefix}.trv_catalog_match.tsv', 'w') as handle:
    handle.write('\t'.join(tsv_columns) + '\n')
    for row in tsv_rows:
        handle.write('\t'.join(str(row[column]) for column in tsv_columns) + '\n')

main.close()
trgt.close()
PY

        rm -f trgt.reordered.vcf.gz trgt.reordered.vcf.gz.tbi trgt.header trgt.fixed.vcf.gz trgt.fixed.vcf.gz.tbi
        tabix -f -p vcf ~{prefix}.vcf.gz
        bcftools view -H ~{prefix}.vcf.gz | wc -l
    >>>

    output {
        File prepared_vcf = "~{prefix}.vcf.gz"
        File prepared_vcf_idx = "~{prefix}.vcf.gz.tbi"
        File replacement_map_tsv = "~{prefix}.map.tsv"
        File trv_catalog_match_tsv = "~{prefix}.trv_catalog_match.tsv"
        Int retained_count = read_int(stdout())
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 2,
        mem_gb: 8,
        disk_gb: 2 * ceil(size(vcf, "GB") + size(trgt_vcf, "GB") + size(trgt_catalog_bed_gz, "GB")) + 20,
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
