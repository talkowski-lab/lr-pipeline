version 1.0

import "../utils/Helpers.wdl"
import "../utils/Structs.wdl"

workflow IntegrateTRs {
    meta {
        description: [
            "This utility integrates tandem-repeat calls into a base VCF for a cohort. It aligns samples between the base and TR VCFs, sets missing filters to pass, tags TR records with their source catalog, assigns TR identifiers and annotates the base VCF with the integrated TR calls. It outputs the TR-annotated VCF."
        ]
    }

    parameter_meta {
        vcf: "Base VCF to integrate TRs into."
        vcf_idx: "Index for the base VCF."
        tr_vcf: "Tandem-repeat VCF to integrate."
        tr_vcf_idx: "Index for the TR VCF."
        contigs: "Contigs to process."
        sample_ids: "Samples shared between the base and TR VCFs."
        tr_catalogs: "Catalogs from which the TR calls were derived."
        tr_catalog_ids: "Identifier for each catalog in `tr_catalogs`."
        tr_annotated_vcf: "Base VCF annotated with integrated TR calls."
        tr_annotated_vcf_idx: "Index for the annotated VCF."
    }

    input {
        File vcf
        File vcf_idx
        File tr_vcf
        File tr_vcf_idx
        Array[String] contigs
        String prefix

        Array[String] sample_ids
        Array[File] tr_catalogs
        Array[String] tr_catalog_ids

        String utils_docker

        RuntimeAttr? runtime_attr_subset_contig_base
        RuntimeAttr? runtime_attr_subset_contig_tr
        RuntimeAttr? runtime_attr_subset_samples_base
        RuntimeAttr? runtime_attr_subset_samples_tr
        RuntimeAttr? runtime_attr_check_sample_consistency
        RuntimeAttr? runtime_attr_set_missing_filters
        RuntimeAttr? runtime_attr_subset_catalog
        RuntimeAttr? runtime_attr_tag_tr_vcf
        RuntimeAttr? runtime_attr_set_tr_ids
        RuntimeAttr? runtime_attr_annotate_vcf
        RuntimeAttr? runtime_attr_concat_vcf
    }

    Boolean single_contig = length(contigs) == 1

    scatter (contig in contigs) {
        if (!single_contig) {
            call Helpers.SubsetVcfToContig as SubsetContigBase {
                input:
                    vcf = vcf,
                    vcf_idx = vcf_idx,
                    contig = contig,
                    prefix = "~{prefix}.~{contig}.base",
                    docker = utils_docker,
                    runtime_attr_override = runtime_attr_subset_contig_base
            }

            call Helpers.SubsetVcfToContig as SubsetContigTr {
                input:
                    vcf = tr_vcf,
                    vcf_idx = tr_vcf_idx,
                    contig = contig,
                    prefix = "~{prefix}.~{contig}.tr",
                    docker = utils_docker,
                    runtime_attr_override = runtime_attr_subset_contig_tr
            }
        }

        File contig_vcf = select_first([SubsetContigBase.subset_vcf, vcf])
        File contig_vcf_idx = select_first([SubsetContigBase.subset_vcf_idx, vcf_idx])

        File contig_tr_vcf = select_first([SubsetContigTr.subset_vcf, tr_vcf])
        File contig_tr_vcf_idx = select_first([SubsetContigTr.subset_vcf_idx, tr_vcf_idx])

        call Helpers.SubsetVcfToSamples as SubsetSamplesBase {
            input:
                vcf = contig_vcf,
                vcf_idx = contig_vcf_idx,
                samples = sample_ids,
                prefix = "~{prefix}.~{contig}.base_subset",
                docker = utils_docker,
                runtime_attr_override = runtime_attr_subset_samples_base
        }

        call Helpers.SubsetVcfToSamples as SubsetSamplesTr {
            input:
                vcf = contig_tr_vcf,
                vcf_idx = contig_tr_vcf_idx,
                samples = sample_ids,
                prefix = "~{prefix}.~{contig}.tr_subset",
                docker = utils_docker,
                runtime_attr_override = runtime_attr_subset_samples_tr
        }

        call Helpers.CheckSampleConsistency {
            input:
                vcfs = [SubsetSamplesBase.subset_vcf, SubsetSamplesTr.subset_vcf],
                vcf_idxs = [SubsetSamplesBase.subset_vcf_idx, SubsetSamplesTr.subset_vcf_idx],
                sample_ids = sample_ids,
                docker = utils_docker,
                runtime_attr_override = runtime_attr_check_sample_consistency
        }

        call Helpers.SetMissingFiltersToPass {
            input:
                vcf = SubsetSamplesTr.subset_vcf,
                vcf_idx = SubsetSamplesTr.subset_vcf_idx,
                prefix = "~{prefix}.~{contig}.tr_subset.refiltered",
                docker = utils_docker,
                runtime_attr_override = runtime_attr_set_missing_filters
        }

        scatter (i in range(length(tr_catalogs))) {
            call Helpers.SubsetTsvToContig {
                input:
                    tsv = tr_catalogs[i],
                    contig = contig,
                    compressed_tsv = true,
                    prefix = "~{prefix}.~{contig}.catalog~{i}",
                    docker = utils_docker,
                    runtime_attr_override = runtime_attr_subset_catalog
            }
        }

        call TagTrVcfWithCatalogs {
            input:
                vcf = SetMissingFiltersToPass.filtered_vcf,
                vcf_idx = SetMissingFiltersToPass.filtered_vcf_idx,
                catalogs = SubsetTsvToContig.subset_tsv,
                catalog_ids = tr_catalog_ids,
                prefix = "~{prefix}.~{contig}.tr_subset.tagged",
                docker = utils_docker,
                runtime_attr_override = runtime_attr_tag_tr_vcf
        }

        call SetTrVariantIds {
            input:
                vcf = TagTrVcfWithCatalogs.tagged_vcf,
                vcf_idx = TagTrVcfWithCatalogs.tagged_vcf_idx,
                prefix = "~{prefix}.~{contig}.tr_subset.ids",
                docker = utils_docker,
                runtime_attr_override = runtime_attr_set_tr_ids
        }

        call AnnotateVcfWithTRs {
            input:
                vcf = SubsetSamplesBase.subset_vcf,
                vcf_idx = SubsetSamplesBase.subset_vcf_idx,
                tr_vcf = SetTrVariantIds.renamed_vcf,
                tr_vcf_idx = SetTrVariantIds.renamed_vcf_idx,
                prefix = "~{prefix}.~{contig}.merged",
                docker = utils_docker,
                runtime_attr_override = runtime_attr_annotate_vcf
        }
    }

    if (!single_contig) {
        call Helpers.ConcatVcfs {
            input:
                vcfs = AnnotateVcfWithTRs.annotated_vcf,
                vcf_idxs = AnnotateVcfWithTRs.annotated_vcf_idx,
                allow_overlaps = false,
                naive = true,
                prefix = "~{prefix}.annotated_trs",
                docker = utils_docker,
                runtime_attr_override = runtime_attr_concat_vcf
        }
    }

    output {
        File tr_annotated_vcf = select_first([ConcatVcfs.concat_vcf, AnnotateVcfWithTRs.annotated_vcf[0]])
        File tr_annotated_vcf_idx = select_first([ConcatVcfs.concat_vcf_idx, AnnotateVcfWithTRs.annotated_vcf_idx[0]])
    }
}

task TagTrVcfWithCatalogs {
    input {
        File vcf
        File vcf_idx
        Array[File] catalogs
        Array[String] catalog_ids
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        python3 <<CODE
from pysam import VariantFile

catalog_files = "~{sep=',' catalogs}".split(',')
catalog_ids_list = "~{sep=',' catalog_ids}".split(',')

lookup = {}
for cat_file, cat_id in zip(catalog_files, catalog_ids_list):
    with open(cat_file, 'r') as f:
        for line in f:
            parts = line.strip().split('\t')
            trid = None
            for field in parts[3].split(';'):
                if field.startswith('ID='):
                    trid = field[3:]
                    break
            if trid:
                lookup[trid] = cat_id

vcf_in = VariantFile("~{vcf}")
if 'allele_type' not in vcf_in.header.info:
    vcf_in.header.add_meta('INFO', items=[('ID', 'allele_type'), ('Number', 1), ('Type', 'String'), ('Description', 'Allele type')])
if 'SOURCE' not in vcf_in.header.info:
    vcf_in.header.add_meta('INFO', items=[('ID', 'SOURCE'), ('Number', 1), ('Type', 'String'), ('Description', 'Source of variant call')])

vcf_out = VariantFile("~{prefix}.vcf.gz", "w", header=vcf_in.header)
for record in vcf_in:
    record.info['allele_type'] = 'trv'
    trid = record.info.get('TRID')
    if trid is not None:
        trid_key = ','.join(trid) if isinstance(trid, tuple) else trid
        source = lookup.get(trid_key)
        if source:
            record.info['SOURCE'] = source
    vcf_out.write(record)

vcf_in.close()
vcf_out.close()
CODE

        tabix -p vcf ~{prefix}.vcf.gz
    >>>

    output {
        File tagged_vcf = "~{prefix}.vcf.gz"
        File tagged_vcf_idx = "~{prefix}.vcf.gz.tbi"
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

task SetTrVariantIds {
    input {
        File vcf
        File vcf_idx
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        python3 <<CODE
from pysam import VariantFile
from collections import defaultdict

vcf_in = VariantFile("~{vcf}")
id_counts = defaultdict(int)
for record in vcf_in:
    new_id = f"{record.chrom}-{record.pos}-TRV-{len(record.ref)-1}"
    id_counts[new_id] += 1
vcf_in.close()

vcf_in = VariantFile("~{vcf}")
vcf_out = VariantFile("~{prefix}.vcf.gz", "w", header=vcf_in.header)
id_seen = defaultdict(int)
for record in vcf_in:
    new_id = f"{record.chrom}-{record.pos}-TRV-{len(record.ref)-1}"
    if id_counts[new_id] > 1:
        id_seen[new_id] += 1
        record.id = f"{new_id}_{id_seen[new_id]}"
    else:
        record.id = new_id
    vcf_out.write(record)
vcf_in.close()
vcf_out.close()
CODE

        tabix -p vcf ~{prefix}.vcf.gz
    >>>

    output {
        File renamed_vcf = "~{prefix}.vcf.gz"
        File renamed_vcf_idx = "~{prefix}.vcf.gz.tbi"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 4,
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

task AnnotateVcfWithTRs {
    input {
        File vcf
        File vcf_idx
        File tr_vcf
        File tr_vcf_idx
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        python3 <<CODE
from pysam import VariantFile

# Load TR intervals from the TR VCF
tr_in = VariantFile("~{tr_vcf}")
tr_intervals = []
for record in tr_in:
    tr_start = record.pos
    tr_end = record.pos + len(record.ref)
    tr_id = record.id if record.id else f"{record.chrom}-{record.pos}-TRV-{len(record.ref)-1}"
    tr_intervals.append((tr_start, tr_end, tr_id))
tr_in.close()

# Annotate the base VCF using a sliding window over TR intervals
vcf_in = VariantFile("~{vcf}")
if 'TR_ENVELOPED' not in vcf_in.header.info:
    vcf_in.header.add_meta('INFO', items=[('ID', 'TR_ENVELOPED'), ('Number', 0), ('Type', 'Flag'), ('Description', 'Variant enveloped by tandem repeat')])
if 'TRID' not in vcf_in.header.info:
    vcf_in.header.add_meta('INFO', items=[('ID', 'TRID'), ('Number', 1), ('Type', 'String'), ('Description', 'ID of enveloping tandem repeat')])

vcf_out = VariantFile("vcf_annotated.vcf.gz", "w", header=vcf_in.header)
annotated_count = 0
tr_idx = 0
active_trs = []

for record in vcf_in:
    v_start = record.pos
    v_end = record.pos + len(record.ref)

    # Add new TR intervals that start at or before this variant
    while tr_idx < len(tr_intervals) and tr_intervals[tr_idx][0] <= v_start:
        active_trs.append(tr_intervals[tr_idx])
        tr_idx += 1

    # Remove expired TR intervals that end before this variant's start
    active_trs = [t for t in active_trs if t[1] >= v_start]

    # Check containment: variant must be fully within a TR interval
    for tr_start, tr_end, tr_id in active_trs:
        if v_start >= tr_start and v_end <= tr_end:
            record.info['TR_ENVELOPED'] = True
            record.info['TRID'] = tr_id
            annotated_count += 1
            break

    vcf_out.write(record)

vcf_in.close()
vcf_out.close()
CODE

        tabix -p vcf vcf_annotated.vcf.gz

        bcftools concat \
            --allow-overlaps \
            -Oz -o ~{prefix}.vcf.gz \
            vcf_annotated.vcf.gz \
            ~{tr_vcf}

        tabix -p vcf ~{prefix}.vcf.gz

        rm -f vcf_annotated.vcf.gz vcf_annotated.vcf.gz.tbi
    >>>

    output {
        File annotated_vcf = "~{prefix}.vcf.gz"
        File annotated_vcf_idx = "~{prefix}.vcf.gz.tbi"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 2,
        mem_gb: 8,
        disk_gb: 15 * ceil(size(vcf, "GB") + size(tr_vcf, "GB")) + 20,
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
