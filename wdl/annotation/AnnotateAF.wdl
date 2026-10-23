# Derived from broadinstitute/gatk-sv (branch kj_project_gnomad_lr) wdl/AnnotateAF.wdl and wdl/ShardedAnnotateVcf.wdl
version 1.0

import "../utils/Helpers.wdl"
import "../utils/Structs.wdl"

workflow AnnotateAF {
    meta {
        description: [
            "This workflow annotates each variant with allele counts, allele frequencies, no-call rates and genotype counts, computed from its genotypes. Frequencies are computed across the whole cohort and, when the matching inputs are given, per sex from a PED and per population from an ancestry table, together with the population with the highest allele frequency (grpmax) at biallelic sites. Variants on chrX and chrY outside the pseudoautosomal regions are counted as hemizygous in XY samples.",
            "At multiallelic sites, per-allele allele purity (`AP_allele`), motif count (`MC_allele`) and, when an LPS table is given, longest polymer sequence (`LPS_allele`) are taken from the first carrier of each allele. Each contig is sharded by record count, annotated in parallel and concatenated into one VCF."
        ]
    }

    parameter_meta {
        vcf: "VCF to annotate."
        vcf_idx: "Index for VCF to annotate."
        contigs: "Contigs to annotate within the input VCF."
        records_per_shard: "Number of variants to keep within a single shard during annotation."
        sample_pop_assignments: "Two-column TSV of sample ID and population label, used for per-population frequencies. Labels must be one of afr, ami, amr, asj, eas, fin, mid, nfe, sas or oth, and '.' leaves a sample unassigned."
        ped: "Cohort pedigree, used for per-sex frequencies and hemizygous counting."
        lps_tsv: "TSV of longest polymer sequence values per TR locus and sample, from `TRGTLPS`, with its header line first."
        strip_info_fields: "INFO fields to remove from each shard before annotation."
        par_bed: "From references."
        af_annotated_vcf: "VCF annotated with allele frequencies."
        af_annotated_vcf_idx: "Index for `af_annotated_vcf`."
    }

    input {
        File vcf
        File vcf_idx
        Array[String] contigs
        String prefix

        Int records_per_shard
        File? sample_pop_assignments
        File? ped
        File? lps_tsv
        Array[String]? strip_info_fields

        File par_bed

        String utils_docker

        RuntimeAttr? runtime_attr_subset_lps_tsv
        RuntimeAttr? runtime_attr_subset
        RuntimeAttr? runtime_attr_shard
        RuntimeAttr? runtime_attr_strip_info_fields
        RuntimeAttr? runtime_attr_compute_afs
        RuntimeAttr? runtime_attr_concat
    }

    Boolean single_contig = length(contigs) == 1

    scatter (c in range(length(contigs))) {
        String contig = contigs[c]

        if (defined(lps_tsv)) {
            call SubsetLpsTsvToContig {
                input:
                    tsv = select_first([lps_tsv]),
                    contig = contig,
                    prefix = "~{prefix}.~{contig}",
                    docker = utils_docker,
                    runtime_attr_override = runtime_attr_subset_lps_tsv
            }
        }

        if (!single_contig) {
            call Helpers.SubsetVcfToContig as SubsetVcf {
                input:
                    vcf = vcf,
                    vcf_idx = vcf_idx,
                    contig = contig,
                    prefix = "~{prefix}.~{contig}",
                    docker = utils_docker,
                    runtime_attr_override = runtime_attr_subset
            }
        }

        call Helpers.ShardVcfByRecords {
            input:
                vcf = select_first([SubsetVcf.subset_vcf, vcf]),
                vcf_idx = select_first([SubsetVcf.subset_vcf_idx, vcf_idx]),
                records_per_shard = records_per_shard,
                prefix = "~{prefix}.~{contig}",
                docker = utils_docker,
                runtime_attr_override = runtime_attr_shard
        }

        scatter (shard in ShardVcfByRecords.shards) {
            Int shard_contig_idx = c
        }
    }

    # Annotate every shard in one flat scatter, since Cromwell can't read unset optional inputs inside a nested scatter
    Array[File] shards = flatten(ShardVcfByRecords.shards)
    Array[File] shard_idxs = flatten(ShardVcfByRecords.shard_idxs)
    Array[Int] shard_contig_idxs = flatten(shard_contig_idx)
    Array[File?] contig_lps_tsvs = SubsetLpsTsvToContig.subset_tsv

    scatter (i in range(length(shards))) {
        Int contig_idx = shard_contig_idxs[i]

        if (defined(strip_info_fields)) {
            call Helpers.StripInfoFields {
                input:
                    vcf = shards[i],
                    vcf_idx = shard_idxs[i],
                    info_fields = select_first([strip_info_fields]),
                    prefix = "~{prefix}.~{i}.stripped",
                    docker = utils_docker,
                    runtime_attr_override = runtime_attr_strip_info_fields
            }
        }

        call ComputeAFs {
            input:
                vcf = select_first([StripInfoFields.stripped_vcf, shards[i]]),
                sample_pop_assignments = sample_pop_assignments,
                ped = ped,
                lps_tsv = contig_lps_tsvs[contig_idx],
                par_bed = par_bed,
                prefix = "~{prefix}.~{i}.AF",
                docker = utils_docker,
                runtime_attr_override = runtime_attr_compute_afs
        }
    }

    call Helpers.ConcatVcfs {
        input:
            vcfs = ComputeAFs.af_vcf,
            vcf_idxs = ComputeAFs.af_vcf_idx,
            allow_overlaps = true,
            naive = false,
            no_version = true,
            prefix = "~{prefix}.annotated",
            docker = utils_docker,
            runtime_attr_override = runtime_attr_concat
    }

    output {
        File af_annotated_vcf = ConcatVcfs.concat_vcf
        File af_annotated_vcf_idx = ConcatVcfs.concat_vcf_idx
    }
}

task SubsetLpsTsvToContig {
    input {
        File tsv
        String contig
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        # Keep the header line and the loci whose TRID contig matches, ignoring any chr prefix
        zcat -f ~{tsv} \
            | awk -F'\t' -v contig="~{contig}" '
                BEGIN {
                    if (contig ~ /^chr/) contig = substr(contig, 4)
                }

                {
                    split($1, arr, ",")
                    split(arr[1], parts, "-")
                    chr = parts[1]
                    if (chr ~ /^chr/) chr = substr(chr, 4)
                    if (chr ~ /^trid/) header = $0
                    else if (chr == contig) data[++n] = $0
                }

                END {
                    print header
                    for (i = 1; i <= n; i++) print data[i]
                }
            ' \
            > ~{prefix}.tsv
    >>>

    output {
        File subset_tsv = "~{prefix}.tsv"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 4,
        disk_gb: 2 * ceil(size(tsv, "GB")) + 5,
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

task ComputeAFs {
    input {
        File vcf
        File? sample_pop_assignments
        File? ped
        File? lps_tsv
        File par_bed
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        python3 /opt/scripts/annotation/compute_AFs.py "~{vcf}" stdout \
            ~{"-p " + sample_pop_assignments} \
            ~{"-f " + ped} \
            ~{"-l " + lps_tsv} \
            --par ~{par_bed} \
            | bgzip -c > "~{prefix}.vcf.gz"

        tabix -p vcf "~{prefix}.vcf.gz"
    >>>

    output {
        File af_vcf = "~{prefix}.vcf.gz"
        File af_vcf_idx = "~{prefix}.vcf.gz.tbi"
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
