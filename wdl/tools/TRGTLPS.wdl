version 1.0

import "../utils/Helpers.wdl"
import "../utils/Structs.wdl"

workflow TRGTLPS {
    meta {
        description: [
            "This tool runs the trgt-lps (https://github.com/PacificBiosciences/trgt-lps) tool per contig to compute the longest polymer sequence (LPS) within each TRGT-genotyped tandem-repeat locus for every sample, concatenating the results into a single TSV. Alongside each contig's LPS table it extracts a small TRID-metadata TSV from the same subset VCF, mapping every `(TRID, motif)` to the LocusIds that record covers, which `CreateTRGTHistograms` needs to resolve variation-cluster records whose TRID names several loci. It outputs the LPS TSV and the per-contig TRID-metadata TSVs."
        ]
    }

    parameter_meta {
        vcf: "TRGT VCF to process."
        vcf_idx: "Index for the TRGT VCF."
        contigs: "Contigs to process."
        normalize_chry_haploid_genotypes: "Whether to normalize haploid chrY genotypes before computing the longest polymer sequence. Applied to the chrY shard only."
        filter_trid_motif_pairs: "`(TRID, motif)` pairs to drop from the concatenated LPS table, given as two-element arrays of the LPS table's `trid` and `motif` column values - e.g. `[['X-149631602-149631617-TCC,X-149631685-149631694-GCT,X-149631723-149631735-CGCCGT', 'CGC']]`. Use it for a row trgt-lps emitted from a spurious `INFO/MOTIFS` value, which `CreateTRGTHistograms` cannot resolve because no LocusId in the TRID carries that motif. Every pair must match at least one LPS row or the task fails. Pass an empty array to filter nothing."
        trgt_lps_tsv: "TSV of per-locus longest polymer sequences."
        vcf_trid_metadata_tsvs: "Per-contig TRID-metadata TSVs, index-aligned with the `contigs` input array, to be passed straight to `CreateTRGTHistograms`."
    }

    input {
        File vcf
        File vcf_idx
        Array[String] contigs
        String prefix

        Boolean normalize_chry_haploid_genotypes
        Array[Array[String]] filter_trid_motif_pairs

        String trgt_lps_docker
        String utils_docker
        String stranalysis_docker

        RuntimeAttr? runtime_attr_subset_vcf
        RuntimeAttr? runtime_attr_normalize_chry_haploid_genotypes
        RuntimeAttr? runtime_attr_add_end
        RuntimeAttr? runtime_attr_trgt_lps
        RuntimeAttr? runtime_attr_extract_trid_metadata
        RuntimeAttr? runtime_attr_concat
        RuntimeAttr? runtime_attr_filter_lps_rows
    }

    scatter (contig in contigs) {
        Boolean normalize_chry_haploid_genotypes_for_contig = normalize_chry_haploid_genotypes && contig == "chrY"

        call Helpers.SubsetVcfToContig {
            input:
                vcf = vcf,
                vcf_idx = vcf_idx,
                contig = contig,
                prefix = "~{prefix}.~{contig}",
                docker = utils_docker,
                runtime_attr_override = runtime_attr_subset_vcf
        }

        if (normalize_chry_haploid_genotypes_for_contig) {
            call Helpers.NormalizeTRGTHaploidGenotypes {
                input:
                    vcf = SubsetVcfToContig.subset_vcf,
                    vcf_idx = SubsetVcfToContig.subset_vcf_idx,
                    prefix = "~{prefix}.~{contig}.haploid_genotypes",
                    docker = utils_docker,
                    runtime_attr_override = runtime_attr_normalize_chry_haploid_genotypes
            }
        }

        File trgt_lps_vcf = select_first([NormalizeTRGTHaploidGenotypes.normalized_vcf, SubsetVcfToContig.subset_vcf])
        File trgt_lps_vcf_idx = select_first([NormalizeTRGTHaploidGenotypes.normalized_vcf_idx, SubsetVcfToContig.subset_vcf_idx])

        call RunTRGTLPS {
            input:
                vcf = trgt_lps_vcf,
                vcf_idx = trgt_lps_vcf_idx,
                prefix = "~{prefix}.~{contig}.lps",
                docker = trgt_lps_docker,
                runtime_attr_override = runtime_attr_trgt_lps
        }

        call Helpers.AddTREndTag {
            input:
                vcf = SubsetVcfToContig.subset_vcf,
                vcf_idx = SubsetVcfToContig.subset_vcf_idx,
                prefix = "~{prefix}.~{contig}.with_end",
                docker = utils_docker,
                runtime_attr_override = runtime_attr_add_end
        }

        call ExtractTridMetadata {
            input:
                vcf = AddTREndTag.vcf_with_end,
                vcf_idx = AddTREndTag.vcf_with_end_idx,
                contig = contig,
                prefix = "~{prefix}.~{contig}.trid_metadata",
                docker = stranalysis_docker,
                runtime_attr_override = runtime_attr_extract_trid_metadata
        }
    }

    call Helpers.ConcatTsvs {
        input:
            tsvs = RunTRGTLPS.lps_tsv,
            sort_output = false,
            preserve_header = true,
            prefix = "~{prefix}.concat",
            docker = utils_docker,
            runtime_attr_override = runtime_attr_concat
    }

    if (length(filter_trid_motif_pairs) > 0) {
        call FilterLpsTsvRows {
            input:
                tsv = ConcatTsvs.concatenated_tsv,
                filter_trid_motif_pairs = filter_trid_motif_pairs,
                prefix = "~{prefix}.concat.filtered",
                docker = utils_docker,
                runtime_attr_override = runtime_attr_filter_lps_rows
        }
    }

    output {
        File trgt_lps_tsv = select_first([FilterLpsTsvRows.filtered_tsv, ConcatTsvs.concatenated_tsv])
        Array[File] vcf_trid_metadata_tsvs = ExtractTridMetadata.trid_metadata_tsv
    }
}

task RunTRGTLPS {
    input {
        File vcf
        File vcf_idx
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        trgt-lps \
            --vcf ~{vcf} \
            --threads ~{select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])} \
            > ~{prefix}.tsv
    >>>

    output {
        File lps_tsv = "~{prefix}.tsv"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 8,
        mem_gb: 4,
        disk_gb: 2 * ceil(size(vcf, "GB")) + 20,
        boot_disk_gb: 10,
        preemptible_tries: 1,
        max_retries: 0
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
        memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GB"
        disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " SSD"
        bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
        docker: docker
        preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
    }
}

task ExtractTridMetadata {
    input {
        File vcf
        File vcf_idx
        String contig
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        ln -s ~{vcf} input.vcf.gz
        ln -s ~{vcf_idx} input.vcf.gz.tbi

        python3 -m str_analysis.extract_trid_metadata_from_TRGT_vcf \
            --input-vcf input.vcf.gz \
            --output-tsv ~{prefix}.tsv.gz \
            --contig ~{contig}
    >>>

    output {
        File trid_metadata_tsv = "~{prefix}.tsv.gz"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 2,
        mem_gb: 8,
        disk_gb: ceil(size(vcf, "GB")) + 20,
        boot_disk_gb: 10,
        preemptible_tries: 1,
        max_retries: 0
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
        memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GB"
        disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
        bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
        docker: docker
        preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
    }
}

# Drop hand-listed (trid, motif) rows, e.g. one trgt-lps emitted from a spurious INFO/MOTIFS value
task FilterLpsTsvRows {
    input {
        File tsv
        Array[Array[String]] filter_trid_motif_pairs
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        awk -F'\t' '
            # Every pair is read before the LPS table, so a pair matching no row can be reported at the end
            NR == FNR {
                if (NF != 2) {
                    printf "Filter pair on line %d has %d field(s), expected a trid and a motif\n", FNR, NF > "/dev/stderr"
                    malformed = 1
                    exit 1
                }
                matched[$1 SUBSEP $2] = 0
                next
            }

            FNR == 1 {
                print
                next
            }

            ($1 SUBSEP $2) in matched {
                matched[$1 SUBSEP $2] += 1
                next
            }

            {
                print
            }

            END {
                if (malformed) {
                    exit 1
                }

                for (key in matched) {
                    split(key, pair, SUBSEP)
                    printf "Filtered %d LPS row(s) for pair (%s, %s)\n", matched[key], pair[1], pair[2] > "/dev/stderr"
                    if (matched[key] == 0) {
                        unmatched = 1
                    }
                }

                if (unmatched) {
                    printf "Every filter pair must match at least one LPS row; the LPS table and the filter list disagree\n" > "/dev/stderr"
                    exit 1
                }
            }
        ' ~{write_tsv(filter_trid_motif_pairs)} ~{tsv} \
            > ~{prefix}.tsv
    >>>

    output {
        File filtered_tsv = "~{prefix}.tsv"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 4,
        disk_gb: 2 * ceil(size(tsv, "GB")) + 10,
        boot_disk_gb: 10,
        preemptible_tries: 1,
        max_retries: 0
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
        memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GB"
        disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
        bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
        docker: docker
        preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
    }
}
