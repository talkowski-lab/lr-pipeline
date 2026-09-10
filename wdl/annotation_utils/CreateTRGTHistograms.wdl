version 1.0

import "../utils/Helpers.wdl"
import "../utils/Structs.wdl"

workflow CreateTRGTHistograms {
    input {
        File lps_tsv
        File metadata_tsv
        Array[String] contigs
        String prefix

        # Per-contig TRID metadata from TRGTLPS.vcf_trid_metadata_tsvs, index-aligned with
        # `contigs`. Without it the converter cannot resolve a compound TRID (a variation cluster
        # record whose INFO/TRID lists several LocusIds) and fails on the first one it sees, so it
        # is required for any callset genotyped against a catalog containing variation clusters.
        # Left empty for a catalog of isolated repeats only, which keeps the previous behavior.
        Array[File] vcf_trid_metadata_tsvs = []

        String stranalysis_docker
        String utils_docker

        RuntimeAttr? runtime_attr_subset_tsv
        RuntimeAttr? runtime_attr_convert
        RuntimeAttr? runtime_attr_concat
    }

    Boolean single_contig = length(contigs) == 1

    # Scatter by index rather than over `contigs` directly so each shard can pick the TRID
    # metadata TSV belonging to its own contig. The converter rejects metadata rows that no LPS
    # row claims, so a whole-genome metadata file paired with one contig's LPS rows would fail.
    scatter (i in range(length(contigs))) {
        String contig = contigs[i]

        # Declared inside the conditional so it is File? outside it, matching the task's optional
        # input. Indexing here is what requires the two arrays to be the same length and order.
        if (length(vcf_trid_metadata_tsvs) > 0) {
            File contig_trid_metadata_tsv = vcf_trid_metadata_tsvs[i]
        }

        if (!single_contig) {
            call SubsetLpsTsvToContig {
                input:
                    tsv = lps_tsv,
                    contig = contig,
                    prefix = "~{prefix}.~{contig}",
                    docker = utils_docker,
                    runtime_attr_override = runtime_attr_subset_tsv
            }
        }

        File contig_lps_tsv = select_first([SubsetLpsTsvToContig.subset_tsv, lps_tsv])

        call ConvertLPSTableToAFHistograms {
            input:
                lps_tsv = contig_lps_tsv,
                metadata_tsv = metadata_tsv,
                vcf_trid_metadata_tsv = contig_trid_metadata_tsv,
                prefix = "~{prefix}.~{contig}.af_histograms",
                docker = stranalysis_docker,
                runtime_attr_override = runtime_attr_convert
        }
    }

    if (!single_contig) {
        call Helpers.ConcatTsvs {
            input:
                tsvs = ConvertLPSTableToAFHistograms.af_histograms_tsv,
                sort_output = false,
                preserve_header = true,
                compressed_tsvs = true,
                prefix = "~{prefix}.af_histograms",
                docker = utils_docker,
                runtime_attr_override = runtime_attr_concat
        }
    }

    output {
        File trgt_histograms_tsv = select_first([ConcatTsvs.concatenated_tsv, ConvertLPSTableToAFHistograms.af_histograms_tsv[0]])
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
        set -eou pipefail

        zcat -f ~{tsv} \
            | awk -F'\t' -v contig="~{contig}" '
                BEGIN {
                    if (contig ~ /^chr/) contig = substr(contig, 4)
                }

                {
                    split($1, arr, ",")
                    # A self-named variation cluster is "VC:{chrom}:{start}-{end}", so its chrom
                    # is the second colon-delimited field. Splitting such a TRID on "-" like an
                    # ordinary LocusId would yield "VC:{chrom}:{start}", matching no contig and
                    # silently dropping the row from every shard.
                    if (arr[1] ~ /^VC:/) {
                        split(arr[1], vcparts, ":")
                        chr = vcparts[2]
                    } else {
                        split(arr[1], parts, "-")
                        chr = parts[1]
                    }
                    if (chr ~ /^chr/) chr = substr(chr, 4)
                    if (chr ~ /^trid/) header = $0
                    else if (chr == contig) data[++n] = $0
                }

                END {
                    print header
                    for (i = 1; i <= n; i++) print data[i]
                }
            ' \
            | gzip -c \
            > ~{prefix}.tsv.gz
    >>>

    output {
        File subset_tsv = "~{prefix}.tsv.gz"
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

task ConvertLPSTableToAFHistograms {
    input {
        File lps_tsv
        File metadata_tsv
        File? vcf_trid_metadata_tsv
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -eou pipefail

        python3 -m str_analysis.convert_multisample_LPS_table_to_allele_frequency_histograms \
            --input-table ~{lps_tsv} \
            --sample-metadata-tsv ~{metadata_tsv} \
            ~{if defined(vcf_trid_metadata_tsv) then "--vcf-trid-metadata-tsv " + vcf_trid_metadata_tsv else ""} \
            --output-format TSV \
            --stratify-by-population \
            --stratify-by-sex
        
        mv "$(dirname ~{lps_tsv})"/*.per_locus_and_motif.*.tsv.gz ~{prefix}.tsv.gz
    >>>

    output {
        File af_histograms_tsv = "~{prefix}.tsv.gz"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 8,
        disk_gb: 2 * ceil(size(lps_tsv, "GB")) + 20,
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
