version 1.0

import "../wdl/utils/Structs.wdl"

# Link methylation outlier haplotypes to rare variants within `window` bp of the outlier region.
#
# `outlier_calls` and `mean_matrix` are RegionHaplotypeMethylation outputs (genome-wide); `vcfs` / `vcf_idxs` /
# `contigs` are parallel per-contig arrays of the phased cohort VCF. All logic lives in outlier_variant_linkage.py
# (passed in as `analysis_script`). Per contig (scattered):
#   1. MakeWindows: merged BED of [start - window, end + window) around the contig's outlier regions.
#   2. ExtractWindowGenotypes: bcftools view of the windows (records with INFO/AF <= prefilter_max_af in any ALT
#      allele, to drop common variants early), then bcftools query of position / ID / FILTER / allele_type /
#      allele_length / cadd_phred and per-sample GT:PS.
#   3. LinkOutlierVariants: rare alleles (AF <= max_af among retained samples = VCF and methylation samples minus
#      excluded_samples) carried by each outlier sample within the window, with zygosity, phase and whether the ALT
#      lies on the outlier haplotype; plus per outlier call the rare-allele burden at each of `distances` and its
#      empirical p-value against all other retained samples.
#   4. Concatenation of the per-contig tables.
# Contigs without outlier calls produce header-only tables.

workflow OutlierVariantLinkage {
    input {
        Array[File] vcfs
        Array[File] vcf_idxs
        Array[String] contigs
        File outlier_calls
        File mean_matrix
        String prefix

        File? excluded_samples
        Int window = 1000000
        Float prefilter_max_af = 0.05
        Float max_af = 0.01
        String distances = "10000 100000 1000000"

        File analysis_script

        String python_docker
        String bcftools_docker

        RuntimeAttr? runtime_attr_make_windows
        RuntimeAttr? runtime_attr_extract_window_genotypes
        RuntimeAttr? runtime_attr_link_outlier_variants
        RuntimeAttr? runtime_attr_concat_tables
    }

    scatter (i in range(length(contigs))) {
        call MakeWindows {
            input:
                outlier_calls = outlier_calls,
                contig = contigs[i],
                window = window,
                analysis_script = analysis_script,
                prefix = "~{prefix}.~{contigs[i]}",
                docker = python_docker,
                runtime_attr_override = runtime_attr_make_windows
        }

        call ExtractWindowGenotypes {
            input:
                vcf = vcfs[i],
                vcf_idx = vcf_idxs[i],
                windows_bed = MakeWindows.windows_bed,
                prefilter_max_af = prefilter_max_af,
                prefix = "~{prefix}.~{contigs[i]}",
                docker = bcftools_docker,
                runtime_attr_override = runtime_attr_extract_window_genotypes
        }

        call LinkOutlierVariants {
            input:
                outlier_calls = outlier_calls,
                mean_matrix = mean_matrix,
                genotypes = ExtractWindowGenotypes.genotypes,
                vcf_samples = ExtractWindowGenotypes.vcf_samples,
                excluded_samples = excluded_samples,
                contig = contigs[i],
                window = window,
                max_af = max_af,
                distances = distances,
                analysis_script = analysis_script,
                prefix = "~{prefix}.~{contigs[i]}",
                docker = python_docker,
                runtime_attr_override = runtime_attr_link_outlier_variants
        }
    }

    call ConcatTables as ConcatOutlierVariants {
        input:
            tables = LinkOutlierVariants.outlier_variants,
            prefix = "~{prefix}.outlier_variants",
            docker = python_docker,
            runtime_attr_override = runtime_attr_concat_tables
    }

    call ConcatTables as ConcatOutlierBurden {
        input:
            tables = LinkOutlierVariants.outlier_burden,
            prefix = "~{prefix}.outlier_burden",
            docker = python_docker,
            runtime_attr_override = runtime_attr_concat_tables
    }

    output {
        File outlier_variants = ConcatOutlierVariants.merged
        File outlier_burden = ConcatOutlierBurden.merged
    }
}

task MakeWindows {
    input {
        File outlier_calls
        String contig
        Int window
        File analysis_script
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        python ~{analysis_script} windows \
            --outlier-calls ~{outlier_calls} \
            --contig ~{contig} \
            --window ~{window} \
            --out-bed ~{prefix}.windows.bed
    >>>

    output {
        File windows_bed = "~{prefix}.windows.bed"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 4,
        disk_gb: 2 * ceil(size(outlier_calls, "GB")) + 10,
        boot_disk_gb: 10,
        preemptible_tries: 2,
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

task ExtractWindowGenotypes {
    input {
        File vcf
        File vcf_idx
        File windows_bed
        Float prefilter_max_af
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        bcftools query -l ~{vcf} > ~{prefix}.vcf_samples.txt
        if [ -s ~{windows_bed} ]; then
            bcftools view \
                -R ~{windows_bed} \
                --regions-overlap pos \
                -i 'MIN(INFO/AF)<=~{prefilter_max_af}' \
                -Ou \
                ~{vcf} \
            | bcftools query \
                -f '%CHROM\t%POS\t%END\t%ID\t%FILTER\t%INFO/allele_type\t%INFO/allele_length\t%INFO/cadd_phred[\t%GT:%PS]\n' \
            | bgzip > ~{prefix}.window_genotypes.tsv.gz
        else
            bgzip < /dev/null > ~{prefix}.window_genotypes.tsv.gz
        fi
    >>>

    output {
        File genotypes = "~{prefix}.window_genotypes.tsv.gz"
        File vcf_samples = "~{prefix}.vcf_samples.txt"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 4,
        disk_gb: 3 * ceil(size(vcf, "GB")) + 10,
        boot_disk_gb: 10,
        preemptible_tries: 2,
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

task LinkOutlierVariants {
    input {
        File outlier_calls
        File mean_matrix
        File genotypes
        File vcf_samples
        File? excluded_samples
        String contig
        Int window
        Float max_af
        String distances
        File analysis_script
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        python ~{analysis_script} link \
            --outlier-calls ~{outlier_calls} \
            --mean-matrix ~{mean_matrix} \
            --genotypes ~{genotypes} \
            --vcf-samples ~{vcf_samples} \
            ~{"--excluded-samples " + excluded_samples} \
            --contig ~{contig} \
            --window ~{window} \
            --max-af ~{max_af} \
            --distances ~{distances} \
            --prefix ~{prefix}
    >>>

    output {
        File outlier_variants = "~{prefix}.outlier_variants.tsv.gz"
        File outlier_burden = "~{prefix}.outlier_burden.tsv.gz"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 2,
        mem_gb: ceil(30 * size(genotypes, "GB")) + 8,
        disk_gb: 3 * ceil(size([genotypes, mean_matrix], "GB")) + 10,
        boot_disk_gb: 10,
        preemptible_tries: 2,
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

task ConcatTables {
    input {
        Array[File] tables
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        first=1
        for f in ~{sep=" " tables}; do
            if [ "$first" -eq 1 ]; then
                gzip -dc "$f"
                first=0
            else
                gzip -dc "$f" | tail -n +2
            fi
        done | gzip > ~{prefix}.tsv.gz
    >>>

    output {
        File merged = "~{prefix}.tsv.gz"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 2,
        disk_gb: 3 * ceil(size(tables, "GB")) + 10,
        boot_disk_gb: 10,
        preemptible_tries: 2,
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
