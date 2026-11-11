version 1.0

# Per-sample count of SVs of one type (default DEL) whose |SVLEN| is within [min_size, max_size], across per-contig VCFs.
# 1. CountPerContig (one job per VCF): selects the SVs (INFO/SVTYPE == svtype, min_size <= |INFO/SVLEN| <= max_size,
#    FILTER PASS when pass_only) and counts, per sample, the selected SVs the sample carries (non-reference GT):
#    n_carried = n_het + n_homalt (haploid alternative calls count as n_homalt); samples carrying none are reported as 0.
# 2. MergeAndPlot: sums the per-contig counts per sample and plots the distribution of n_carried across samples
#    (histogram with median and mean, PDF in Arial when arial_ttf is given, otherwise the fontconfig Arial substitute).
# Outputs: per-sample table (sample n_het n_homalt n_carried), the selected sites (chrom pos id svtype svlen filter),
# per-contig site counts and the PDF.

import "Structs.wdl"

workflow PerSampleSvCountInSizeRange {
    input {
        Array[File] vcfs
        String prefix

        String svtype = "DEL"
        Int min_size = 400
        Int max_size = 1000
        Boolean pass_only = true

        File? arial_ttf

        String docker

        RuntimeAttr? runtime_attr_count_per_contig
        RuntimeAttr? runtime_attr_merge_and_plot
    }

    parameter_meta {
        vcfs: "Per-contig VCFs (bgzipped), all with the same samples"
        prefix: "Prefix for output files"
        svtype: "INFO/SVTYPE value to select (e.g. DEL, DUP, INS)"
        min_size: "Minimum |INFO/SVLEN| (inclusive)"
        max_size: "Maximum |INFO/SVLEN| (inclusive)"
        pass_only: "Keep only FILTER == PASS records"
        arial_ttf: "Arial TrueType font for the PDF (optional)"
        docker: "Image with bcftools, R (ggplot2, data.table) and fontconfig (e.g. gd-breakpoint-profiles)"
    }

    scatter (vcf in vcfs) {
        call CountPerContig {
            input:
                vcf = vcf,
                svtype = svtype,
                min_size = min_size,
                max_size = max_size,
                pass_only = pass_only,
                prefix = prefix + "." + basename(vcf, ".vcf.gz"),
                docker = docker,
                runtime_attr_override = runtime_attr_count_per_contig
        }
    }

    call MergeAndPlot {
        input:
            counts = CountPerContig.counts,
            sites = CountPerContig.sites,
            svtype = svtype,
            min_size = min_size,
            max_size = max_size,
            pass_only = pass_only,
            arial_ttf = arial_ttf,
            prefix = prefix,
            docker = docker,
            runtime_attr_override = runtime_attr_merge_and_plot
    }

    output {
        File per_sample_counts = MergeAndPlot.per_sample_counts
        File selected_sites = MergeAndPlot.selected_sites
        File sites_per_contig = MergeAndPlot.sites_per_contig
        File count_distribution_pdf = MergeAndPlot.count_distribution_pdf
    }
}

task CountPerContig {
    input {
        File vcf
        String svtype
        Int min_size
        Int max_size
        Boolean pass_only
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        bcftools view ~{if pass_only then "-f PASS" else ""} -Oz -o selected.vcf.gz \
            -i 'INFO/SVTYPE="~{svtype}" && ((INFO/SVLEN>=~{min_size} && INFO/SVLEN<=~{max_size}) || (INFO/SVLEN<=-~{min_size} && INFO/SVLEN>=-~{max_size}))' \
            ~{vcf}
        bcftools query -f '%CHROM\t%POS\t%ID\t%INFO/SVTYPE\t%INFO/SVLEN\t%FILTER\n' selected.vcf.gz > ~{prefix}.sites.tsv
        echo "~{basename(vcf)}: $(wc -l < ~{prefix}.sites.tsv) selected SVs" >&2
        bcftools query -l selected.vcf.gz > samples.txt
        bcftools query -i 'GT="alt"' -f '[%SAMPLE\t%GT\n]' selected.vcf.gz \
            | awk -F'\t' -v OFS='\t' '
                FNR == NR {order[++n] = $1; het[$1] = 0; hom[$1] = 0; next}
                $2 ~ /[1-9]/ {if ($2 ~ /0/) het[$1]++; else hom[$1]++}
                END {print "sample", "n_het", "n_homalt", "n_carried"; for (i = 1; i <= n; i++) {s = order[i]; print s, het[s], hom[s], het[s] + hom[s]}}
            ' samples.txt - > ~{prefix}.counts.tsv
    >>>

    output {
        File counts = "~{prefix}.counts.tsv"
        File sites = "~{prefix}.sites.tsv"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 4,
        disk_gb: ceil(1.5 * size(vcf, "GB")) + 10,
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

task MergeAndPlot {
    input {
        Array[File] counts
        Array[File] sites
        String svtype
        Int min_size
        Int max_size
        Boolean pass_only
        File? arial_ttf
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        if [[ -n "~{arial_ttf}" ]]; then
            mkdir -p ~/.fonts
            cp ~{arial_ttf} ~/.fonts/Arial.ttf
            fc-cache -f > /dev/null
        fi
        cat ~{sep=" " sites} | sort -k1,1V -k2,2n > ~{prefix}.selected_sites.tsv
        printf 'chrom\tpos\tid\tsvtype\tsvlen\tfilter\n' | cat - ~{prefix}.selected_sites.tsv | gzip -c > ~{prefix}.selected_sites.tsv.gz
        printf 'chrom\tn_sites\n' > ~{prefix}.sites_per_contig.tsv
        cut -f1 ~{prefix}.selected_sites.tsv | uniq -c | awk -v OFS='\t' '{print $2, $1}' >> ~{prefix}.sites_per_contig.tsv
        N_SITES=$(wc -l < ~{prefix}.selected_sites.tsv)

        Rscript - ~{prefix} "~{svtype}" ~{min_size} ~{max_size} ~{if pass_only then "PASS" else "all"} "${N_SITES}" ~{sep=" " counts} <<'RS'
suppressPackageStartupMessages({library(data.table); library(ggplot2)})
a <- commandArgs(trailingOnly = TRUE)
prefix <- a[1]; svtype <- a[2]; min_size <- as.integer(a[3]); max_size <- as.integer(a[4]); filt <- a[5]
n_sites <- as.integer(a[6]); files <- a[-(1:6)]
d <- rbindlist(lapply(files, fread, colClasses = c(sample = "character")))
d <- d[, .(n_het = sum(n_het), n_homalt = sum(n_homalt), n_carried = sum(n_carried)), by = sample]
fwrite(d, paste0(prefix, ".per_sample_counts.tsv.gz"), sep = "\t")

fams <- unlist(strsplit(tryCatch(system("fc-list : family", intern = TRUE), error = function(e) "Arial"), ","))
font_family <- if ("Arial" %in% trimws(fams)) "Arial" else "sans"
x <- d$n_carried
med <- median(x); mu <- mean(x)
bw <- max(1, ceiling(diff(range(x)) / 60))
stats <- data.table(v = c(med, mu), lab = c(sprintf("median %g", med), sprintf("mean %.1f", mu)), lt = c("dashed", "dotted"))
p <- ggplot(d, aes(n_carried)) +
  geom_histogram(binwidth = bw, boundary = -0.5, fill = "grey60", colour = "white", linewidth = 0.2) +
  geom_vline(data = stats, aes(xintercept = v, linetype = lt), colour = "firebrick", linewidth = 1) +
  scale_linetype_identity() +
  geom_text(data = stats, aes(x = v, y = Inf, label = lab), vjust = c(1.5, 3.2), hjust = -0.05, size = 7.76,
            colour = "firebrick", family = font_family) +
  labs(x = sprintf("%s (%d to %d bp) per sample", svtype, min_size, max_size), y = "Samples",
       title = sprintf("%s: %s SVs, %d to %d bp (%s)", prefix, svtype, min_size, max_size, filt),
       subtitle = sprintf("%s samples, %s selected SVs; per-sample count of carried SVs (het + homalt)",
                          format(nrow(d), big.mark = ","), format(n_sites, big.mark = ","))) +
  theme_bw(base_size = 22, base_family = font_family) +
  theme(panel.grid.minor = element_blank())
out <- paste0(prefix, ".count_distribution.pdf")
if (capabilities("cairo")) cairo_pdf(out, width = 14, height = 9, family = font_family) else pdf(out, width = 14, height = 9)
print(p)
invisible(dev.off())
RS
    >>>

    output {
        File per_sample_counts = "~{prefix}.per_sample_counts.tsv.gz"
        File selected_sites = "~{prefix}.selected_sites.tsv.gz"
        File sites_per_contig = "~{prefix}.sites_per_contig.tsv"
        File count_distribution_pdf = "~{prefix}.count_distribution.pdf"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 8,
        disk_gb: 20,
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
