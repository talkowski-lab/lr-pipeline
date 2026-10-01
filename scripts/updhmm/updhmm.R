#!/usr/bin/env Rscript
# updhmm.R — run the Bioconductor UPDhmm HMM on a single trio VCF and write a
# tab-separated table of uniparental-disomy (UPD) events.
#
# Usage:
#   updhmm.R <trio_vcf> <proband_id> <mother_id> <father_id> <family_id> <out_tsv> <ncpu> \
#            [min_mendelian_error] [min_size]
#
# The trio VCF is a merged, biallelic-SNP multisample VCF containing exactly the proband,
# mother and father (autosomes only), with GT plus DP/AD where available. Sample columns
# must be named with the ids passed here so vcfCheck() can locate each member.
#
# min_mendelian_error/min_size are collapseEvents()'s own confidence thresholds (an event
# must have MORE Mendelian errors than min_mendelian_error and span MORE than min_size bp to
# be retained); they default to UPDhmm's own package defaults (2, 500000) if omitted.

suppressPackageStartupMessages({
    library(VariantAnnotation)
    library(UPDhmm)
    library(BiocParallel)
    library(Rsamtools)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 7) {
    stop("usage: updhmm.R <trio_vcf> <proband_id> <mother_id> <father_id> <family_id> <out_tsv> <ncpu> [min_mendelian_error] [min_size]")
}
vcf_file  <- args[[1]]
proband   <- args[[2]]
mother    <- args[[3]]
father    <- args[[4]]
family_id <- args[[5]]
out_tsv   <- args[[6]]
ncpu      <- suppressWarnings(as.integer(args[[7]]))
min_me    <- if (length(args) >= 8) suppressWarnings(as.integer(args[[8]])) else 2L
min_size  <- if (length(args) >= 9) suppressWarnings(as.numeric(args[[9]])) else 500e3
if (is.na(min_me)) min_me <- 2L
if (is.na(min_size)) min_size <- 500e3

# Canonical output schema for UPDhmm 1.8.0 (Bioconductor 3.23, matches the
# bioconductor_docker:RELEASE_3_23 base image below). collapseEvents() -- the last step
# in this pipeline -- returns a *different*, coarser schema than calculateEvents():
# raw per-block calls (ID, chromosome, start, end, group, n_snps, ratio_*,
# n_mendelian_error) are filtered (n_mendelian_error > 2 and span > 500kb, both
# collapseEvents defaults) and rolled up per (ID, chromosome, group). We pin the
# collapsed schema so every trio's TSV shares an identical header, which the
# downstream ConcatTsvs (header taken from the first shard) relies on.
base_cols <- c("ID", "chromosome", "start", "end", "group", "n_events",
               "total_mendelian_error", "total_size", "total_snps", "prop_covered",
               "ratio_father", "ratio_mother", "ratio_proband", "collapsed_events")

bp <- if (!is.na(ncpu) && ncpu > 1) MulticoreParam(workers = ncpu) else SerialParam()

# Process one chromosome at a time so peak memory is bounded by the largest
# single contig rather than the whole genome. UPD events never span chromosomes
# (the HMM runs along one contig), so per-chromosome results are identical to a
# genome-wide run once concatenated. Only GT/DP/AD are read, and the trio VCF is
# tabix-indexed, so each ranged read touches only that contig's records.
contigs <- seqnamesTabix(vcf_file)
seq_len <- seqlengths(seqinfo(scanVcfHeader(vcf_file)))

# Rsamtools/htslib's tabix range query has a hard ceiling of 2^29 - 1 (536,870,911);
# .Machine$integer.max exceeds that and errors out ("'end' must be <= 536870912").
# Fall back to a value comfortably above the largest human chromosome (~249 Mb) but
# safely under that ceiling for VCFs whose header omits contig length=.
tabix_range_max <- 536870911L
fallback_end <- 260000000L

per_chrom <- list()
for (chr in contigs) {
    end <- if (chr %in% names(seq_len) && !is.na(seq_len[[chr]])) {
        min(seq_len[[chr]], tabix_range_max)
    } else {
        fallback_end
    }
    which_gr <- GRanges(chr, IRanges(1L, end))
    # Request DP/AD alongside GT for add_ratios' depth-ratio computation. readVcf tolerates
    # requesting geno fields absent from the header (warns and simply omits them), so this
    # is safe whether or not a given trio's VCF actually carries DP/AD.
    param <- ScanVcfParam(info = NA, geno = c("GT", "DP", "AD"), which = which_gr)
    vcf <- readVcf(vcf_file, genome = "hg38", param = param)
    if (nrow(vcf) == 0) { rm(vcf); next }

    vcf <- vcfCheck(vcf, proband = proband, mother = mother, father = father)
    ev <- calculateEvents(vcf, add_ratios = TRUE, BPPARAM = bp)
    if (!is.null(ev) && nrow(ev) > 0) per_chrom[[chr]] <- ev

    rm(vcf)
    gc(verbose = FALSE)
}

if (length(per_chrom) > 0) {
    events <- collapseEvents(do.call(rbind, per_chrom), min_ME = min_me, min_size = min_size)
} else {
    events <- data.frame(setNames(rep(list(character(0)), length(base_cols)), base_cols),
                         stringsAsFactors = FALSE)
}

# Guarantee a fixed column set/order regardless of what calculateEvents/collapseEvents
# emit for this input (e.g. zero-event trios), so headers are identical across shards.
for (col in base_cols) {
    if (!(col %in% colnames(events))) events[[col]] <- rep(NA, nrow(events))
}
extra <- sort(setdiff(colnames(events), base_cols))
events <- events[, c(base_cols, extra), drop = FALSE]

# rep(..., nrow(events)) rather than a bare scalar: when collapseEvents() filters every
# event out (0 rows -- a valid, common outcome, not just the "no chromosomes" case above),
# a scalar here would build a 1-row frame and cbind against a 0-row `events` would error
# ("differing number of rows: 1, 0").
out <- cbind(
    data.frame(family_id = rep(family_id, nrow(events)),
               proband_id = rep(proband, nrow(events)),
               stringsAsFactors = FALSE),
    events
)

write.table(out, file = out_tsv, sep = "\t", quote = FALSE,
            row.names = FALSE, col.names = TRUE)

cat(sprintf("updhmm.R: wrote %d UPD event(s) for family=%s proband=%s -> %s\n",
            nrow(out), family_id, proband, out_tsv))
