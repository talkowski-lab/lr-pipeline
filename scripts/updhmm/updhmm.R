#!/usr/bin/env Rscript
# updhmm.R — run the Bioconductor UPDhmm HMM on a single trio VCF and write a
# tab-separated table of uniparental-disomy (UPD) events.
#
# Usage:
#   updhmm.R <trio_vcf> <proband_id> <mother_id> <father_id> <family_id> <out_tsv> <ncpu>
#
# The trio VCF is a merged, biallelic-SNP, GT-only multisample VCF containing exactly
# the proband, mother and father (autosomes only). Sample columns must be named with
# the ids passed here so vcfCheck() can locate each member.

suppressPackageStartupMessages({
    library(VariantAnnotation)
    library(UPDhmm)
    library(BiocParallel)
    library(Rsamtools)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 7) {
    stop("usage: updhmm.R <trio_vcf> <proband_id> <mother_id> <father_id> <family_id> <out_tsv> <ncpu>")
}
vcf_file  <- args[[1]]
proband   <- args[[2]]
mother    <- args[[3]]
father    <- args[[4]]
family_id <- args[[5]]
out_tsv   <- args[[6]]
ncpu      <- suppressWarnings(as.integer(args[[7]]))

# Canonical calculateEvents() output schema. The ratio_* columns are always present
# (filled with NA when add_ratios = FALSE). We pin this order so every trio's TSV shares
# an identical header, which the downstream ConcatTsvs (header taken from the first
# shard) relies on.
base_cols <- c("ID", "chromosome", "start", "end", "group", "n_snps",
               "ratio_father", "ratio_mother", "ratio_proband", "n_mendelian_error")

bp <- if (!is.na(ncpu) && ncpu > 1) MulticoreParam(workers = ncpu) else SerialParam()

# Process one chromosome at a time so peak memory is bounded by the largest
# single contig rather than the whole genome. UPD events never span chromosomes
# (the HMM runs along one contig), so per-chromosome results are identical to a
# genome-wide run once concatenated. Only GT is read, and the trio VCF is tabix-
# indexed, so each ranged read touches only that contig's records.
contigs <- seqnamesTabix(vcf_file)
seq_len <- seqlengths(seqinfo(scanVcfHeader(vcf_file)))

per_chrom <- list()
for (chr in contigs) {
    end <- if (chr %in% names(seq_len) && !is.na(seq_len[[chr]])) seq_len[[chr]] else .Machine$integer.max
    which_gr <- GRanges(chr, IRanges(1L, end))
    param <- ScanVcfParam(info = NA, geno = "GT", which = which_gr)
    vcf <- readVcf(vcf_file, genome = "hg38", param = param)
    if (nrow(vcf) == 0) { rm(vcf); next }

    vcf <- vcfCheck(vcf, proband = proband, mother = mother, father = father)
    ev <- calculateEvents(vcf, add_ratios = TRUE, BPPARAM = bp)
    if (!is.null(ev) && nrow(ev) > 0) per_chrom[[chr]] <- ev

    rm(vcf)
    gc(verbose = FALSE)
}

if (length(per_chrom) > 0) {
    events <- collapseEvents(do.call(rbind, per_chrom))
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

out <- cbind(
    data.frame(family_id = family_id, proband_id = proband, stringsAsFactors = FALSE),
    events
)

write.table(out, file = out_tsv, sep = "\t", quote = FALSE,
            row.names = FALSE, col.names = TRUE)

cat(sprintf("updhmm.R: wrote %d UPD event(s) for family=%s proband=%s -> %s\n",
            nrow(out), family_id, proband, out_tsv))
