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

# Read only GT to keep memory bounded on whole-genome trios.
param <- ScanVcfParam(info = NA, geno = "GT")
vcf <- readVcf(vcf_file, genome = "hg38", param = param)

vcf <- vcfCheck(vcf, proband = proband, mother = mother, father = father)

bp <- if (!is.na(ncpu) && ncpu > 1) MulticoreParam(workers = ncpu) else SerialParam()
events <- calculateEvents(vcf, add_ratios = TRUE, BPPARAM = bp)

if (!is.null(events) && nrow(events) > 0) {
    events <- collapseEvents(events)
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
