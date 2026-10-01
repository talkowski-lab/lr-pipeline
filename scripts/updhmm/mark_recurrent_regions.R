#!/usr/bin/env Rscript
# mark_recurrent_regions.R — flag UPD calls that recur across unrelated trios in the cohort.
#
# Real germline UPD is essentially never recurrent across unrelated individuals. A region
# called in several different (unrelated) trios is the classic signature of a technical
# artifact -- a segmental duplication, centromeric/low-complexity region, or systematic
# mapping/genotyping bias -- rather than true UPD. UPDhmm ships identifyRecurrentRegions()/
# markRecurrentRegions() specifically for this: they cluster overlapping events across the
# whole cohort and flag any cluster supported by at least min_support distinct probands.
#
# This does not drop rows -- it only adds Recurrent ("Yes"/"No") and n_samples columns to
# the cohort-wide events table, so callers can filter Recurrent == "Yes" out downstream
# without losing the raw calls.
#
# Usage:
#   mark_recurrent_regions.R <cohort_events_tsv> <min_support> <error_threshold> \
#                             <max_dist> <min_overlap> <out_tsv>

suppressPackageStartupMessages({
    library(UPDhmm)
    library(GenomicRanges)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 6) {
    stop("usage: mark_recurrent_regions.R <cohort_events_tsv> <min_support> <error_threshold> <max_dist> <min_overlap> <out_tsv>")
}
in_tsv          <- args[[1]]
min_support     <- as.integer(args[[2]])
error_threshold <- as.numeric(args[[3]])
max_dist        <- as.numeric(args[[4]])
min_overlap     <- as.numeric(args[[5]])
out_tsv         <- args[[6]]

df <- read.delim(in_tsv, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE)

if (nrow(df) == 0) {
    df$Recurrent <- character(0)
    df$n_samples <- integer(0)
    write.table(df, file = out_tsv, sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
    cat("mark_recurrent_regions.R: 0 events; wrote empty annotated table.\n")
    quit(status = 0)
}

# Work around a bug in UPDhmm 1.8.0: identifyRecurrentRegionsByChr() detects whether the
# input has 'n_mendelian_error' or 'total_mendelian_error' (our collapseEvents()-derived
# table only has the latter), but its actual error_threshold filter hard-codes
# mcols(gr)$n_mendelian_error regardless -- so on a 'total_mendelian_error'-only table that
# column is NULL, `any(NULL < error_threshold)` is FALSE, and it silently returns NULL for
# every chromosome (recurrence detection never fires). Alias the column (on a copy, so the
# alias doesn't leak into the written output) so the intended behavior (documented and
# supported per its own error_threshold parameter) actually runs.
df_for_clustering <- df
if (!("n_mendelian_error" %in% names(df)) && ("total_mendelian_error" %in% names(df))) {
    df_for_clustering$n_mendelian_error <- df_for_clustering$total_mendelian_error
}

# identifyRecurrentRegions() looks for 'ID' by default; our cohort table's ID column holds
# each event's proband id (set in updhmm.R), which is exactly the per-sample identifier
# recurrence should be counted over.
recurrent_gr <- identifyRecurrentRegions(
    df_for_clustering,
    ID_col = "ID",
    error_threshold = error_threshold,
    min_support = min_support,
    max_dist = max_dist
)
if (is.null(recurrent_gr)) recurrent_gr <- GRanges()

out <- markRecurrentRegions(df, recurrent_gr, min_overlap = min_overlap)

write.table(out, file = out_tsv, sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)

n_recurrent <- sum(out$Recurrent == "Yes")
cat(sprintf(
    "mark_recurrent_regions.R: flagged %d/%d event(s) as recurrent (min_support=%d, %d distinct region(s))\n",
    n_recurrent, nrow(out), min_support, length(recurrent_gr)
))
