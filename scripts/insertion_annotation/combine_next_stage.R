#!/usr/bin/env Rscript
# Combine the immune-locus and mCNV annotations of the insertions left after the duplication stage.
#   mCNV            multi-copy local tandem gain (classify_partial_structure.R + multicopy_gain_table.R): the same nearby
#                   reference segment appears >= 1.5 times in the insertion; detail = segment, copies added, number of
#                   alleles at the locus (insertions whose segments overlap) and the locus copy range
#   immune          annotate_immune_loci.R call (site in an immune locus and / or sequence from one)
#   next_call       "mCNV", "mCNV in <locus>", "<immune call>", or "unresolved (next stage)"
suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
})
opt <- parse_args(OptionParser(
  usage = "%prog --ids <ids.txt> --structure <partial_structure.tsv> --mcnv <multicopy_gain.tsv> --immune <immune.tsv> --out-prefix <prefix>",
  description = "Combine mCNV and immune-locus calls for the remaining insertions.",
  option_list = list(
    make_option("--ids", help = "File of insertion IDs (the remaining set)"),
    make_option("--structure", help = "classify_partial_structure.R per-insertion output"),
    make_option("--mcnv", help = "multicopy_gain_table.R output"),
    make_option("--immune", help = "annotate_immune_loci.R output"),
    make_option("--out-prefix", dest = "out_prefix", help = "Output prefix")
  )
))
for (a in c("ids", "structure", "mcnv", "immune", "out_prefix")) if (is.null(opt[[a]])) stop("--", a, " is required")
x <- data.table(ID = readLines(opt$ids))
s <- fread(opt$structure, na.strings = c("", "NA"))[, .(ID, len, structure, remap_class = NULL)]
m <- fread(opt$mcnv, na.strings = c("", "NA"))[, .(ID, mcnv_locus = locus, mcnv_alleles = locus_n_alleles, mcnv_copy_range = locus_copy_range,
                                                  mcnv_segment = aligned_region, mcnv_segment_len = region_len, mcnv_copies_added = copies_by_len)]
im <- fread(opt$immune, na.strings = c("", "NA"))[, .(ID, immune_call)]
x <- Reduce(function(a, b) merge(a, b, by = "ID", all.x = TRUE), list(x, s, m, im))
x[, mcnv := !is.na(mcnv_segment)]
x[, next_call := fcase(mcnv & !is.na(immune_call), paste0("mCNV; ", immune_call),
                       mcnv, "mCNV",
                       !is.na(immune_call), immune_call,
                       default = "unresolved (next stage)")]
x[, next_detail := fifelse(mcnv, sprintf("segment %s (%s bp), +%s copies; %d allele(s) at the locus, copies added %s",
                                         mcnv_segment, mcnv_segment_len, mcnv_copies_added, mcnv_alleles, mcnv_copy_range), structure)]
setorder(x, next_call, -len)
fwrite(x, paste0(opt$out_prefix, ".next_stage_calls.tsv"), sep = "\t")
cnt <- x[, .N, by = .(next_call)][order(-N)]
fwrite(cnt, paste0(opt$out_prefix, ".next_stage_counts.tsv"), sep = "\t")
print(cnt)
cat(sprintf("mCNV: %d insertions at %d loci (%d loci with >= 2 alleles, %d insertions)\n", x[mcnv == TRUE, .N],
            uniqueN(x[mcnv == TRUE, mcnv_locus]), uniqueN(x[mcnv_alleles >= 2, mcnv_locus]), x[mcnv_alleles >= 2, .N]))
