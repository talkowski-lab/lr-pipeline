#!/usr/bin/env Rscript
# Duplication stage for the non-repeat insertions: one call per insertion from
#   call_tandem_dups.R            tandem duplication (>= 90% adjacent), inverted tandem duplication, tandem duplication
#                                 with the insertion point inside the unit
#   plot_interspersed_coverage.R  interspersed duplication = the best non-adjacent alignment (the interspersed origin)
#                                 covers >= --min-cov of the inserted sequence; the source location is kept in the detail
#                                 (same chromosome + distance, other chromosome, unplaced / alt contig)
# Everything else (incl. remap tandem at 80-90%, rotated multi-copy pieces, partial, novel) goes to the next stage.
suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
})
opt <- parse_args(OptionParser(
  usage = "%prog --tandem <tandem_dups.tsv> --interspersed <interspersed_dup_coverage.tsv> --out-prefix <prefix>",
  description = "Combine tandem and interspersed duplication calls for non-repeat insertions.",
  option_list = list(
    make_option("--tandem", help = "call_tandem_dups.R per-insertion output"),
    make_option("--interspersed", help = "plot_interspersed_coverage.R per-insertion output"),
    make_option("--min-cov", dest = "min_cov", type = "double", default = 0.9, help = "Min fraction covered by the interspersed origin [default %default]"),
    make_option("--out-prefix", dest = "out_prefix", help = "Output prefix")
  )
))
for (a in c("tandem", "interspersed", "out_prefix")) if (is.null(opt[[a]])) stop("--", a, " is required")
td <- fread(opt$tandem, na.strings = c("", "NA"))
it <- fread(opt$interspersed, na.strings = c("", "NA"))
x <- merge(td, it[, .(ID, best_cov, union_cov, source, source_strand, source_identity, source_mapq, location, chrom_dist)], by = "ID", all.x = TRUE)
called <- c("tandem duplication", "inverted tandem duplication", "tandem duplication, insertion point inside the unit")
x[, dup_call := fcase(td_call %in% called, td_call,
                      !is.na(best_cov) & best_cov >= opt$min_cov, "interspersed duplication",
                      default = "next stage")]
# a tandem duplication whose insertion point the caller placed inside the unit is the same event (the inserted sequence
# is a rotation of the unit): one class, with the insertion point kept in the detail
x[dup_call == "tandem duplication, insertion point inside the unit", dup_call := "tandem duplication"]
x[, dup_detail := fcase(
  dup_call %like% "^tandem|^inverted", sprintf("unit %s (%s, %s bp, identity %.3f), insertion / unit %s, insertion point %s (%.2f)",
                                               unit, unit_strand, unit_len, unit_identity, len_over_unit, ip_where, ip_position),
  dup_call == "interspersed duplication", sprintf("source %s (%s), %s%s, covers %.0f%%, identity %.3f, mapq %s", source, source_strand, location,
                                                  fifelse(location == "same chromosome", sprintf(" %s bp away", trimws(format(chrom_dist, big.mark = ","))), ""),
                                                  100 * best_cov, source_identity, source_mapq),
  default = fifelse(is.na(td_call), remap_class, td_call))]
setorder(x, dup_call, -len)
fwrite(x[, .(ID, CHROM, POS, len, remap_class, dup_call, dup_detail)], paste0(opt$out_prefix, ".dup_stage.tsv"), sep = "\t")
s <- x[, .N, by = .(dup_call, location = fifelse(dup_call == "interspersed duplication", location, ""))][order(dup_call, -N)]
fwrite(s, paste0(opt$out_prefix, ".dup_stage_counts.tsv"), sep = "\t")
print(s)
cat(sprintf("next stage: %d (%s)\n", x[dup_call == "next stage", .N],
            paste(x[dup_call == "next stage", .N, by = .(k = fifelse(is.na(td_call), remap_class, td_call))][order(-N), sprintf("%s %d", k, N)], collapse = ", ")))
