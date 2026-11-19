#!/usr/bin/env Rscript
# Call T2T matches from compare_reference_coverage.R output: an insertion whose best single T2T alignment covers
# >= --min-perc % of the inserted sequence is "present in T2T-CHM13" (sequence missing from hg38 at that site), with the
# T2T location; the rest stay unresolved. Counts by size range.
suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
})
opt <- parse_args(OptionParser(
  usage = "%prog --coverage <reference_coverage.tsv> --out-prefix <prefix>",
  description = "Insertions explained by a single T2T alignment covering >= 90% of the inserted sequence.",
  option_list = list(
    make_option("--coverage", help = "compare_reference_coverage.R per-insertion output (reference B = T2T)"),
    make_option("--min-perc", dest = "min_perc", type = "double", default = 90, help = "Min %% of the insertion in one T2T alignment [default %default]"),
    make_option("--out-prefix", dest = "out_prefix", help = "Output prefix")
  )
))
for (a in c("coverage", "out_prefix")) if (is.null(opt[[a]])) stop("--", a, " is required")
x <- fread(opt$coverage, na.strings = c("", "NA"))
x[, len := as.numeric(sub(".*-", "", ID))]
x[, size := cut(len, c(0, 100, 250, 1000, 5000, Inf), right = FALSE, labels = c("50-99", "100-249", "250-999", "1-5k", ">=5k"))]
x[, t2t_chrom := sub(":.*", "", T2T_best_hit)]
x[, t2t_call := fifelse(!is.na(T2T_best_perc) & T2T_best_perc >= opt$min_perc,
                        fifelse(t2t_chrom == sub("-.*", "", ID), "present in T2T-CHM13 (same chromosome)", "present in T2T-CHM13 (other chromosome)"),
                        "unresolved")]
setorder(x, t2t_call, -len)
fwrite(x[, .(ID, len, t2t_call, T2T_best_perc, T2T_best_hit, T2T_cov, hg38_best_perc, hg38_cov)], paste0(opt$out_prefix, ".t2t_calls.tsv"), sep = "\t")
s <- dcast(x[, .N, by = .(size, t2t_call)], size ~ t2t_call, value.var = "N", fill = 0)
fwrite(s, paste0(opt$out_prefix, ".t2t_calls_by_size.tsv"), sep = "\t")
print(s)
