#!/usr/bin/env Rscript
# Tandem-duplication stage for insertions without a repeat annotation (repeat_class "non-repeat (next stage)" in
# build_repeat_classes.py), from the remap step-3 hits (replay_remap_step3.R).
#   1. tandem duplication: an alignment on the insertion contig, + strand, within --adjacent bp of the site (distance 0
#      when the interval touches or contains the site, else the gap to its nearer end) that covers >= --min-perc % of
#      the insertion (remap itself uses the insertion length as the distance and 80%);
#      inverted tandem duplication: the same on the - strand
#   2. rescued from remap partial: the insertion point lies inside the duplicated unit, so the inserted sequence is
#      [unit right of the site][unit left of the site] and aligns as two + pieces meeting at the site: one starts the
#      query (within --tol) and maps starting at the site, one ends the query and maps ending at the site
#      ("tandem duplication, insertion point inside the unit"; insertion length within --len-tol of the unit, else
#      "multi-copy / partial unit")
# Insertions with no call have td_call NA (empty in the TSV) and go to the next stage.
# For every call: the duplicated unit (hg38 interval), its length, insertion length / unit length, and where the
# insertion point sits relative to the unit: offset = bases of the unit left of the site, position = offset / unit
# length (0 = at the left edge, 1 = at the right edge, in between = inside the unit).
suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
})
opt <- parse_args(OptionParser(
  usage = "%prog --classes <repeat_classes.tsv> --hits <step3_hits.tsv.gz> --decisions <step3_decisions.tsv.gz> --out-prefix <prefix>",
  description = "Call tandem duplications (>= 90% adjacent alignment) and rotated tandem duplications among non-repeat insertions.",
  option_list = list(
    make_option("--classes", help = "build_repeat_classes.py output"),
    make_option("--hits", help = "Per-hit table from replay_remap_step3.R"),
    make_option("--decisions", help = "Per-insertion table from replay_remap_step3.R (CHROM, POS, alt_len, remap_class)"),
    make_option("--select", default = "non-repeat (next stage)", help = "repeat_class to analyse [default %default]"),
    make_option("--min-perc", dest = "min_perc", type = "double", default = 90, help = "Min %% of the insertion covered [default %default]"),
    make_option("--adjacent", type = "integer", default = 5, help = "Max site-to-copy distance for a tandem duplication, bp [default %default]"),
    make_option("--tol", type = "integer", default = 200, help = "bp tolerance for the rotated-unit pieces [default %default]"),
    make_option("--len-tol", dest = "len_tol", type = "double", default = 0.15, help = "|insertion / unit - 1| for a single-unit rotated dup [default %default]"),
    make_option("--out-prefix", dest = "out_prefix", help = "Output prefix")
  )
))
for (a in c("classes", "hits", "decisions", "out_prefix")) if (is.null(opt[[a]])) stop("--", a, " is required")

cl <- fread(opt$classes, na.strings = c("", "NA"))[repeat_class == opt$select, .(ID)]
d <- fread(opt$decisions)[ID %in% cl$ID, .(ID, CHROM, POS, alt_len, remap_class)]
h <- merge(fread(opt$hits)[ID %in% cl$ID], d, by = "ID")
h[, `:=`(rchrom = sub(":.*", "", coords), rs = as.numeric(sub("^[^:]+:([0-9]+)-.*", "\\1", coords)), re = as.numeric(sub(".*-", "", coords)))]
h[, site_dist := fifelse(rchrom == CHROM & rs - 1 <= POS & POS <= re, 0, chrom_dist)]
h[, adjacent := rchrom == CHROM & site_dist <= opt$adjacent]

res <- h[, {
  ad_p <- which(adjacent & strand == "+"); ad_m <- which(adjacent & strand == "-")
  bp <- if (length(ad_p)) ad_p[which.max(perc[ad_p])] else NA_integer_
  bm <- if (length(ad_m)) ad_m[which.max(perc[ad_m])] else NA_integer_
  # rotated: A = starts the query and maps starting at the site; B = ends the query and maps ending at the site
  A <- which(strand == "+" & rchrom == CHROM & q_start <= opt$tol + 1 & abs(rs - (POS + 1)) <= opt$tol)
  B <- which(strand == "+" & rchrom == CHROM & q_end >= qlen - opt$tol & abs(re - POS) <= opt$tol)
  rot <- length(A) > 0 && length(B) > 0 && max(q_end[A]) >= min(q_start[B]) - opt$tol
  .(best_adj_plus = if (is.na(bp)) 0 else perc[bp], best_adj_minus = if (is.na(bm)) 0 else perc[bm],
    plus_unit = if (is.na(bp)) NA_character_ else coords[bp], minus_unit = if (is.na(bm)) NA_character_ else coords[bm],
    plus_ident = if (is.na(bp)) NA_real_ else identity[bp], minus_ident = if (is.na(bm)) NA_real_ else identity[bm],
    rot_ident = if (rot) min(max(identity[A]), max(identity[B])) else NA_real_,
    rotated = rot, rot_start = if (rot) min(rs[B]) else NA_real_, rot_end = if (rot) max(re[A]) else NA_real_)
}, by = .(ID, CHROM, POS, alt_len, remap_class)]
res <- merge(d[, .(ID)], res, by = "ID", all.x = TRUE)
res <- merge(res, d, by = "ID", suffixes = c(".x", ""))[, c("CHROM.x", "POS.x", "alt_len.x", "remap_class.x") := NULL]
for (k in c("best_adj_plus", "best_adj_minus")) set(res, which(is.na(res[[k]])), k, 0)
res[is.na(rotated), rotated := FALSE]
res[, len := alt_len - 1L]

res[, td_call := fcase(
  best_adj_plus >= opt$min_perc, "tandem duplication",
  best_adj_minus >= opt$min_perc, "inverted tandem duplication",
  rotated & abs(len / (rot_end - rot_start + 1) - 1) <= opt$len_tol, "tandem duplication, insertion point inside the unit",
  rotated, "rotated pieces, multi-copy / partial unit",
  best_adj_plus >= 80, "remap tandem at 80-90% (not called)",
  default = NA_character_)]
res[, unit := fcase(td_call == "tandem duplication", plus_unit, td_call == "inverted tandem duplication", minus_unit,
                    rotated, sprintf("%s:%d-%d", CHROM, as.integer(rot_start), as.integer(rot_end)), default = NA_character_)]
res[, unit_strand := fcase(td_call == "inverted tandem duplication", "-", !is.na(unit), "+", default = NA_character_)]
res[, unit_identity := fcase(td_call == "tandem duplication", plus_ident, td_call == "inverted tandem duplication", minus_ident,
                             rotated, rot_ident, default = NA_real_)]
res[, `:=`(us = as.numeric(sub("^[^:]+:([0-9]+)-.*", "\\1", unit)), ue = as.numeric(sub(".*-", "", unit)))]
res[, unit_len := ue - us + 1]
res[, len_over_unit := round(len / unit_len, 3)]
res[, ip_offset := POS - us + 1]
res[, ip_position := round(ip_offset / unit_len, 3)]
res[, ip_where := fcase(is.na(unit), NA_character_,
                        ip_position <= 0.05, "at / left of the unit's left edge",
                        ip_position >= 0.95, "at / right of the unit's right edge",
                        default = "inside the unit")]
res[, c("us", "ue") := NULL]
setorder(res, td_call, -len)
fwrite(res[, .(ID, CHROM, POS, len, remap_class, td_call, unit, unit_strand, unit_identity, unit_len, len_over_unit, ip_offset, ip_position, ip_where,
               best_adj_plus, best_adj_minus, rotated)], paste0(opt$out_prefix, ".tandem_dups.tsv"), sep = "\t")

o <- function(x, name) { fwrite(x, sprintf("%s.%s.tsv", opt$out_prefix, name), sep = "\t"); cat("\n==", name, "\n"); print(x) }
cat(sprintf("insertions analysed (%s): %d\n", opt$select, nrow(res)))
o(res[, .N, by = .(td_call = fifelse(is.na(td_call), "(next stage)", td_call), remap_class)][order(td_call, remap_class)], "tandem_dup_counts")
o(res[!is.na(td_call) & !is.na(unit), .(n = .N, med_len = as.numeric(median(len)), med_len_over_unit = median(len_over_unit),
                                      left_edge = sum(ip_where %like% "left"), inside = sum(ip_where == "inside the unit"),
                                      right_edge = sum(ip_where %like% "right")), by = td_call], "tandem_dup_insertion_point")
