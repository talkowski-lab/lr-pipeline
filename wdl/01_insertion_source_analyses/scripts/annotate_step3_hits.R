#!/usr/bin/env Rscript
# Add the existing annotations to the per-hit table of the remap step-3 replay (replay_remap_step3.R).
#
# Per hit (from annotate_hit_repeats.R): fraction of the hit's reference interval in each hg38 repeat
#   category (f_STR, f_VNTR, f_MEI_*, f_satellite, f_simple_lowcomp, f_other_TE), hit_label, hit_group;
#   hit_at_svan_dup: the hit overlaps a SVAN DUP_COORD segment of its insertion;
#   hit_in_rotated_unit: the hit lies within the unit of a rotated tandem dup (find_rotated_tandem_dups.R).
# Per insertion (repeated on every hit row):
#   remap_class; VCF FILTER, AF, REGION, TRID, TR_ENVELOPED, SUB_FAMILY, allele_subtype;
#   SVAN svan_type, svan_mei_fam, DUP_COORD; merged tandem_dup status and merged_label (merge_svan_remap_mei.R);
#   rotated_tandem_dup, rotated_unit, rotated_len_over_unit; ins_best_label (label of the best-scoring hit);
#   ins_annotation: one call per insertion, first match wins:
#     VCF/override MEI > solved tandem dup > rotated tandem dup (insertion length within --rot-len-tol of the unit)
#     > STR / VNTR / satellite / simple_lowcomp by best hit > "aligns to ref <MEI/other_TE>" > "no repeat" / "no alignment".
suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
})

opt <- parse_args(OptionParser(
  usage = "%prog --hits <step3_hits.tsv.gz> --hit-repeats <hit_repeats.tsv.gz> --decisions <step3_decisions.tsv.gz> --ins-labels <insertion_repeat_labels.tsv.gz> --merged <merged_labels.tsv.gz> --rotated <rotated_tandem_dups.tsv> --info <info.tsv.gz> --out <annotated_hits.tsv.gz>",
  description = "Annotate remap step-3 hits with MEI / STR / VNTR / satellite / tandem-dup annotations.",
  option_list = list(
    make_option("--hits", help = "Per-hit table from replay_remap_step3.R"),
    make_option("--hit-repeats", dest = "hit_repeats", help = "Per-hit repeat table from annotate_hit_repeats.R"),
    make_option("--decisions", help = "Per-insertion table from replay_remap_step3.R"),
    make_option("--ins-labels", dest = "ins_labels", help = "Per-insertion table from annotate_hit_repeats.R"),
    make_option("--merged", help = "Per-insertion table from merge_svan_remap_mei.R"),
    make_option("--rotated", help = "Per-insertion table from find_rotated_tandem_dups.R"),
    make_option("--info", help = "INFO table from extract_ins_info.sh"),
    make_option("--rot-len-tol", dest = "rot_tol", type = "double", default = 0.15,
                help = "Max |insertion length / unit length - 1| for a rotated tandem dup in ins_annotation [default %default]"),
    make_option("--out", help = "Output TSV (.gz)")
  )
))
for (a in c("hits", "hit_repeats", "decisions", "ins_labels", "merged", "rotated", "info", "out")) if (is.null(opt[[a]])) stop("--", a, " is required")

keys <- c("ID", "coords", "strand", "q_start", "q_end")
h <- fread(opt$hits)
rp <- fread(opt$hit_repeats)
fcols <- grep("^f_", names(rp), value = TRUE)
h <- merge(h, rp[, c(keys, fcols, "hit_label", "hit_group"), with = FALSE], by = keys, all.x = TRUE, sort = FALSE)
if (h[is.na(hit_label), .N] > 0) stop("hits without a repeat annotation: ", h[is.na(hit_label), .N])

dec <- fread(opt$decisions)[, .(ID, remap_class)]
info <- fread(opt$info, sep = "\t", quote = "", na.strings = ".")[, .(ID, FILTER, AF, REGION, TRID, TR_ENVELOPED, SUB_FAMILY, allele_subtype)]
mg <- fread(opt$merged, na.strings = c("", "NA"))[, .(ID, svan_type, svan_mei_fam, DUP_COORD, tandem_dup, mei_family, merged_label)]
rot <- fread(opt$rotated, na.strings = c("", "NA"))[rotated == TRUE, .(ID, rotated_tandem_dup = TRUE,
                                                                    rotated_unit = sprintf("%s:%d-%d", sub("-.*", "", ID), unit_start, unit_end),
                                                                    rotated_len_over_unit = len_over_unit)]
il <- fread(opt$ins_labels)[, .(ID, ins_best_label = best_label)]

ins <- Reduce(function(a, b) merge(a, b, by = "ID", all.x = TRUE), list(dec, info, mg, rot, il))
ins[is.na(rotated_tandem_dup), rotated_tandem_dup := FALSE]
ins[, ins_annotation := fcase(
  merged_label %like% "^MEI ", merged_label,
  merged_label == "tandem duplication (solved)", "tandem duplication (solved)",
  rotated_tandem_dup & abs(rotated_len_over_unit - 1) <= opt$rot_tol, "tandem duplication (rotated, site inside unit)",
  ins_best_label %in% c("STR", "VNTR", "satellite", "simple_lowcomp"), ins_best_label,
  ins_best_label %like% "^MEI_" | ins_best_label == "other_TE", paste("aligns to ref", ins_best_label),
  ins_best_label == "none", "no repeat",
  default = "no alignment")]

# per-hit flags against SVAN DUP_COORD segments and rotated units
parse_segs <- function(id, x) {
  s <- unlist(strsplit(x, ","))
  data.table(ID = id, chrom = sub(":.*", "", s), s = as.numeric(sub("^[^:]+:([0-9]+)-.*", "\\1", s)),
             e = as.numeric(sub("^[^:]+:[0-9]+-([0-9]+).*", "\\1", s)))
}
dup <- ins[!is.na(DUP_COORD), parse_segs(ID, DUP_COORD), by = seq_len(ins[!is.na(DUP_COORD), .N])][, seq_len := NULL]
unit <- ins[rotated_tandem_dup == TRUE, parse_segs(ID, rotated_unit), by = seq_len(ins[rotated_tandem_dup == TRUE, .N])][, seq_len := NULL]
h[, `:=`(hchrom = sub(":.*", "", coords), hs = as.numeric(sub("^[^:]+:([0-9]+)-.*", "\\1", coords)), he = as.numeric(sub(".*-", "", coords)))]
h[, row := .I]
ov <- function(segs, within) {
  x <- merge(h[, .(row, ID, hchrom, hs, he)], segs, by = "ID", allow.cartesian = TRUE)[hchrom == chrom]
  x <- if (within) x[hs >= s & he <= e] else x[hs <= e & he >= s]
  h$row %in% x$row
}
h[, hit_at_svan_dup := ov(dup, FALSE)]
h[, hit_in_rotated_unit := ov(unit, TRUE)]
h[, c("hchrom", "hs", "he", "row") := NULL]

out <- merge(h, ins, by = "ID", all.x = TRUE, sort = FALSE)[order(ID, -score)]
fwrite(out, opt$out, sep = "\t")

cat(sprintf("hits: %d; insertions with hits: %d\n", nrow(out), uniqueN(out$ID)))
cat("\n== insertions by ins_annotation (aligned insertions)\n")
print(dcast(unique(out[, .(ID, remap_class, ins_annotation)])[, .N, by = .(ins_annotation, remap_class)],
            ins_annotation ~ remap_class, value.var = "N", fill = 0))
cat("\n== hits by hit_label\n"); print(out[, .N, by = hit_label][order(-N)])
cat(sprintf("\nhits at a SVAN DUP_COORD segment: %d; hits inside a rotated tandem-dup unit: %d\n",
            out[hit_at_svan_dup == TRUE, .N], out[hit_in_rotated_unit == TRUE, .N]))
