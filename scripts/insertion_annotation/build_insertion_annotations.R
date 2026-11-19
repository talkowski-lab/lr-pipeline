#!/usr/bin/env Rscript
# One row per insertion (every insertion in the remap TSV, aligned or not) with its final call and,
# for complex insertions, where each piece sits on the inserted sequence and where it aligns in hg38.
#
# Complex insertion = multi_piece from tile_insertion_pieces.R (pieces from >= 2 loci, little overlap,
# tiles cover most of the insertion, no single piece explains it).
# ins_call, first match wins:
#   VCF/override MEI (VCF SUB_FAMILY for every insertion, aligned or not) > composite satellite (optional --composite, unit span covers >= --composite-min-cov)
#   > solved tandem dup > complex insertion > rotated tandem dup > copy of a reference TE (--ref-te-min of the best
#   alignment in one TE family and of the insertion; failing insertions get the best alignment's non-TE label
#   or "no repeat") > ins_annotation from
#   annotate_step3_hits.R (STR / VNTR / satellite / ... / no repeat) > "no alignment".
# Piece columns (complex insertions only, pieces in order along the inserted sequence):
#   PIECE_POS: start-end of each piece on the inserted sequence, 1-based, anchor base removed
#              (ALT minus its first base), comma-separated
#   PIECE_REF: hg38 alignment of each piece as chrom:start-end_strand (1-based, inclusive; strand of the piece
#              relative to hg38 +), comma-separated, same order as PIECE_POS
#   PIECE_IDENTITY, PIECE_MAPQ: per piece, same order
# Composite satellite columns (call_composite_repeat_copies.R; filled for every insertion with >= 1 anchor copy):
#   COMPOSITE_REPEAT, COMPOSITE_ORIENTATION, COMPOSITE_COPIES, COMPOSITE_FULL_COPIES (all three elements present),
#   COMPOSITE_COV (fraction of the insertion spanned by the copies), COMPOSITE_COPY_LEN (bp per copy, in order),
#   COMPOSITE_COPY_LAYOUT (element lengths per copy, "<before>:bp|<anchor>:bp|<after>:bp", copies separated by ";")
suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
})

opt <- parse_args(OptionParser(
  usage = "%prog --decisions <step3_decisions.tsv.gz> --hits <step3_hits.annotated.tsv.gz> --tiled <tiled_pieces.tsv> --pieces <multi_piece.pieces.tsv> --info <info.tsv.gz> --out <insertion_annotations.tsv.gz>",
  description = "Per-insertion annotation table with complex-insertion piece positions and hg38 alignments.",
  option_list = list(
    make_option("--decisions", help = "Per-insertion table from replay_remap_step3.R (all insertions)"),
    make_option("--hits", help = "Annotated hit table from annotate_step3_hits.R"),
    make_option("--tiled", help = "Per-insertion table from tile_insertion_pieces.R"),
    make_option("--pieces", help = "Per-piece table from tile_insertion_pieces.R"),
    make_option("--info", help = "INFO table from extract_ins_info.sh"),
    make_option("--composite", default = NULL, help = "Optional per-insertion table from call_composite_repeat_copies.R"),
    make_option("--composite-min-cov", dest = "comp_cov", type = "double", default = 0.8,
                help = "Min COMPOSITE_COV to call the insertion a composite satellite [default %default]"),
    make_option("--ref-te-min", dest = "ref_te_min", type = "double", default = 0.9,
                help = "Copy of a reference transposon: min fraction of the best alignment in one TE family AND min fraction of the insertion that is that TE [default %default]"),
    make_option("--keep-anchor", dest = "keep_anchor", type = "logical", default = TRUE,
                help = "Query sequences in the alignments include the ALT anchor base [default %default]"),
    make_option("--out", help = "Output TSV (.gz)")
  )
))
for (a in c("decisions", "hits", "tiled", "pieces", "info", "out")) if (is.null(opt[[a]])) stop("--", a, " is required")

dec <- fread(opt$decisions)[, .(ID, CHROM, POS, alt_len, remap_class, remap_coords, remap_ori, remap_perc, n_hits, n_full)]
info <- fread(opt$info, sep = "\t", quote = "", na.strings = ".")[, .(ID, FILTER, AF, AC, REGION, TRID, TR_ENVELOPED, SUB_FAMILY, allele_subtype)]
h <- fread(opt$hits, na.strings = c("", "NA"))
ins_cols <- c("ID", "svan_type", "svan_mei_fam", "DUP_COORD", "tandem_dup", "merged_label",
              "rotated_tandem_dup", "rotated_unit", "rotated_len_over_unit", "ins_best_label", "ins_annotation")
ia <- unique(h[, ..ins_cols])
# reference transposon copy from the best-scoring alignment: one TE family covers >= ref_te_min of the alignment's
# reference interval and >= ref_te_min of the insertion (family fraction x alignment coverage)
te_cols <- c(MEI_Alu = "f_MEI_Alu", MEI_L1 = "f_MEI_L1", MEI_SVA = "f_MEI_SVA", MEI_ERVK = "f_MEI_ERVK", other_TE = "f_other_TE")
nte_cols <- c(STR = "f_STR", VNTR = "f_VNTR", satellite = "f_satellite", simple_lowcomp = "f_simple_lowcomp")
bh <- h[order(ID, -score, -perc)][, .SD[1], by = ID]
tm <- as.matrix(bh[, ..te_cols]); nm <- as.matrix(bh[, ..nte_cols])
bh[, `:=`(ref_te_family = names(te_cols)[max.col(tm, ties.method = "first")], ref_te_frac = apply(tm, 1, max))]
bh[, ref_te_frac_of_ins := ref_te_frac * pmin(perc, 100) / 100]
bh[, nonte_label := fifelse(apply(nm, 1, max) >= 0.5, names(nte_cols)[max.col(nm, ties.method = "first")], "no repeat")]
ia <- merge(ia, bh[, .(ID, ref_te_family, ref_te_frac = round(ref_te_frac, 4), ref_te_frac_of_ins = round(ref_te_frac_of_ins, 4), nonte_label)],
            by = "ID", all.x = TRUE)
ia[, ref_te_copy := !is.na(ref_te_frac) & ref_te_frac >= opt$ref_te_min & ref_te_frac_of_ins >= opt$ref_te_min]
til <- fread(opt$tiled, na.strings = c("", "NA"))[, .(ID, complex_insertion = multi_piece, n_pieces, n_loci, tile_cov, max_piece_frac)]

off <- if (opt$keep_anchor) 1L else 0L
pc <- fread(opt$pieces)[order(ID, q_start)]
pc[, `:=`(ps = pmax(1L, q_start - off), pe = q_end - off)]
pc[, `:=`(rchrom = sub(":.*", "", coords), rs = sub("^[^:]+:([0-9]+)-.*", "\\1", coords), re = sub(".*-", "", coords))]
pcs <- pc[, .(PIECE_POS = paste(sprintf("%d-%d", ps, pe), collapse = ","),
              PIECE_REF = paste(sprintf("%s:%s-%s_%s", rchrom, rs, re, strand), collapse = ","),
              PIECE_IDENTITY = paste(sprintf("%.4f", identity), collapse = ","),
              PIECE_MAPQ = paste(mapq, collapse = ",")), by = ID]

x <- Reduce(function(a, b) merge(a, b, by = "ID", all.x = TRUE), list(dec, info, ia, til, pcs))
if (!is.null(opt$composite)) {
  cm <- fread(opt$composite)[, .(ID, COMPOSITE_REPEAT = composite_repeat, COMPOSITE_ORIENTATION = composite_orientation,
                                 COMPOSITE_COPIES = composite_copies, COMPOSITE_FULL_COPIES = composite_full_copies,
                                 COMPOSITE_COV = composite_cov, COMPOSITE_COPY_LEN = composite_copy_len,
                                 COMPOSITE_COPY_LAYOUT = composite_copy_layout)]
  x <- merge(x, cm, by = "ID", all.x = TRUE)
} else {
  x[, `:=`(COMPOSITE_REPEAT = NA_character_, COMPOSITE_ORIENTATION = NA_character_, COMPOSITE_COPIES = NA_integer_,
           COMPOSITE_FULL_COPIES = NA_integer_, COMPOSITE_COV = NA_real_, COMPOSITE_COPY_LEN = NA_character_,
           COMPOSITE_COPY_LAYOUT = NA_character_)]
}
x[is.na(complex_insertion), complex_insertion := FALSE]
x[is.na(rotated_tandem_dup), rotated_tandem_dup := FALSE]
x[complex_insertion == FALSE, c("PIECE_POS", "PIECE_REF", "PIECE_IDENTITY", "PIECE_MAPQ") := NA_character_]
# MEI family from the VCF SUB_FAMILY for every insertion, aligned or not (merged_label, which carries manual
# overrides, exists only for insertions with alignments and takes precedence)
x[, vcf_mei := fcase(grepl("^Alu", SUB_FAMILY), "MEI Alu", grepl("^L1", SUB_FAMILY), "MEI L1", grepl("^SVA", SUB_FAMILY), "MEI SVA",
                     default = NA_character_)]
x[, ins_call := fcase(
  merged_label %like% "^MEI ", merged_label,
  !is.na(vcf_mei), vcf_mei,
  !is.na(COMPOSITE_COV) & COMPOSITE_COV >= opt$comp_cov, paste("satellite", COMPOSITE_REPEAT),
  merged_label %in% "tandem duplication (solved)", "tandem duplication (solved)",
  complex_insertion == TRUE, "complex insertion",
  ins_annotation %like% "^tandem duplication \\(rotated", ins_annotation,
  ref_te_copy %in% TRUE, paste("copy of ref", ref_te_family),
  ins_annotation %like% "^aligns to ref ", nonte_label,
  !is.na(ins_annotation), ins_annotation,
  default = "no alignment")]
setcolorder(x, intersect(c("ID", "CHROM", "POS", "alt_len", "remap_class", "ins_call", "complex_insertion", "n_pieces", "n_loci",
                           "PIECE_POS", "PIECE_REF", "PIECE_IDENTITY", "PIECE_MAPQ", "COMPOSITE_REPEAT", "COMPOSITE_ORIENTATION",
                           "COMPOSITE_COPIES", "COMPOSITE_FULL_COPIES", "COMPOSITE_COV", "COMPOSITE_COPY_LEN", "COMPOSITE_COPY_LAYOUT"),
                         names(x)))
setorder(x, CHROM, POS, ID)
fwrite(x, opt$out, sep = "\t")

cat(sprintf("insertions: %d; complex insertions: %d (called complex: %d)\n", nrow(x), x[complex_insertion == TRUE, .N], x[ins_call == "complex insertion", .N]))
cat("\n== ins_call x remap class\n")
print(dcast(x[, .N, by = .(ins_call, remap_class)], ins_call ~ remap_class, value.var = "N", fill = 0))
cat("\n== composite satellite calls\n")
print(x[ins_call %like% "^satellite ", .(ID, remap_class, COMPOSITE_COPIES, COMPOSITE_COV, COMPOSITE_COPY_LEN)])
cat("\n== complex insertions\n")
print(x[complex_insertion == TRUE, .(ID, ins_call, PIECE_POS, PIECE_REF, PIECE_MAPQ)])
