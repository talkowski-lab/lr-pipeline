#!/usr/bin/env Rscript
# One final annotation per insertion, combining the stages (first match wins):
#   repeat class (build_repeat_classes.py), unless "non-repeat (next stage)"
#   duplication call (combine_dup_stage.R): tandem / inverted tandem / interspersed (where the caller placed the
#   insertion point inside the unit is kept in the detail)
#   CNV call (annotate_cnv.R): CNV_<copies>copy[_<locus>]
#   CPX_INS / CPX_SV (--cpx, optional): CPX_INS = pieces from several regions in insertion order with orientation
#   (make_cpx_annotation.py); CPX_SV = curated complex SV structure from the recurrence plots
#   T2T label (name_t2t_matches.py): T2T_<genes or location> = sequence present in T2T-CHM13
#   local interspersed duplication (--local-dups, optional; align_to_local_origin.py): >= 90% of the insertion aligned
#   to a nearby origin that does not contain the insertion site
#   BWA call (--bwa, optional; summarize_bwa_alignments.py): >= 90% of the insertion in one BWA alignment; NUMT when on
#   chrM, else tandem duplication (copy within 5 bp of the site; "insertion point inside the unit" in the detail when
#   it spans the site) or interspersed duplication; detail = origin interval, strand, identity, mapq, distance, origin cytoband
#   iterative-alignment call (--iterative, optional; classify_iterative_layouts.py): "rotated tandem duplication" ->
#   tandem duplication (insertion point inside the unit in the detail); "tandem duplication (pieces)" -> tandem duplication
#   diverged copy (--iterative low-identity matches; with --site-context the class is "repeat region, diverged copy" when
#   the site has a context (SD block / SD / subtelomeric / pericentromeric / acrocentric arm / TR / simple repeat),
#   which goes in the detail with the weighted identity and pieces)
#   SD link (--sd-links, optional): when the insertion site and its origin lie in the two copies of one SD pair, the link
#   family ("<origin copy> inserted into <site copy>", family, pair identity) is appended to the detail of any class
#   last step (--last-step, optional; assemble_local_tandem.py): tandem duplication, CPX_SV (local pieces out of tandem
#   order) or CPX_INS (chimeric: pieces from several places, origins listed)
#   unresolved (with the immune-locus call from combine_next_stage.R when there is one)
suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
})
opt <- parse_args(OptionParser(
  usage = "%prog --classes <repeat_classes.tsv> --dups <dup_stage.tsv> --cnv <cnv_calls.tsv> --next <next_stage_calls.tsv> --t2t <t2t_labels.tsv> --out-prefix <prefix>",
  description = "Combine all stages into one annotation per insertion.",
  option_list = list(
    make_option("--classes", help = "build_repeat_classes.py output"),
    make_option("--dups", help = "combine_dup_stage.R <prefix>.dup_stage.tsv"),
    make_option("--cnv", help = "annotate_cnv.R output"),
    make_option("--next", dest = "nxt", help = "combine_next_stage.R <prefix>.next_stage_calls.tsv"),
    make_option("--t2t", help = "name_t2t_matches.py output"),
    make_option("--cpx", default = NULL,
                help = "Optional comma-separated tables (ID, cpx_call, cpx_detail): make_cpx_annotation.py CPX_INS and curated CPX_SV"),
    make_option("--local-dups", dest = "local_dups", default = NULL,
                help = "Optional align_to_local_origin.py output; rows called 'interspersed duplication' are used"),
    make_option("--bwa", default = NULL,
                help = "Optional summarize_bwa_alignments.py output (full-length calls: NUMT, tandem / interspersed duplication)"),
    make_option("--iterative", default = NULL,
                help = "Optional classify_iterative_layouts.py output (rotated / in-pieces tandem duplications are used)"),
    make_option("--site-context", dest = "site_context", default = NULL,
                help = "Optional site_context.py output (ID, context) for the low-identity matches"),
    make_option("--sd-links", dest = "sd_links", default = NULL,
                help = "Optional link_sd_insertions.py output; between-copy SD links are appended to the detail"),
    make_option("--last-step", dest = "last_step", default = NULL,
                help = "Optional assemble_local_tandem.py output: tandem duplication / complex SV (CPX_SV) / chimeric (CPX_INS)"),
    make_option("--out-prefix", dest = "out_prefix", help = "Output prefix")
  )
))
for (a in c("classes", "dups", "cnv", "nxt", "t2t", "out_prefix")) if (is.null(opt[[a]])) stop("--", a, " is required")
na <- c("", "NA")
# read a stage table; an empty table (header only) gets typed columns (numeric for counts / fractions, character for the
# rest) so merges and fcase still work
rd <- function(f) {
  d <- fread(f, na.strings = na)
  if (nrow(d) == 0) {
    for (k in names(d)) set(d, j = k, value = if (grepl("perc|cov|identity|len|dist|mapq|frac|n_", k)) numeric() else character())
  }
  if ("ID" %in% names(d)) d[, ID := as.character(ID)]
  d
}
x <- rd(opt$classes)[, .(ID, len, repeat_class, repeat_detail)]
d <- rd(opt$dups)[, .(ID, dup_call, dup_detail)]
cn <- rd(opt$cnv)[, .(ID, cnv_call, cnv_detail)]
nx <- rd(opt$nxt)[, .(ID, next_call)]
t2 <- rd(opt$t2t)[, .(ID, t2t_label, T2T_best_hit, T2T_best_perc)]
cp <- if (!is.null(opt$cpx)) rbindlist(lapply(strsplit(opt$cpx, ",")[[1]], function(f) rd(f)[, .(ID, cpx_call, cpx_detail)])) else data.table(ID = character(), cpx_call = character(), cpx_detail = character())
ld <- data.table(ID = character(), local_call = character(), local_detail = character())
if (!is.null(opt$local_dups)) {
  ld <- rd(opt$local_dups)[call == "interspersed duplication"]
  ld <- ld[, .(ID, local_call = call,
               local_detail = sprintf("origin %s (%s), %.0f%% of the insertion aligned at %.0f%% identity, %d bp from the site",
                                      origin_aligned, strand, 100 * aligned_frac, 100 * identity, dist_to_site))]
}
nm <- data.table(ID = character(), bwa_call = character(), bwa_detail = character())
if (!is.null(opt$bwa)) {
  nm <- rd(opt$bwa)[!call %like% "^partial|^unmapped"]
  nm <- nm[, .(ID, bwa_call = fifelse(aligned_location %like% "^chrM:", "NUMT", fifelse(call %like% "^interspersed", "interspersed duplication", call)),
               bwa_detail = sprintf("origin %s (%s)%s, %.0f%% of the insertion at %.1f%% identity, mapq %d, %s%s; BWA",
                                    aligned_location, strand, fifelse(is.na(origin_band) | origin_band == "", "", paste0(" [", origin_band, "]")),
                                    100 * best_qcov, 100 * best_identity, best_mapq,
                                    fifelse(is.na(dist_to_site), "another contig", paste0(trimws(format(dist_to_site, big.mark = ",")), " bp from the site")),
                                    fifelse(relation == "spans site", ", insertion point inside the unit", "")))]
}
it <- data.table(ID = character(), it_call = character(), it_detail = character())
if (!is.null(opt$iterative)) {
  itall <- rd(opt$iterative)
  sc <- if (!is.null(opt$site_context)) rd(opt$site_context) else data.table(ID = character(), context = character())
  lo <- merge(itall[class == "low-identity match"], sc, by = "ID", all.x = TRUE)
  lo <- lo[, .(ID, it_call = fifelse(!is.na(context) & context != "", "repeat region, diverged copy", "diverged copy"),
               it_detail = sprintf("%s%.1f%% weighted identity over %d pieces on %s: %s; iterative alignment",
                                   fifelse(!is.na(context) & context != "", paste0(context, "; "), ""), 100 * weighted_identity,
                                   n_segments, fifelse(ref_used == "site", "hg38", ref_used), layout))]
  it <- itall[class %in% c("rotated tandem duplication", "tandem duplication (pieces)")]
  it <- it[, .(ID, it_call = "tandem duplication",
               it_detail = sprintf("%s%s; %d pieces at %.1f%% weighted identity: %s; iterative alignment",
                                   fifelse(class == "rotated tandem duplication", "insertion point inside the unit, ", ""), note, n_segments,
                                   100 * weighted_identity, layout))]
  it <- rbind(it, lo)
}
ls_ <- data.table(ID = character(), ls_call = character(), ls_detail = character())
if (!is.null(opt$last_step)) {
  ls_ <- rd(opt$last_step)[class != "unresolved"]
  ls_ <- ls_[, .(ID, ls_call = fcase(class == "tandem duplication", "tandem duplication", class == "complex SV", "CPX_SV",
                                     default = "CPX_INS"),
                 ls_detail = fcase(class == "tandem duplication",
                                   sprintf("unit %s; %s; local %.0f%% at %.1f%%; last step", unit, note, 100 * local_cov,
                                           100 * local_identity),
                                   class == "complex SV",
                                   sprintf("local pieces out of tandem order: %s; %.0f%% at %.1f%%; last step", local_layout,
                                           100 * local_cov, 100 * local_identity),
                                   default = sprintf("pieces: %s; %.0f%% at %.1f%%; last step", genome_layout, 100 * genome_cov,
                                                     100 * genome_identity)))]
}
x <- Reduce(function(a, b) merge(a, b, by = "ID", all.x = TRUE), list(x, d, cn, nx, t2, cp, ld, nm, it, ls_))
x[, final_class := fcase(
  repeat_class != "non-repeat (next stage)", repeat_class,
  !is.na(dup_call) & dup_call != "next stage", dup_call,
  !is.na(cnv_call) & cnv_call %like% "^CNV_", cnv_call,
  !is.na(cpx_call), cpx_call,
  !is.na(t2t_label), t2t_label,
  !is.na(local_call), local_call,
  !is.na(bwa_call), bwa_call,
  !is.na(it_call), it_call,
  !is.na(ls_call), ls_call,
  !is.na(next_call) & next_call %like% "IGL", paste0("unresolved; ", next_call),
  default = "unresolved")]
x[, final_detail := fcase(
  repeat_class != "non-repeat (next stage)", repeat_detail,
  !is.na(dup_call) & dup_call != "next stage", dup_detail,
  !is.na(cnv_call) & cnv_call %like% "^CNV_", cnv_detail,
  !is.na(cpx_call), cpx_detail,
  !is.na(t2t_label), sprintf("present in T2T-CHM13: %s, %.0f%% of the insertion", T2T_best_hit, T2T_best_perc),
  !is.na(local_call), local_detail,
  !is.na(bwa_call), bwa_detail,
  !is.na(it_call), it_detail,
  !is.na(ls_call), ls_detail,
  default = NA_character_)]
x[, final_group := fcase(final_class %like% "^MEI |^VNTR in MEI|insertion into|insertion, no local|^MEI nested", "MEI",
                         final_class %in% c("STR", "VNTR", "STR_novel", "VNTR_novel", "other TR"), "STR / VNTR",
                         final_class %like% "^satellite", "satellite",
                         final_class %like% "duplication|paralog", "duplication",
                         final_class %like% "^CNV_", "CNV",
                         final_class %in% c("CPX_INS", "CPX_SV"), "complex insertion / SV",
                         final_class %like% "^T2T_", "present in T2T",
                         final_class == "NUMT", "NUMT",
                         final_class %like% "diverged copy", "diverged copy",
                         default = "unresolved")]
if (!is.null(opt$sd_links)) {
  sl <- rd(opt$sd_links)[best_link %like% "inserted into",
                                              .(ID, sd_link = sprintf("SD link: %s (family %s, %s)", best_link, family,
                                                                      fifelse(is.na(best_identity), "same family",
                                                                              sprintf("%.1f%%", 100 * best_identity))))]
  x <- merge(x, sl, by = "ID", all.x = TRUE)
  x[!is.na(sd_link), final_detail := fifelse(is.na(final_detail) | final_detail == "", sd_link, paste0(final_detail, "; ", sd_link))]
}
if (anyDuplicated(x$ID)) stop("duplicate IDs after merging the stage tables (an input lists an ID twice): ",
                               paste(head(unique(x$ID[duplicated(x$ID)]), 5), collapse = ", "))
setorder(x, ID)
fwrite(x[, .(ID, len, final_group, final_class, final_detail)], paste0(opt$out_prefix, ".final_annotation.tsv"), sep = "\t")
s <- x[, .N, by = final_group][order(-N)]
fwrite(s, paste0(opt$out_prefix, ".final_annotation_counts.tsv"), sep = "\t")
print(s)
