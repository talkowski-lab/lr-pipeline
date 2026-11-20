#!/usr/bin/env Rscript
# Annotate the reference intervals that inserted sequences align to (truvari remap step-3 hits) with
# hg38 repeat tracks, and label each insertion by what its hits land on.
#
# Categories (per hit: fraction of the hit's reference interval covered, union per category):
#   from the RepeatMasker-style track (class / family columns):
#     MEI_Alu (SINE/Alu), MEI_L1 (LINE/L1), MEI_SVA (Retroposon/SVA), MEI_ERVK (LTR/ERVK),
#     satellite (Satellite/*), simple_lowcomp (Simple_repeat, Low_complexity), other_TE (any other class)
#   from the TR catalog (BED, 4th column with MOTIFS=...): STR (longest motif 1-6 bp), VNTR (longest motif >= 7 bp)
# A hit's label is the category covering the largest fraction of it if that fraction is >= --min-frac,
# preferring TR categories on ties; otherwise "none".
# An insertion is labelled from its best-scoring hit (alignment score), and the fraction of the insertion
# covered by hits of that label is reported.
# Needs bedtools on PATH.
suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
})

opt <- parse_args(OptionParser(
  usage = "%prog --hits <step3_hits.tsv.gz> --decisions <step3_decisions.tsv.gz> --repeats <REPEATS_hg38.bed.gz> --tr-catalog <TR.bed.gz> --out-prefix <prefix>",
  description = "Label remap hits and insertions as MEI / STR / VNTR / other repeat by hg38 tracks.",
  option_list = list(
    make_option("--hits", help = "Per-hit table from replay_remap_step3.R"),
    make_option("--decisions", help = "Per-insertion table from replay_remap_step3.R"),
    make_option("--repeats", help = "RepeatMasker-style BED (header line starting with #; columns ref beg end name class family repBeg repEnd strand milliDiv)"),
    make_option("--tr-catalog", dest = "tr", help = "TR catalog BED (chrom start end ID=...;MOTIFS=m1,m2;...), e.g. TRExplorer"),
    make_option("--info", help = "Optional INFO table from extract_ins_info.sh (TRID, SUB_FAMILY, allele_subtype) for a cross-tab against the VCF annotation"),
    make_option("--min-frac", dest = "min_frac", type = "double", default = 0.5, help = "Minimum fraction of a hit covered for a category label [default %default]"),
    make_option("--str-max-motif", dest = "str_max", type = "integer", default = 6, help = "Longest motif length counted as STR [default %default]"),
    make_option("--out-prefix", dest = "out_prefix", help = "Output prefix")
  )
))
for (a in c("hits", "decisions", "repeats", "tr", "out_prefix")) if (is.null(opt[[a]])) stop("--", a, " is required")
if (Sys.which("bedtools") == "") stop("bedtools not found on PATH")

tmp <- tempfile("hitrep_"); dir.create(tmp)

# ---- hits as BED ----
h <- fread(opt$hits)
h[, hit_id := .I]
h[, chrom := sub(":.*", "", coords)]
h[, start1 := as.numeric(sub("^[^:]+:([0-9]+)-.*", "\\1", coords))]
h[, end1 := as.numeric(sub(".*-", "", coords))]
hb <- h[, .(chrom, s = start1 - 1, e = end1, hit_id)][order(chrom, s)]
fwrite(hb, file.path(tmp, "hits.bed"), sep = "\t", col.names = FALSE, scipen = 50)

# ---- category BED from the two tracks ----
rep <- fread(cmd = paste("gzip -dc", shQuote(opt$repeats), "| grep -v '^#'"), header = FALSE, sep = "\t",
             select = c(1:3, 5, 6), col.names = c("chrom", "s", "e", "class", "family"))
rep[, cat := fcase(class == "SINE" & family == "Alu", "MEI_Alu",
                   class == "LINE" & family == "L1", "MEI_L1",
                   class == "Retroposon" & family == "SVA", "MEI_SVA",
                   class == "LTR" & family == "ERVK", "MEI_ERVK",
                   class == "Satellite", "satellite",
                   class %in% c("Simple_repeat", "Low_complexity"), "simple_lowcomp",
                   default = "other_TE")]
tr <- fread(cmd = paste("gzip -dc", shQuote(opt$tr)), header = FALSE, sep = "\t", select = 1:4, col.names = c("chrom", "s", "e", "info"))
tr[, motifs := sub(";.*", "", sub(".*MOTIFS=", "", info))]
tr[, max_motif := vapply(strsplit(motifs, ","), function(m) max(nchar(m)), integer(1))]
tr[, cat := fifelse(max_motif <= opt$str_max, "STR", "VNTR")]
cats <- rbind(rep[, .(chrom, s, e, cat)], tr[, .(chrom, s, e, cat)])
cats <- cats[chrom %in% unique(hb$chrom)][order(chrom, s)]
fwrite(cats, file.path(tmp, "cats.bed"), sep = "\t", col.names = FALSE, scipen = 50)
rm(rep, tr)

# ---- per hit, per category: union bp covered ----
# merge each category separately, then intersect
cat_levels <- c("STR", "VNTR", "MEI_Alu", "MEI_L1", "MEI_SVA", "MEI_ERVK", "satellite", "simple_lowcomp", "other_TE")
cov <- rbindlist(lapply(cat_levels, function(k) {
  f <- file.path(tmp, paste0(k, ".bed"))
  system(sprintf("awk -F'\\t' '$4==\"%s\"' %s | bedtools merge -i - > %s", k, shQuote(file.path(tmp, "cats.bed")), shQuote(f)))
  if (file.info(f)$size == 0) return(NULL)
  x <- fread(cmd = sprintf("bedtools intersect -a %s -b %s -wo", shQuote(file.path(tmp, "hits.bed")), shQuote(f)),
             header = FALSE, sep = "\t", col.names = c("chrom", "s", "e", "hit_id", "c2", "s2", "e2", "ov"))
  if (nrow(x) == 0) return(NULL)
  x[, .(bp = sum(ov)), by = hit_id][, cat := k]
}))
cov <- merge(cov, h[, .(hit_id, hit_len = end1 - start1 + 1)], by = "hit_id")
cov[, frac := pmin(1, bp / hit_len)]
wide <- dcast(cov, hit_id ~ cat, value.var = "frac", fill = 0)
for (k in setdiff(cat_levels, names(wide))) wide[, (k) := 0]
setnames(wide, cat_levels, paste0("f_", cat_levels))
h <- merge(h, wide, by = "hit_id", all.x = TRUE)
for (k in paste0("f_", cat_levels)) set(h, which(is.na(h[[k]])), k, 0)

fcols <- paste0("f_", cat_levels)
h[, top_frac := do.call(pmax, .SD), .SDcols = fcols]
h[, hit_label := {
  m <- as.matrix(.SD); lab <- cat_levels[max.col(m, ties.method = "first")]   # ties -> TR first (column order)
  fifelse(top_frac >= opt$min_frac, lab, "none")
}, .SDcols = fcols]
h[, hit_group := fcase(hit_label %in% c("STR", "VNTR"), hit_label,
                       hit_label %like% "^MEI_", "MEI",
                       hit_label == "none", "none",
                       default = "other_repeat")]
fwrite(h[, c("ID", "coords", "strand", "q_start", "q_end", "perc", "identity", "mapq", "score", "type", "full", fcols, "hit_label", "hit_group"), with = FALSE],
       paste0(opt$out_prefix, ".hit_repeats.tsv.gz"), sep = "\t")

# ---- per insertion: label from the best-scoring hit; query coverage by hits with the same label ----
union_len <- function(s, e) {
  o <- order(s); s <- s[o]; e <- e[o]; tot <- 0; cs <- s[1]; ce <- e[1]
  for (i in seq_along(s)[-1]) { if (s[i] > ce + 1) { tot <- tot + ce - cs + 1; cs <- s[i]; ce <- e[i] } else ce <- max(ce, e[i]) }
  tot + ce - cs + 1
}
best <- h[order(ID, -score, -perc)][, .SD[1], by = ID][, .(ID, best_hit = coords, best_label = hit_label, best_group = hit_group,
                                                             best_top_frac = round(top_frac, 3), best_perc = round(perc, 2),
                                                             best_identity = identity, best_mapq = mapq)]
qcov <- merge(h, best[, .(ID, best_label)], by = "ID")[hit_label == best_label & best_label != "none",
            .(label_query_cov = round(union_len(q_start, q_end) / qlen[1], 3), n_label_hits = .N), by = ID]
grp <- h[, .(n_hits = .N, n_MEI_hits = sum(hit_group == "MEI"), n_TR_hits = sum(hit_group %in% c("STR", "VNTR")),
             labels_seen = paste(sort(unique(hit_label)), collapse = ",")), by = ID]

dec <- fread(opt$decisions)[, .(ID, alt_len, remap_class)]
out <- merge(dec, merge(merge(best, qcov, by = "ID", all.x = TRUE), grp, by = "ID"), by = "ID", all.x = TRUE)
out[is.na(best_label), `:=`(best_label = "no alignment", best_group = "no alignment")]
fwrite(out, paste0(opt$out_prefix, ".insertion_repeat_labels.tsv.gz"), sep = "\t")

o <- function(x, name) { fwrite(x, sprintf("%s.%s.tsv", opt$out_prefix, name), sep = "\t"); cat("\n==", name, "\n"); print(x) }
o(h[, .N, by = hit_label][order(-N)], "hit_label_counts")
o(dcast(out[, .N, by = .(remap_class, best_label)], best_label ~ remap_class, value.var = "N", fill = 0), "insertion_label_by_remap_class")
if (!is.null(opt$info)) {
  info <- fread(opt$info, sep = "\t", quote = "", na.strings = ".")[, .(ID, TRID, SUB_FAMILY, allele_subtype)]
  v <- merge(out[best_label != "no alignment"], info, by = "ID")
  v[, vcf := fcase(grepl("^Alu", SUB_FAMILY), "MEI Alu", grepl("^L1", SUB_FAMILY), "MEI L1", grepl("^SVA", SUB_FAMILY), "MEI SVA",
                   !is.na(allele_subtype), allele_subtype, !is.na(TRID), "in TR (TRID)", default = "plain ins")]
  o(dcast(v[, .N, by = .(best_label, vcf)], best_label ~ vcf, value.var = "N", fill = 0), "insertion_label_by_vcf")
}
unlink(tmp, recursive = TRUE)
