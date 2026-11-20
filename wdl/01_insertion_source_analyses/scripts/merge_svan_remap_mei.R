#!/usr/bin/env Rscript
# Merge SVAN, truvari remap and the final gnomAD-LR VCF MEI annotation into one label per insertion.
#  1. Tandem duplication "solved": SVAN ITYPE_N == DUP, remap class == tandem, a SVAN DUP_COORD
#     segment at the remap hit locus, and that segment's length within --len-tol of the insertion length.
#  2. MEI: the family comes from the VCF SUB_FAMILY (final callset annotation), with manual overrides;
#     SVAN FAM_N is compared against it, and the remap class is reported for each MEI.
suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
})

opt <- parse_args(OptionParser(
  usage = "%prog --joined <svan_remap.joined.tsv.gz> --out-prefix <prefix> [--overrides <tsv>] [--len-tol 0.1]",
  description = "Label insertions as solved tandem duplications or MEIs from SVAN + remap + VCF SUB_FAMILY.",
  option_list = list(
    make_option("--joined", help = "Joined table from compare_svan_remap.R (carries VCF SUB_FAMILY)"),
    make_option("--overrides", help = "Optional TSV with header ID, family, subfamily, note: manual MEI family calls that replace the VCF annotation"),
    make_option("--len-tol", dest = "len_tol", type = "double", default = 0.1,
                help = "Max |source length / insertion length - 1| for a solved tandem dup [default %default]"),
    make_option("--out-prefix", dest = "out_prefix", help = "Output prefix")
  )
))
for (a in c("joined", "out_prefix")) if (is.null(opt[[a]])) stop("--", a, " is required")

d <- fread(opt$joined, na.strings = c("", "NA", "."))
for (col in names(d)[sapply(d, is.character)]) set(d, which(d[[col]] %in% c("", ".")), col, NA_character_)

# ---- 1. tandem duplications ----
seg_start <- function(x) as.numeric(sub("-.*", "", sub("^.*:", "", x)))
seg_end <- function(x) as.numeric(sub("_.*", "", sub("^[^-]*-", "", x)))
hit_start <- seg_start(d$coords)
d[, hit_len := seg_end(coords) - seg_start(coords) + 1]
# SVAN DUP_COORD segment matching the remap hit locus (same contig, start within max(len, 1 kb)); NA if none
d[, svan_src_len := mapply(function(s, hc, hs, l) {
  if (is.na(s) || is.na(hc)) return(NA_real_)
  s <- strsplit(s, ",")[[1]]
  ok <- sub(":.*", "", s) == hc & abs(seg_start(s) - hs) <= max(l, 1000)
  if (!any(ok)) return(NA_real_)
  len_s <- seg_end(s[ok]) - seg_start(s[ok]) + 1
  len_s[which.min(abs(len_s - l))]
}, DUP_COORD, hit_chr, hit_start, len)]
d[, src_len_ratio := round(svan_src_len / len, 3)]

d[, tandem_dup := fcase(
  svan_type != "DUP", NA_character_,
  class == "tandem" & !is.na(svan_src_len) & abs(src_len_ratio - 1) <= opt$len_tol, "solved",
  class == "tandem" & !is.na(svan_src_len), "type+source agree, length differs",
  class == "tandem", "type agrees, source differs",
  class == "partial" & !is.na(svan_src_len), "remap partial at SVAN source",
  class == "novel", "remap no hit",
  default = "remap other class")]

# ---- 2. MEIs ----
norm_fam <- function(x) fcase(is.na(x), NA_character_,
                              grepl("^ALU|^Alu", x), "Alu",
                              grepl("^LINE|^L1", x), "L1",
                              grepl("^SVA", x), "SVA",
                              default = x)
d[, vcf_subfamily := SUB_FAMILY]
d[, vcf_fam := norm_fam(SUB_FAMILY)]
d[, mei_override := NA_character_]
if (!is.null(opt$overrides)) {
  ov <- fread(opt$overrides, sep = "\t", quote = "")
  stopifnot(all(c("ID", "family", "subfamily", "note") %in% names(ov)), !anyDuplicated(ov$ID))
  cat(sprintf("overrides: %d; in this insertion set: %d\n", nrow(ov), sum(ov$ID %in% d$ID)))
  d[ov, on = "ID", `:=`(vcf_fam = i.family, vcf_subfamily = i.subfamily, mei_override = i.note)]
}
d[, svan_fam := fifelse(svan_mei_fam == "none", NA_character_, norm_fam(svan_mei_fam))]
d[, mei_call := fcase(
  !is.na(vcf_fam) & !is.na(svan_fam) & vcf_fam == svan_fam, "VCF + SVAN agree",
  !is.na(vcf_fam) & !is.na(svan_fam), "VCF + SVAN family differ",
  !is.na(vcf_fam), "VCF only",
  !is.na(svan_fam), "SVAN only",
  default = NA_character_)]
d[, mei_family := fcoalesce(vcf_fam, svan_fam)]

# ---- merged label ----
# The VCF MEI annotation is trusted and takes precedence over the remap: a VCF MEI is called MEI even
# when remap says tandem (or SVAN says DUP). SVAN-only MEI calls are reported, not adopted.
d[, merged_label := fcase(
  !is.na(vcf_fam), paste0("MEI ", vcf_fam),
  tandem_dup %in% "solved", "tandem duplication (solved)",
  mei_call %in% "SVAN only", paste0("SVAN-only MEI call (", svan_fam, "), not in VCF"),
  default = "unresolved")]

fwrite(d, paste0(opt$out_prefix, ".merged_labels.tsv.gz"), sep = "\t")

out <- function(x, name) {
  fwrite(x, sprintf("%s.%s.tsv", opt$out_prefix, name), sep = "\t")
  cat("\n==", name, "\n"); print(x)
}
out(d[!is.na(tandem_dup), .N, by = tandem_dup][order(-N)], "tandem_dup_status")
out(d[!is.na(src_len_ratio) & svan_type == "DUP" & class == "tandem",
      .(n = .N, ratio_min = min(src_len_ratio), ratio_q25 = quantile(src_len_ratio, .25),
        ratio_median = median(src_len_ratio), ratio_q75 = quantile(src_len_ratio, .75), ratio_max = max(src_len_ratio))],
    "tandem_dup_source_length_ratio")
out(dcast(d[!is.na(mei_call), .N, by = .(mei_call, mei_family)], mei_call ~ mei_family, value.var = "N", fill = 0), "mei_concordance")
out(dcast(d[!is.na(mei_call), .N, by = .(mei_call, class)], mei_call ~ class, value.var = "N", fill = 0), "mei_by_remap_class")
out(d[mei_call %in% "VCF + SVAN family differ", .(ID, len, vcf_subfamily, svan_mei_fam, svan_type, class)], "mei_family_conflicts")
out(d[mei_call %in% "SVAN only", .(ID, len, svan_mei_fam, svan_type, class, coords, ori, perc)], "svan_only_mei_calls")
out(d[mei_call %in% "SVAN only" & class == "tandem", .(ID, coords, ori, perc, len, svan_mei_fam, svan_type, DUP_COORD, REGION, AF)],
    "svan_only_mei_tandem_hits")
out(d[!is.na(vcf_fam) & (class %in% c("tandem", "tandem_inverted", "tandem_complex") | svan_type %in% "DUP"),
      .(ID, len, vcf_subfamily, svan_type, class, coords, tandem_dup)], "vcf_mei_vs_tandem")
out(d[, .N, by = merged_label][order(-N)], "merged_label_counts")
