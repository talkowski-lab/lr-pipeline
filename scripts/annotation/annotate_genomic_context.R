#!/usr/bin/env Rscript

# Assign genomic context labels to variants from breakpoint and CNV body coverage tables
library("optparse")

option_list = list(
    make_option(c('-i', "--input"), type="character", default=NULL, help="input sites file"),
    make_option(c('-o', "--output"), type="character", default=NULL, help="output TSV file"),
    make_option(c('-p', "--path"), type="character", default=NULL, help="temp directory with coverage files")
)

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

read_safe_no_header <- function(file, col.names = NULL, sep = "\t", ...) {
    if (!file.exists(file)) {
        stop("File does not exist: ", file)
    }
    line_count <- length(readLines(file, n = 1))
    if (line_count == 0) {
        if (is.null(col.names)) {
            return(data.frame())
        } else {
            return(as.data.frame(setNames(replicate(length(col.names), character(0), simplify = FALSE), col.names)))
        }
    } else {
        return(read.table(file, header = FALSE, sep = sep, ...))
    }
}

# Read input sites
dat = read.table(opt$input, sep='\t')

# Read breakpoint coverage results
sr_le = read_safe_no_header(paste(opt$path, 'le_bp.vs.SR', sep='/'))
sd_le = read_safe_no_header(paste(opt$path, 'le_bp.vs.SD', sep='/'))
rm_le = read_safe_no_header(paste(opt$path, 'le_bp.vs.RM', sep='/'))

sr_ri = read_safe_no_header(paste(opt$path, 'ri_bp.vs.SR', sep='/'))
sd_ri = read_safe_no_header(paste(opt$path, 'ri_bp.vs.SD', sep='/'))
rm_ri = read_safe_no_header(paste(opt$path, 'ri_bp.vs.RM', sep='/'))

# Default all variants to "US"
dat[, ncol(dat)+1] = 'US'

# Assign regions with priority: RM < SD < SR
if (nrow(rm_le) > 0) {
    if (nrow(dat[dat[,4] %in% rm_le[,4],]) > 0) {
        dat[dat[,4] %in% rm_le[,4],][, ncol(dat)] = 'RM'
    }
}
if (nrow(rm_ri) > 0) {
    if (nrow(dat[dat[,4] %in% rm_ri[,4],]) > 0) {
        dat[dat[,4] %in% rm_ri[,4],][, ncol(dat)] = 'RM'
    }
}

if (nrow(sd_le) > 0) {
    if (nrow(dat[dat[,4] %in% sd_le[,4],]) > 0) {
        dat[dat[,4] %in% sd_le[,4],][, ncol(dat)] = 'SD'
    }
}
if (nrow(sd_ri) > 0) {
    if (nrow(dat[dat[,4] %in% sd_ri[,4],]) > 0) {
        dat[dat[,4] %in% sd_ri[,4],][, ncol(dat)] = 'SD'
    }
}

if (nrow(sr_le) > 0) {
    if (nrow(dat[dat[,4] %in% sr_le[,4],]) > 0) {
        dat[dat[,4] %in% sr_le[,4],][, ncol(dat)] = 'SR'
    }
}
if (nrow(sr_ri) > 0) {
    if (nrow(dat[dat[,4] %in% sr_ri[,4],]) > 0) {
        dat[dat[,4] %in% sr_ri[,4],][, ncol(dat)] = 'SR'
    }
}

colnames(dat)[ncol(dat)] = 'GC'
colnames(dat)[1] = 'CHROM'
colnames(dat)[4] = 'ID'
colnames(dat)[7] = 'POS'
colnames(dat)[8] = 'REF'
colnames(dat)[9] = 'ALT'

# Override with large CNV body coverage if available
if (file.exists(paste(opt$path, 'lg_cnv.vs.SR', sep='/'))) {
    lg_cnv_sr = read.table(paste(opt$path, 'lg_cnv.vs.SR', sep='/'))
    lg_cnv_sd = read.table(paste(opt$path, 'lg_cnv.vs.SD', sep='/'))
    lg_cnv_rm = read.table(paste(opt$path, 'lg_cnv.vs.RM', sep='/'))
    lg_cnv = cbind(lg_cnv_sr[, c(4, ncol(lg_cnv_sr))], lg_cnv_sd[, ncol(lg_cnv_sd)], lg_cnv_rm[, ncol(lg_cnv_rm)])
    lg_cnv[, ncol(lg_cnv)+1] = 1 - rowSums(lg_cnv[, c(2:4)])
    colnames(lg_cnv) = c('ID', 'SR', 'SD', 'RM', 'US')
    lg_cnv[, ncol(lg_cnv)+1] = 'US'

    if (nrow(lg_cnv[lg_cnv$RM > .5,]) > 0) {
        lg_cnv[lg_cnv$RM > .5,][, ncol(lg_cnv)] = 'RM'
    }
    if (nrow(lg_cnv[lg_cnv$SD > .5,]) > 0) {
        lg_cnv[lg_cnv$SD > .5,][, ncol(lg_cnv)] = 'SD'
    }
    if (nrow(lg_cnv[lg_cnv$SR > .5,]) > 0) {
        lg_cnv[lg_cnv$SR > .5,][, ncol(lg_cnv)] = 'SR'
    }

    dat = merge(dat, lg_cnv[, c(1,6)], by='ID', all=T)
    dat[!is.na(dat[, ncol(dat)]),][, ncol(dat)-1] = dat[!is.na(dat[, ncol(dat)]),][, ncol(dat)]
}

write.table(dat[, c('CHROM', 'POS', 'REF', 'ALT', 'ID', 'GC')], opt$output, quote=F, sep='\t', col.names=F, row.names=F)
