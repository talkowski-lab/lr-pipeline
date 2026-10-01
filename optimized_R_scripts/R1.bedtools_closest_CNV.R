#!R
## Compare CNVs
#!/usr/bin/env Rscript
library("optparse")

option_list = list(
        make_option(c("-i", "--inputbed"), type="character", default=NULL,
              help="input bed from bedtools closest", metavar="character"),
        make_option(c("-o", "--outputbed"), type="character", default=NULL,
              help="name of output bed for each query SVID with closest reference svid and information", metavar="character"),
        make_option(c("-p","--population"), type="character", default=NULL,
              help="file indluding population AF to be added", metavar="character" )
               );


#input_bed = 'bedtools_closest.DEL.bed'
#af_bed = 'CMC_V2.cleaned_filters_qual_recalibrated.AF'
#output_bed = 'bedtools_closest.DEL.selected.bed'
#output_svid = 'bedtools_closest.DEL.selected.svid'

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

input_bed = opt$inputbed
#af_bed = opt$afinfo
output_bed = opt$outputbed
pop_file = opt$population

modify_gnomad_AF<-function(af){
	af_list = strsplit(as.character(af), ',')[[1]]
	out_af=0
	if(length(af_list)>3){
		for(i in af_list[c(1:2,4:length(af_list))]){
			out_af=out_af+as.double(i)
		}
	}
	else if(length(af_list)==3){
		for(i in af_list[c(1:2)]){
			out_af=out_af+as.double(i)
		}
	}
	else if(length(af_list)==2){
		out_af	= out_af+as.double(af_list[1])
	}
	else if(length(af_list)==1){
		out_af	= af
	}
	return(out_af)
	}
AF_cor_vs_RO<-function(out){
	stat_out=data.frame('RO_cff'=0,'AF_cor'=0,'count_SV'=0, 'total_SV'=0)
	RO_range = seq(0,1, by=.05)
	rec=0
	for(i in RO_range){
		tmp = out[!out[,20]<i,]
		af_cor = cor(tmp[,9],tmp[,21])
		rec=rec+1
		stat_out[rec, 1] = i
		stat_out[rec, 2] = af_cor
		stat_out[rec, 3] = nrow(tmp)
		stat_out[rec, 4] = nrow(out)
	}
	return(stat_out)
	}
AF_cor_vs_bp<-function(out){
	stat_out=data.frame('bp_cff'=0,'AF_cor'=0,'count_SV'=0,'total_SV'=0)
	RO_range = seq(0,1000, by=50)
	rec=0
	for(i in RO_range){
		tmp = out[!out$max_bp_dis>i,]
		af_cor = cor(tmp[,9],tmp[,21])
		rec=rec+1
		stat_out[rec, 1] = i
		stat_out[rec, 2] = af_cor
		stat_out[rec, 3] = nrow(tmp)
		stat_out[rec, 4] = nrow(out)
	}
	return(stat_out)
	}
AF_cor_vs_RO_and_bp<-function(out){
	stat_out=data.frame('RO_cff'=0,'bp_cff'=0,'AF_cor'=0,'count_SV'=0,'total_SV'=0)
	RO_range = seq(0,1, by=.05)
	bp_range = seq(0,1000,by=50)
	rec=0
	for(i in RO_range){
		for(j in bp_range){
			#print(c(i,j))
			rec=rec+1
			tmp = out[!out$RO<i & !out$max_bp_dis>j,]
			af_cor = cor(tmp[,9],tmp[,21])
			stat_out[rec, 1] = i
			stat_out[rec, 2] = j
			stat_out[rec, 3] = af_cor
			stat_out[rec, 4] = nrow(tmp)
			stat_out[rec, 5] = nrow(out)
			}}
	return(stat_out)
	}
add_SV_Size<-function(chs){
  chs[,ncol(chs)+1]=0
  colnames(chs)[ncol(chs)]='size_cate'
  chs[!chs[,3]-chs[,2]<5000,][,ncol(chs)]='>5Kb'
  chs[chs[,3]-chs[,2]<5000 & !chs[,3]-chs[,2]<1000,][,ncol(chs)]='1Kb-5Kb'
  chs[chs[,3]-chs[,2]<1000 & !chs[,3]-chs[,2]<500,][,ncol(chs)]='500bp-1Kb'
  chs[chs[,3]-chs[,2]<500,][,ncol(chs)]='<500bp'
  return(chs)
	}

pop=read.table(pop_file)
pop_colname = paste(pop[,1],'AF',sep='_')
pop_colname[pop_colname=='ALL_AF']='AF'

out_columns <- c('name','name.1',pop_colname,'Reciprocal_Overlap')

dat=read.table(input_bed,sep='\t', header=T, fill = T)
# if there's no data write an empty table and exit
if (nrow(dat) == 0) {
	out_columns[c(1,2)]=c('query_svid','ref_svid')
	out2 <- data.frame(matrix(ncol = length(out_columns), nrow = 0))
	names(out2) <- out_columns
	write.table(out2, output_bed, quote=F, sep='\t', col.names=T, row.names=F)
	quit()
}

dat[,c("SVLEN", "SVLEN.1")] <- apply(dat[,c("SVLEN", "SVLEN.1")], 2, function(col) ifelse(col == ".", as.numeric(NA), as.numeric(col)))

# Vectorized replacement for the three per-row apply() calls below. We reuse the
# exact same as.matrix() character coercion apply() performed (so integer parsing,
# including scientific-notation edge cases on the numeric SVLEN columns, is
# byte-identical), applied column-wise instead of row-by-row.
mm = as.matrix(dat)
dat[,ncol(dat)+1] = pmax(as.integer(mm[,2]), as.integer(mm[,8]))
dat[,ncol(dat)+1] = pmin(as.integer(mm[,3]), as.integer(mm[,9]))
dat[,ncol(dat)+1] = dat[,ncol(dat)]-dat[,ncol(dat)-1]
colnum = ncol(dat)
{
	num = as.integer(dat[,colnum])
	a = abs(as.integer(mm[,6])); b = abs(as.integer(mm[,13]))
	den = pmax(a, b, na.rm=TRUE)         # returns the non-NA value when exactly one is NA
	den[is.na(a) & is.na(b)] = -Inf      # matches max(c(NA,NA), na.rm=TRUE) == -Inf
	dat[,ncol(dat)+1] = num/den
}
colnames(dat)[ncol(dat)]='Reciprocal_Overlap'
dat[!is.finite(dat[,"Reciprocal_Overlap"]),"Reciprocal_Overlap"] = NA
dat[,ncol(dat)+1] = pmax(abs(dat[,2] - dat[,8]), abs(dat[,3] - dat[,9]), na.rm=TRUE)
colnames(dat)[ncol(dat)]='Breakpoint_Distance'
dat[,ncol(dat)+1] = abs(log(abs(dat[,"SVLEN.1"])/abs(dat[,"SVLEN"])))
colnames(dat)[ncol(dat)]='Size_Ratio_Dist'
dat[!is.finite(dat[,"Size_Ratio_Dist"]),"Size_Ratio_Dist"] = Inf
dat[,ncol(dat)+1] = !is.na(dat[,"Reciprocal_Overlap"]) & dat[,"Reciprocal_Overlap"]>.5
colnames(dat)[ncol(dat)]='match_pass'

# Vectorized replacement for the per-query do.call(rbind, lapply(unique(...)))
# selection, which was O(N * #queries). Order all rows once by the same composite
# key; the stable radix sort breaks full ties by original row order, exactly like
# the per-query order()[1,]. The first row of each query group is then its winner.
# Restore the original unique(dat[,4]) appearance order so downstream sorting (and
# any coordinate-tie ordering) is identical.
o = order(!dat[,"match_pass"], -dat[,"Reciprocal_Overlap"], dat[,"Breakpoint_Distance"], dat[,"Size_Ratio_Dist"])
dat_o = dat[o,]
out = dat_o[!duplicated(dat_o[,4]),]
out = out[order(match(out[,4], unique(dat[,4]))),]
out=out[order(out[,3]),]
out=out[order(out[,2]),]
out=out[order(out[,1]),]

out2 = out[,out_columns]
colnames(out2)[c(1,2)]=c('query_svid','ref_svid')
out2=out2[out2$Reciprocal_Overlap>.5,]
write.table(out2, output_bed, quote=F, sep='\t', col.names=T, row.names=F)

