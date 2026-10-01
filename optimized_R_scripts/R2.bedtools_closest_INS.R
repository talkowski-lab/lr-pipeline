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
		af_cor = cor(tmp[,9],tmp[,"INS_dis"])
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
		af_cor = cor(tmp[,9],tmp[,""])
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
			af_cor = cor(tmp[,9],tmp[,"INS_dis"])
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

out_columns <- c('name','name.1',pop_colname,'INS_dis','INS_ratio')

dat=read.table(input_bed,sep='\t', header=T, fill = T)
# if there's no data write an empty table and exit
if (nrow(dat) == 0) {
	out_columns[c(1,2)]=c('query_svid','ref_svid')
	out2 <- data.frame(matrix(ncol = length(out_columns), nrow = 0))
	names(out2) <- out_columns
	write.table(out2, output_bed, quote=F, sep='\t', col.names=T, row.names=F)
	quit()
}

#Alba: Added this line so it doesn't break when SVLEN.1 equals to '.'
dat[,c("SVLEN", "SVLEN.1")] <- apply(dat[,c("SVLEN", "SVLEN.1")], 2, function(col) ifelse(col == ".", as.numeric(-1), as.numeric(col)))

dat[,ncol(dat)+1] =abs(dat[,8]-dat[,2])
colnames(dat)[ncol(dat)]='INS_dis'
dat[,ncol(dat)+1] = dat[,13]/dat[,6]
colnames(dat)[ncol(dat)]='INS_ratio'
dat[!is.finite(dat[,"INS_ratio"]),"INS_ratio"] = NA
dat[,ncol(dat)+1] = abs(log(dat[,"INS_ratio"]))
colnames(dat)[ncol(dat)]='INS_ratio_dist'
dat[!is.finite(dat[,"INS_ratio_dist"]),"INS_ratio_dist"] = Inf
# Treat malformed size ratios as non-passing so they cannot win tie-breaking.
dat[,ncol(dat)+1] = !is.na(dat[,"INS_ratio"]) & dat[,"INS_dis"]<100 & dat[,"INS_ratio"]<10 & dat[,"INS_ratio"]>.1
colnames(dat)[ncol(dat)]='match_pass'

# Vectorized replacement for the per-query do.call(rbind, lapply(unique(...)))
# selection, which was O(N * #queries). Order all rows once by the same composite
# key; the stable radix sort breaks full ties by original row order, exactly like
# the per-query order()[1,]. The first row of each query group is then its winner.
# Restore the original unique(dat[,4]) appearance order so downstream sorting (and
# any coordinate-tie ordering) is identical.
o = order(!dat[,"match_pass"], dat[,"INS_dis"], dat[,"INS_ratio_dist"])
dat_o = dat[o,]
out = dat_o[!duplicated(dat_o[,4]),]
out = out[order(match(out[,4], unique(dat[,4]))),]
out=out[order(out[,3]),]
out=out[order(out[,2]),]
out=out[order(out[,1]),]

out2 = out[,out_columns]
colnames(out2)[c(1,2)]=c('query_svid','ref_svid')

out2=out2[out2$INS_dis<100 & out2$INS_ratio<10 & out2$INS_ratio>.1,]
write.table(out2, output_bed, quote=F, sep='\t', col.names=T, row.names=F)
