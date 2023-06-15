library('data.table')
library('reshape2')
library('ggplot2')

gtf=fread('/home/shared/hg_align_db/GRCh38_gencode_primary/gencode.v38.primary_assembly.annotation.ercc.phix.gtf')
gtf=gtf[gtf$V3=='gene',]
gtf=data.frame(
    ENSG=gsub('.*gene_id "', '', gsub('";.*', '', as.character(gtf$V9))),
    gene=gsub('gene_id "', '', gsub('";.*', '',  gsub('.*gene_name "', '', as.character(gtf$V9)))),
    gtype=gsub('gene_id "', '',gsub('";.*', '',  gsub('.*gene_type "', '', as.character(gtf$V9)))),
    length=gtf$V5-gtf$V4+1
)
gtf$gtype=gsub('(ERCC)-[0-9]*', '\\1', gtf$gtype)


WD = '/home/patrick/SF11949_run/SF11949_run_1_0.8x_t2'
exprDir=paste0(WD, '/05_feature_count/')
outDir=paste0(WD, '/07_reports/')
plan(multicore, workers=15)
# combine into one expression matrix
inL=list()
i=1
for(expr in list.files(exprDir)[grep('count_matrix.tsv$', list.files(exprDir))]){
	inL[[i]]=fread(paste0(exprDir, expr))
	i=i+1
}
genes=c()
genes=unique(unlist(lapply(inL, function(x) genes=c(genes, x[[1]]))))
out=data.frame(genes)
for(ea in inL){
	out=cbind(out, ea[match(genes, ea[[1]]),7])
}
#out=out[,-which(colnames(out)=='Chr')]
out[is.na(out)]=0
colnames(out)=c('ENSG', gsub('04_deduplicated/', '', colnames(out)[-1]))
colnames(out)=c('ENSG', gsub('\\..*', '', colnames(out)[-1]))

out=data.frame(Gene=gtf[match(out[[1]], gtf[[1]]), c(2,3,4)], ENSG=out[[1]], out[,-1])
colnames(out)[1:3]=c('Gene', 'type', 'length')
out=data.frame(out)
outAg=data.frame(length=out$length, counts=apply(out[,-seq(1,4)], 1, sum))
outAg$lengthT=as.factor(cut(outAg$length, seq(0, 25000, 500)))
outPlot=aggregate(outAg$counts, by=list(outAg$lengthT), FUN=sum)
outPlot$x=outPlot$x/sum(outPlot$x)
WD = '/home/patrick/SF11949_run/SF11949_run_1_1.0x'
exprDir=paste0(WD, '/05_feature_count/')
outDir=paste0(WD, '/07_reports/')
plan(multicore, workers=15)
# combine into one expression matrix
inL=list()
i=1
for(expr in list.files(exprDir)[grep('count_matrix.tsv$', list.files(exprDir))]){
	inL[[i]]=fread(paste0(exprDir, expr))
	i=i+1
}
genes=c()
genes=unique(unlist(lapply(inL, function(x) genes=c(genes, x[[1]]))))
out=data.frame(genes)
for(ea in inL){
	out=cbind(out, ea[match(genes, ea[[1]]),7])
}
#out=out[,-which(colnames(out)=='Chr')]
out[is.na(out)]=0
colnames(out)=c('ENSG', gsub('04_deduplicated/', '', colnames(out)[-1]))
colnames(out)=c('ENSG', gsub('\\..*', '', colnames(out)[-1]))

out=data.frame(Gene=gtf[match(out[[1]], gtf[[1]]), c(2,3,4)], ENSG=out[[1]], out[,-1])
colnames(out)[1:3]=c('Gene', 'type', 'length')
out=data.frame(out)
outAg1=data.frame(length=out$length, counts=apply(out[,-seq(1,4)], 1, sum))
outAg1$lengthT=as.factor(cut(outAg1$length, seq(0, 25000, 500)))
outPlot1=aggregate(outAg1$counts, by=list(outAg1$lengthT), FUN=sum)
outPlot1$x=outPlot1$x/sum(outPlot1$x)

outAll=data.frame(outPlot, ratio1=outPlot$x)
colnames(outAll)=c('range','ratio0.85', 'ratio1.0')
p=ggplot(outAll, aes(x=ratio0.85, y=ratio1.0))+
    geom_point() +
    labs(x='Proportional counts - 0.85x ratio', y='Proportional counts - 1.00x ratio', title='No evidence of size bias in 0.85 v 1.00 SPRI bead ratio ')

q=ggplot(outPlot, aes(x=Group.1, y=x)) +
    geom_bar(stat='identity') +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
    labs(x='Size bins (bp)', y='Proportional count', title='Distributions of sizes for 0.85x ratio', subtitle='601430 reads')

r=ggplot(outPlot1, aes(x=Group.1, y=x)) +
    geom_bar(stat='identity') +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
    labs(x='Size bins (bp)', y='Proportional count', title='Distributions of sizes for 1.00x ratio', subtitle='5664931 reads')

pdf('size_compare.pdf')
print(q)
print(r)
print(p)
dev.off()






setwd(WD)
