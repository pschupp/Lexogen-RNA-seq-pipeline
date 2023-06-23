# TODO: replace geneSet with dir from yaml file
# TODO: replace gtfFile with dir from yaml file
library('future')
library('data.table')
library('reshape2')
library('ggplot2')
WD = getwd()
exprDir=paste0(WD, '/07_feature_count/')
outDir=paste0(WD, '/11_reports/')
plan(multicore, workers=15)

test_cor=function(expr){
	if(!(exists('infileDf'))){
		i=1
		infile=list()
		for(file in list.files(geneSet)){
			infile[[i]]=read.csv(paste0(geneSet,file), header=F)
			names(infile)[i]=gsub('_.*', '', file)
			i=i+1
		}
		infileDf=as.data.frame(infile)
		colnames(infileDf)=names(infile)
	}
	outL=list()
	for(colN in seq(1,ncol(infileDf))){
		outL[[colN]]=cor(t(expr[which(expr$Gene %in% infileDf[,colN]),-c(1,2,3)]))
	}
	out=data.frame(matrix(unlist(lapply(outL, function(x) summary(unlist(x[upper.tri(x)])))), nrow=6))
	colnames(out)=colnames(infileDf)
	rownames(out)=c('Min.', '1st Qu.', 'Median', 'Mean', '3rd Qu.', 'Max.')
	print(out)
}
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

# convert ENSG ids to gene names
gtf=fread(gtfFile)
gtf=gtf[gtf$V3=='gene',]
gtf=data.frame(
    ENSG=gsub('.*gene_id "', '', gsub('";.*', '', as.character(gtf$V9))),
    gene=gsub('gene_id "', '', gsub('";.*', '',  gsub('.*gene_name "', '', as.character(gtf$V9)))),
    gtype=gsub('gene_id "', '',gsub('";.*', '',  gsub('.*gene_type "', '', as.character(gtf$V9))))
)
gtf$gtype=gsub('(ERCC)-[0-9]*', '\\1', gtf$gtype)
out=data.frame(Gene=gtf[match(out[[1]], gtf[[1]]), c(2,3)], ENSG=out[[1]], out[,-1])
colnames(out)[1:2]=c('Gene', 'type')
out=data.frame(out)
countBrake=t(data.frame(
        rRNA=apply(out[grep('^rRNA_*', out$type), -seq(1,3)], 2, sum), 
        mtRNA=apply(out[grep('^Mt_*', out$type), -seq(1,3)], 2, sum), 
        other_RNA=apply(out[grep('^scaRNA|^scRNA|^snoRNA|^sRNA|^miRNA|^lncRNA|^misc_RNA', out$type), -seq(1,3)], 2, sum), 
        protein=apply(out[grep('protein_coding', out$type), -seq(1,3)], 2, sum), 
        all=apply(out[, -seq(1,3)], 2, sum)))

countBrake=rbind(countBrake, other=(apply(out[, -seq(1,3)], 2, sum)-apply(countBrake[1:4,], 2, sum)))
countBrake=t(t(countBrake[-5,])/countBrake[5,])
plot=reshape2::melt(countBrake[c(1,2,3,5,4),])


out$type[is.na(out$type)]='PhiX'
setwd(WD)
pdf('features.pdf')
print(ggplot(plot, aes(x=Var2, y=value, fill=Var1))+
	geom_bar(stat="identity", width=0.8, size=1)+
	scale_y_continuous(expand=c(0,0), breaks=seq(0,1,0.2))+
	labs(fill='Feature type', title="Alignment by feature", subtitle='not including no feature alignments', x="Sample", y="Percentage")+
	theme_classic()+
	theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)))
dev.off()
write.csv(countBrake, file='feature_alignment_breakdown.csv')
library('RUVSeq')
genes=out[,c(1,2,3)]
expr=RUVg(round(as.matrix(out[,-c(1,2,3)])), which(out$Gene %in% c('HSPB1', 'PPDPF', 'RAB13')), k=3)
expr=data.frame(genes, expr$normalizedCounts)
write.csv(expr, 'counts_ruvg_all_features.csv', row.names=F, quote=F)
zeros=apply(expr[,-c(1,2)], 1, function(x) length(which(x==0)))
expr=expr[(which(zeros<6)),]
expr=expr[!is.na(expr$Gene),]
expr=aggregate(expr[,-c(1,2,3)], by=list(expr$Gene), FUN=sum)
expr=data.frame(Gene=expr[,1], genes[match(expr[,1], genes$Gene), c(2,3)], expr[,-1])
write.csv(expr, 'counts_ruvg_expressed_features.csv', row.names=F, quote=F)
expr=expr[expr$type=='protein_coding',]
write.csv(expr, 'counts_ruvg_protein_coding.csv', row.names=F, quote=F)
expr[,-c(1,2,3)]=normalizeQuantiles(expr[,-c(1,2,3,4,5,6)])
expr[,-c(1,2,3)]=log2(expr[,-c(1,2,3,4,5,6)]+1)
expr=read.csv('counts_ruvg_protein_coding.csv')
exprPC=prcomp(t(expr[,-c(1,2,3,4,5,6)]))
m1=apply(expr[,-c(1,2,3,4,5,6)], 1, function(x) lm(x ~ exprPC$x[,1]+exprPC$x[,2]))
m2=t(data.frame(lapply(m1, function(x) resid(x))))
rownames(m2)=seq(1,nrow(m2))
m2=data.frame(expr[,c(1,2,3)], m2)
m1=apply(expr[,-c(1,2,3,4,5,6)], 1, function(x) lm(x ~ exprPC$x[,1]+exprPC$x[,2]))
m2=t(data.frame(lapply(m1, function(x) resid(x))))
rownames(m2)=seq(1,nrow(m2))
m2=data.frame(expr[,c(1,2,3)], m2)
setwd('/home/patrick/opentrons_testing_run/06_expression_matrix')
pdf('cors.pdf')
expr=read.csv('counts_ruvg_protein_coding.csv')
hist(unlist(as.numeric(cor(t(expr[sample(seq(1,nrow(expr)), 10000),-c(1,2,3)])))))
expr=read.csv('counts_protein_coding_features.csv')
hist(unlist(as.numeric(cor(t(expr[sample(seq(1,nrow(expr)), 10000),-c(1,2,3,4,5)])))))
expr=read.csv('counts_ruvg_all_features.csv')
hist(unlist(as.numeric(cor(t(expr[sample(seq(1,nrow(expr)), 10000),-c(1,2,3,4,5,6)])))))
expr=read.csv('counts_ruvg_expressed_features.csv')
hist(unlist(as.numeric(cor(t(expr[sample(seq(1,nrow(expr)), 10000),-c(1,2,3,4,5,6)])))))
exprPC=prcomp(t(expr[,-c(1,2,3,4,5,6)]))
m1=apply(expr[,-c(1,2,3,4,5,6)], 1, function(x) lm(x ~ exprPC$x[,1]))
m2=t(data.frame(lapply(m1, function(x) resid(x))))
rownames(m2)=seq(1,nrow(m2))
m2=data.frame(expr[,c(1,2,3)], m2)
hist(unlist(as.numeric(cor(t(m2[sample(seq(1,nrow(m2)), 10000),-c(1,2,3)])))))
dev.off()

expr=read.csv('counts_ruvg_expressed_features.csv')
expr=read.csv('counts_ruvg_protein_coding.csv')
source('~/code/git/FindModules/FindModules.R')
expr=data.frame(expr[,c(1,2,3)], normalizeQuantiles(expr[,-seq(1,6)]))
pdf('quant_ruvg_prot.pdf')
hist(unlist(as.numeric(cor(t(expr[sample(seq(1,nrow(expr)), 10000),-c(1,2,3)])))))
dev.off()

FindModules(
 projectname="ruvg_prot_qn",
 expr=expr,
 geneinfo=c(1),
 sampleindex=seq(4,ncol(expr)),
 samplegroups=as.factor(seq(1,23)),
 subset=NULL,
 simMat=NULL,
 saveSimMat=FALSE,
 simType="Bicor",
 beta=1,
 overlapType="None",
 TOtype="signed",
 TOdenom="min",
 MIestimator="mi.mm",
 MIdisc="equalfreq",
 signumType="rel",
 iterate=TRUE,
 signumvec=c(0.99, 0.97, 0.95, 0.9, .8,.7,0.6),
 minsizevec=c(20, 15, 10, 8),
 signum=NULL,
 minSize=NULL,
 merge.by="ME",
 merge.param=0.8,
 export.merge.comp=T,
 ZNCcut=2,
 calcSW=FALSE,
 loadTree=FALSE,
 writeKME=TRUE,
 calcBigModStat=FALSE,
 writeModSnap=TRUE
)


expr=read.csv('~/opentrons_testing_run/06_expression_matrix/rpm_top50_25981_features.csv')
expr=expr[,-seq(3,5)]
sif=read.csv('~/opentrons_testing_run/07_sample_networks/Sequencing Master Table.csv')
sif=sif[-c(21, seq(25,30)), c(1,2,4,5,6,8)]
sif=data.frame(name=colnames(expr)[-c(1,2)], sif)
colnames(sif)=c('name', 'section', 'nd_conc', 'r260_280', 'r260_230', 'rin', 'lib_conc')
sif$section[21:23]=seq(49,51)
sif$nd_conc[21:23]=c(27.1, 21.0, 17.7)
sif=data.frame(sif, group=rep('all', nrow(sif)))
sif[,c(1,2,8)]=apply(sif[,c(1,2,8)], 2, as.factor)
sif[,seq(3,7)]=apply(sif[,seq(3,7)], 2, as.numeric)
setwd('~/opentrons_testing_run/07_sample_networks')
# sample networks
source('~/code/git/SampleNetworks/SampleNetwork_1.07_PGS_mod.R') # dev branch
SampleNetwork(
	datExprT=expr,
	method1="correlation",
	impute1=FALSE,
	exclude1="any",
	subset1=NULL,
	skip1=c(2),
	indices1=list(seq(3, ncol(expr))),
	sampleinfo1=sif,
	subgroup1=1,
	samplelabels1=1,
	grouplabels1=8,
	fitmodels1=TRUE,
	whichmodel1="univariate",
	whichfit1="pc1",
	btrait1=c(2,3,4,5,6,7),
	trait1=NULL, 
	asfactors1=c(2),
	projectname1='opentrons_sn',
	cexlabels1=1,
	normalize1=TRUE,
	replacenegs1=FALSE,
	exportfigures1=TRUE,
	verbose=TRUE
)
# consider excluding E5_S37
expr=read.csv('~/opentrons_testing_run/07_sample_networks/opentrons_sn_SampleNetworks/all_11-59-11/opentrons_sn_all_23_Qnorm.csv')
setwd('~/opentrons_testing_run/07_sample_networks')
expr=aggregate(expr[,-c(1,2)], by=list(expr$Gene), FUN=sum)
colnames(expr)[1]='Gene'
zeros=apply(expr[,-1], 1, function(x) length(which(x==0)))
expr=expr[(which(zeros<6)),]
source('~/code/git/FindModules/FindModules.R')
setwd('~/opentrons_testing_run/07_sample_networks')
FindModules(
 projectname="qn_top50",
 expr=expr,
 geneinfo=c(1),
 sampleindex=seq(2,ncol(expr)),
 samplegroups=as.factor(seq(1,23)),
 subset=NULL,
 simMat=NULL,
 saveSimMat=FALSE,
 simType="Bicor",
 beta=1,
 overlapType="None",
 TOtype="signed",
 TOdenom="min",
 MIestimator="mi.mm",
 MIdisc="equalfreq",
 signumType="rel",
 iterate=TRUE,
 signumvec=c(0.999, 0.99, 0.97, 0.95, 0.9, .8,.7),
 minsizevec=c(20, 15, 10, 8),
 signum=NULL,
 minSize=NULL,
 merge.by="ME",
 merge.param=0.85,
 export.merge.comp=T,
 ZNCcut=2,
 calcSW=FALSE,
 loadTree=FALSE,
 writeKME=TRUE,
 calcBigModStat=FALSE,
 writeModSnap=TRUE
)
