expr=read.csv('counts_all_features.csv')
expr=expr[,-c(2,3,4,5)]
sums=apply(expr[,-1], 2, sum)
cors=cor(t(expr[,-1]), sums)
names(cors)=expr$Gene
sort(cors, decreasing=T)[1:10]

zeros=apply(expr[,-1], 1, function(x) length(which(x==0)))
expr=expr[(which(zeros<6)),]
vars=apply(expr[,-1], 1, var)
means=apply(expr[,-1], 1, mean)

out=data.frame(means, vars)
out=out[(out$vars<1E4),]
pdf('var_mean.pdf')
plot(out)
dev.off()


# FindModules
expr=read.csv('~/opentrons_testing_run/06_expression_matrix/rpm_protein_coding_features.csv')
setwd('~/opentrons_testing_run/07_sample_networks')
expr=aggregate(expr[,-c(1,2)], by=list(expr$Gene), FUN=sum)
colnames(expr)[1]='Gene'
expr=expr[,-c(2,3,4)]
zeros=apply(expr[,-1], 1, function(x) length(which(x==0)))
expr=expr[(which(zeros<6)),]
source('~/code/git/FindModules/FindModules.R')
setwd('~/opentrons_testing_run/07_sample_networks')
FindModules(
 projectname="sample_networks_ruvg",
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
# Our sets
modDir='~/opentrons_testing_run/07_sample_networks/no_sample_networks_Modules'
setwd(modDir)
source('~/code/GSEA/GSEA_fxs/GSEAfxsV3.r')
MyGSHGloop(kmecut1='topmodposbc',exclude="none",pvalcut1=NULL)
setwd(modDir)
MyGSHGloop(kmecut1='topmodposfdr',exclude="none",pvalcut1=NULL)
# Broad sets
source('~/code/git/GSEA/GSEAfxsV3.r')
broadSets=getBroadSets("/home/shared/genesets/Broad_GSEA/v7/msigdb_v7.4.xml")
print("Success reading in 'broadSets'")
modDir='~/opentrons_testing_run/07_sample_networks/no_sample_networks_Modules'
setwd(modDir)
## To run enrichment analysis for broad gene sets in all networks:
### Note that kmecut1 can equal "seed", "topmodposbc" (recommended), or "topmodposfdr".
print("Beginning loop BC")
BroadGSHGloop(kmecut1="topmodposbc",pvalcut1=NULL)
setwd(modDir)
print("Beginning loop FD")
BroadGSHGloop(kmecut1="topmodposfdr",pvalcut1=NULL)
setwd(modDir)

source('~/code/git/GSEA/GSEAfxsV3.r')
modDir='/home/patrick/opentrons_testing_run/07_sample_networks/no_sample_networks_protein_coding_Modules'
MyGSHGloop(kmecut1='topmodposbc',exclude="none",pvalcut1=NULL, moduleDir=modDir, geneSets='~/code/GSEA/genesets_slim', outputName='with_sn_bc')
MyGSHGloop(kmecut1='topmodposfdr',exclude="none",pvalcut1=NULL, moduleDir=modDir, geneSets='~/code/GSEA/genesets_slim', outputName='with_sn_fdr')

# in R
library('future')
alignDir='~/opentrons_testing_run/04_feature_count/'
outDir='~/opentrons_testing_run/05_no_dedup_counts/'
plan(multicore, workers=15)
for(name in list.files(alignDir)[grep('001.align.bamAligned.sortedByCoord.out.bam.featureCounts.sam', list.files(alignDir))]){
	temp=future({
	alig=fread(paste0(alignDir,name), skip='@CO', sep='\t', select=c(1,2,3,4,17,18), fill=T)
	colnames(alig)=paste0('V', c(1,2,3,4,17,18))
	alig$V1=gsub('.*:', '', alig$V1)
	alig$V17=gsub('XN:i:', '', alig$V17, fixed=T)
	alig=alig[alig$V17>0,]
	alig=alig[alig$V17<5,]
	alig$V18=gsub('XT:Z:', '', alig$V18, fixed=T)
	
	doubles=alig[alig$V17==2,]
	doubles$V17=as.numeric(doubles$V17)
	doubles$V17=0.5
	doublesNames=strsplit(doubles$V18, ',')
	doublesNames=unlist(c(lapply(doublesNames, function(x) x[[1]]),lapply(doublesNames, function(x) x[[2]])))
	doubles=rbind(doubles,doubles)
	doubles$V18=doublesNames

	triples=alig[alig$V17==3,]
	triples$V17=as.numeric(triples$V17)
	triples$V17=0.33
	triplesNames=strsplit(triples$V18, ',')
	triplesNames=unlist(c(lapply(triplesNames, function(x) x[[1]]),lapply(triplesNames, function(x) x[[2]]),lapply(triplesNames, function(x) x[[3]])))
	triples=rbind(triples,triples,triples)
	triples$V18=triplesNames
	quads=alig[alig$V17==4,]
	quads$V17=as.numeric(quads$V17)
	quads$V17=0.25
	quadsNames=strsplit(quads$V18, ',')
	quadsNames=unlist(c(lapply(quadsNames, function(x) x[[1]]),lapply(quadsNames, function(x) x[[2]]),lapply(quadsNames, function(x) x[[3]]),lapply(quadsNames, function(x) x[[4]])))
	quads=rbind(quads,quads,quads,quads)
	quads$V18=quadsNames

	aligN=rbind(alig[alig$V17==1,], doubles, triples, quads)
#	aligN=aligN[!duplicated(paste0(aligN$V1,aligN$V3, aligN$V4, aligN$V18)),]
	aligN$V17=as.numeric(aligN$V17)
	aligA=aggregate(aligN$V17, by=list(aligN$V18), FUN=sum)
	write.table(aligA, file=paste0(outDir, name, '.counts.tsv'), row.names=F, col.names=F, quote=F, sep=',')
	})
}
exprDir='~/opentrons_testing_run/05_no_dedup_counts/'
inL=list()
i=1
for(expr in list.files(exprDir)[grep('_R1_001.align.bamAligned.sortedByCoord.out.bam.featureCounts.sam.counts.tsv', list.files(exprDir))]){
	inL[[i]]=fread(paste0(exprDir, expr))
	i=i+1
}
genes=c()
genes=unique(unlist(lapply(inL, function(x) genes=c(genes, x[[1]]))))
out=data.frame(genes)
for(ea in inL){
	out=cbind(out, ea[match(genes, ea[[1]]),2])
}
out[is.na(out)]=0
colnames(out)=c('ENSG', gsub('_R1.*', '', list.files(exprDir)))
gtf=fread('/home/shared/hg_align_db/GRCh38_gencode_primary/gencode.v38.primary_assembly.annotation.gtf')
gtf=gtf[gtf$V3=='transcript',]
gtf=data.frame(ENSG=gsub('.*gene_id "', '', gsub('";.*', '', as.character(gtf$V9))),gene=gsub('";.*', '',  gsub('.*gene_name "', '', as.character(gtf$V9))),gtype=gsub('";.*', '',  gsub('.*gene_type "', '', as.character(gtf$V9))))
gtf=gtf[gtf$gtype=='protein_coding',]
out=data.frame(Gene=gtf[match(out[[1]], gtf[[1]]), 2], ENSG=out[[1]], out[,-1])
out=out[!is.na(out$Gene),]
write.csv(out, file='~/opentrons_testing_run/06_expression_matrix_no_dedup/counts_protein_coding_features.csv', row.names=F, quote=F)
cSum=apply(out[,-c(1,2)], 2, sum)/1E6
out[,-c(1,2)]=out[,-c(1,2)]/cSum
write.csv(out, file='~/opentrons_testing_run/06_expression_matrix_no_dedup/rpm_protein_coding_features.csv', row.names=F, quote=F)

# FindModules
expr=read.csv('~/opentrons_testing_run/06_expression_matrix_no_dedup/rpm_protein_coding_features.csv')
setwd('~/opentrons_testing_run/07_sample_networks')
expr=aggregate(expr[,-c(1,2)], by=list(expr$Gene), FUN=sum)
colnames(expr)[1]='Gene'
expr=expr[,-c(2,3,4)]
zeros=apply(expr[,-1], 1, function(x) length(which(x==0)))
expr=expr[(which(zeros<6)),]
source('~/code/git/FindModules/FindModules.R')
setwd('~/opentrons_testing_run/07_sample_networks')
FindModules(
 projectname="no_sample_networks_protein_coding_no_dedup",
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
