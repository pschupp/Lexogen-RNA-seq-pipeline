setwd('~/@patrick/opentrons_seq/nov_2021_run/05_feature_count')
i=1
for(sample in list.files('~/@patrick/opentrons_seq/nov_2021_run/05_feature_count')[grep('count_matrix.tsv$', list.files())]){
    temp=fread(sample, sep='\t')
    if(i==1){out=temp} else{
        out=cbind(out, temp[,ncol(temp), with=F])
    }
    i=i+1
}

gtf=fread('/home/shared/hg_align_db/GRCm39_gencode_primary/gencode.vM27.primary_assembly.annotation.ercc.phix.gtf', skip=5)
colnames(gtf)=c("chr", "db", "type", "start", "end", "score", "strand", "frame", "extra")
gtf=gtf[which(gtf$type=="gene"),]
gtf.ENSM=gsub(".*gene_id\\s\"\\\"*|\\..*", "", gtf$extra)
gtf.type=gsub(".*gene_type\\s\"\\\"*|\".*", "", gtf$extra)
gtf.name=gsub(".*gene_name\\s\"\\\"*|\".*", "", gtf$extra)
gtf.out=data.frame(gtf[,1], gtf.ENSM, gtf.type, gtf.name)
gtf.out$gtf.ENSM=gsub("\";.*", "", gtf.out$gtf.ENSM)
gtf.out$gtf.name[grep('^ERCC-[0-9]',gtf.out$gtf.ENSM)]=gtf.out$gtf.ENSM[grep('^ERCC-[0-9]',gtf.out$gtf.ENSM)]

out=data.table(Gene=gtf.out$gtf.name[match(gsub('\\..*', '', out$Geneid), gtf.out$gtf.ENSM)], out)
colnames(out)=gsub('.*\\/', '', gsub('_001\\..*', '', colnames(out)))
#summary(apply(out[,-seq(1,7)], 2, sum))
#
# perform ruvg regression out for 0 up to 20 factors
# will measure before and after regression of the first n factors:
#   a) measure correlation of AOMN genesets from Ben Barres mouse-specific gensets
#   b) plot RLE for 20 samples, evenly distributed
#   c) overall correlation structure of the data
# will make a, b, and c a function that can be applied after each modification

# remove genes that are not expressed in at least 10% sections
out=out[-which(apply(out[,-seq(1,7)], 1, function(x) length(which(x==0)))>5),]
outMin=data.frame(out[,-seq(2,7)])

reg_eval=function(df, plotName){
    cellMarkers=list(as.character(read.table('/home/shared/genesets/OurSets/MOSET7.csv')[1:10,1]),as.character(read.table('/home/shared/genesets/OurSets/MOSET8.csv')[1:10,1]), as.character(read.table('/home/shared/genesets/OurSets/MOSET6805.csv')[1:10,1]),as.character(read.table('/home/shared/genesets/OurSets/MOSET9.csv')[1:10,1]))
    names(cellMarkers)=c('astro', 'oligo', 'micro', 'neuron')
    cellCors=lapply(cellMarkers, function(x) cor(t(df[grep(paste(x, collapse='|'), df$Gene, ignore.case=T),])))
    cellCorsDistro=lapply(cellCors, function(x) summary(x[upper.tri(x)]))
    return(cellCorsDistro)
    pdf(plotName)
    for(i in seq(1, length(cellCors))){
        plot=data.frame(Cor=cellCors[[i]][upper.tri(cellCors[[i]])])
        print(ggplot(plot, aes(x=Cor))+
              geom_histogram()+
              labs(title=names(cellCors)[i]))
    }
    dfRLE=apply(df[,-1], 2, function(x) log10((x/(median(x, na.rm=T))+1)))
    dfRLEPlot=melt(dfRLE[, seq(1, ncol(dfRLE), length=20)])[,-1]
	print(ggplot(dfRLEPlot, aes(x=Var2, y=value))+
          geom_boxplot(outlier.shape=NA, notch=F)+
          #labs(title=paste('RUVg RLE plot, k='), subtitle=paste('Median:', signif(median(unlist(dfin)), 4), ', Median variance: ', signif(median(apply(dfin, 2, var)),4)), x='Nuclei', y='RLE (log count/median)')+
          theme_classic()+
          theme(legend.position="none", axis.text.x=element_blank())
      )
    dfCor=cor(t(df[sample(seq(1,nrow(df)), 2000),-1]))
    dfCorM=melt(dfCor[upper.tri(dfCor)])
    print(ggplot(dfCorM, aes(x=value))+
        geom_histogram()+
        labs(title="Correlation distribution"))
    dev.off()
}
library('RUVSeq')
library('RColorBrewer')
eval0=reg_eval(outMin, 'input_data.pdf')
eval0=t(do.call(rbind.data.frame, eval0))
colnames(eval0)=c('astro', 'oligo', 'micro', 'neuron')
rownames(eval0)=c('min', '1quart', 'median', 'mean', '3quart', 'max')

spikes=rep(FALSE, nrow(outMin))
spikes[grep("^ERCC-", outMin$Gene)]=TRUE

exprRun_1=RUVg(as.matrix(round(outMin[,-1])), spikes, k=1, round=F, center=TRUE, isLog=F)$normalizedCounts
exprRun_1=data.frame(Gene=outMin$Gene, exprRun_1)
eval1=reg_eval(exprRun_1, 'k1.pdf')
eval1=t(do.call(rbind.data.frame, eval1))
colnames(eval1)=c('astro', 'oligo', 'micro', 'neuron')
rownames(eval1)=c('min', '1quart', 'median', 'mean', '3quart', 'max')


exprRun_2=RUVg(as.matrix(round(outMin[,-1])), spikes, k=2, round=F, center=TRUE, isLog=F)$normalizedCounts
exprRun_2=data.frame(Gene=outMin$Gene, exprRun_2)
eval2=reg_eval(exprRun_2, 'k2.pdf')
eval2=t(do.call(rbind.data.frame, eval2))
colnames(eval2)=c('astro', 'oligo', 'micro', 'neuron')
rownames(eval2)=c('min', '1quart', 'median', 'mean', '3quart', 'max')

exprRun_3=RUVg(as.matrix(round(outMin[,-1])), spikes, k=3, round=F, center=TRUE, isLog=F)$normalizedCounts
exprRun_3=data.frame(Gene=outMin$Gene, exprRun_3)
eval3=reg_eval(exprRun_3, 'k3.pdf')
eval3=t(do.call(rbind.data.frame, eval3))
colnames(eval3)=c('astro', 'oligo', 'micro', 'neuron')
rownames(eval3)=c('min', '1quart', 'median', 'mean', '3quart', 'max')

exprRun_4=RUVg(as.matrix(round(outMin[,-1])), spikes, k=4, round=F, center=TRUE, isLog=F)$normalizedCounts
exprRun_4=data.frame(Gene=outMin$Gene, exprRun_4)
eval4=reg_eval(exprRun_4, 'k4.pdf')
eval4=t(do.call(rbind.data.frame, eval4))
colnames(eval4)=c('astro', 'oligo', 'micro', 'neuron')
rownames(eval4)=c('min', '1quart', 'median', 'mean', '3quart', 'max')

exprRun_5=RUVg(as.matrix(round(outMin[,-1])), spikes, k=5, round=F, center=TRUE, isLog=F)$normalizedCounts
exprRun_5=data.frame(Gene=outMin$Gene, exprRun_5)
eval5=reg_eval(exprRun_5, 'k5.pdf')
eval5=t(do.call(rbind.data.frame, eval5))
colnames(eval5)=c('astro', 'oligo', 'micro', 'neuron')
rownames(eval5)=c('min', '1quart', 'median', 'mean', '3quart', 'max')

# remove either 1 or 2 factors
exprRun_1=exprRun_1[,-1]
reads.sum=apply(exprRun_1, 2, sum)
for(ea in seq(1,ncol(exprRun_1))){
	exprRun_1[,ea]=exprRun_1[,ea]/(reads.sum[ea]/1E6)
}
exprRun_1=data.frame(Gene=outMin$Gene, exprRun_1)
exprRun_1=exprRun_1[-grep('^ERCC-', exprRun_1$Gene),]
# exprRun_1=log2(exprRun_1+1)
exprRun_1Bak=aggregate(exprRun_1[,-1], by=list(exprRun_1$Gene), FUN=sum)
colnames(exprRun_1Bak)[1]='Gene'
source("/home/shared/code/FindModules/FindModules094.R")
setwd('/mnt/bdata/@patrick/opentrons_seq/nov_2021_run/06_expression_analysis/')
FindModules(
	projectname="opentrons_modules_k=1",
	expr=exprRun_1Bak,
	geneinfo=c(1),
	sampleindex=c(2:ncol(exprRun_1)),
	samplegroups=as.factor(c(substr(colnames(exprRun_1Bak)[-1], 1,2))),
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
	signumvec=c(0.999, 0.99, 0.95, 0.9, 0.8, 0.7,0.6),
	minsizevec=c(15, 10, 8, 5, 3),
	minMEcor=0.85,
	signum=NULL,
	minSize=NULL,
	ZNCcut=2,
	calcSW=FALSE,
	loadTree=FALSE,
	writeKME=TRUE,
	calcBigModStat=FALSE
)
# next FindModules

~/@patrick/opentrons_seq/nov_2021_run/06_expression_analysis/opentrons_modules_k=1_Modules/Bicor-None_signum0.354_minSize8_minMEcor0.85_18537

source('~/code/git/GSEA/GSEAfxsV3.r')
MyGSHGloop(kmecut1='fdr',exclude="none",pvalcut1=NULL, moduleDir='~/@patrick/opentrons_seq/nov_2021_run/06_expression_analysis/opentrons_modules_k=1_Modules', geneSets='/home/patrick/code/GSEA/genesets_slim', outputName='slim_sets_fdr')
MyGSHGloop(kmecut1='bc',exclude="none",pvalcut1=NULL, moduleDir='~/@patrick/opentrons_seq/nov_2021_run/06_expression_analysis/opentrons_modules_k=1_Modules', geneSets='/home/patrick/code/GSEA/genesets_slim', outputName='slim_sets_bc')

source('~/code/git/GSEA/GSEAfxsV3.r')
MyGSHGloop(kmecut1='fdr',exclude="none",pvalcut1=NULL, moduleDir='~/@patrick/opentrons_seq/nov_2021_run/06_expression_analysis/opentrons_modules_k=1_Modules', geneSets='/home/patrick/code/GSEA/genesets_slim', outputName='broad_fdr')
MyGSHGloop(kmecut1='bc',exclude="none",pvalcut1=NULL, moduleDir='~/@patrick/opentrons_seq/nov_2021_run/06_expression_analysis/opentrons_modules_k=1_Modules', geneSets='/home/patrick/code/GSEA/genesets_slim', outputName='broad_bc')
