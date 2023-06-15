i=1
infile=list()
for(file in list.files('/home/shared/genesets/fidelity_cell_types/hifi100')){
	infile[[i]]=read.csv(paste0('/home/shared/genesets/fidelity_cell_types/hifi100/',file), header=F)
	names(infile)[i]=gsub('_.*', '', file)
	i=i+1
}
infileDf=as.data.frame(infile)
colnames(infileDf)=names(infile)
expr=read.csv('~/opentrons_testing_run/06_expression_matrix/rpm_log2_top50_25981_features.csv')
expr=data.frame(Gene=expr[,1], log2(expr[,-c(1,2)]+1))
outL=list()
for(colN in seq(1,ncol(infileDf))){
	outL[[colN]]=cor(t(expr[which(expr$Gene %in% infileDf[,colN]),-1]))
}
summary(unlist(outL[[1]][upper.tri(outL[[1]])]))
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
-0.7206  0.1855  0.4414  0.3941  0.6683  0.9545 

i=1
infile=list()
for(file in list.files('/home/shared/genesets/fidelity_cell_types/hifi100')){
	infile[[i]]=read.csv(paste0('/home/shared/genesets/fidelity_cell_types/hifi100/',file), header=F)
	names(infile)[i]=gsub('_.*', '', file)
	i=i+1
}
infileDf=as.data.frame(infile)
colnames(infileDf)=names(infile)
expr=read.csv('~/opentrons_testing_run/07_sample_networks/opentrons_sn_SampleNetworks/all_11-59-11/opentrons_sn_all_23_Qnorm.csv')
expr=data.frame(Gene=expr[,1], log2(expr[,-c(1,2)]+1))
outL=list()
for(colN in seq(1,ncol(infileDf))){
	outL[[colN]]=cor(t(expr[which(expr$Gene %in% infileDf[,colN]),-1]))
}
names(outL)=colnames(infileDf)
summary(unlist(outL[[7]][upper.tri(outL[[7]])]))


expr=read.csv('~/bdata/SF10711/rna.seq/Sams_RNA_Seq/Normalized_read_counts_using_RUVg_ERCC_K20Factors.csv')
expr=read.csv('~/bdata/SF10711/rna.seq/raw.counts/RNAseq_non-normalized-hg19-ERCC92_trimmed_reads_stranded_8lanes_featurecounts_Q1_s2.csv')
expr=data.frame(Gene=expr[,1], expr[,-seq(1,6)])
cSum=apply(expr[,-c(1)], 2, sum)/1E6
expr[,-c(1)]=expr[,-c(1)]/cSum
expr=data.frame(Gene=expr[,1], log2(expr[,-1]+1))
expr=expr[,c(1, sample(seq(2,ncol(expr)), 26, replace=F))]
outL=list()
for(colN in seq(1,ncol(infileDf))){
	outL[[colN]]=cor(t(expr[which(expr$Gene %in% infileDf[,colN]),-1]))
}
summary(unlist(outL[[7]][upper.tri(outL[[7]])]))
Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
-0.6171  0.3523  0.5858  0.5258  0.7487  0.9791 

# no ERCC cor etc.
    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
-0.5518  0.1935  0.3680  0.3469  0.5227  0.9604 


expr=read.csv('~/bdata/SF10711/rna.seq/raw.counts/RNAseq_non-normalized-hg19-ERCC92_trimmed_reads_stranded_8lanes_featurecounts_Q1_s2.csv')
summary(unlist(as.numeric(cor(expr[grep('ERCC-00', expr$Geneid), seq(7,102)]))))
erccPC1=prcomp(t(expr[grep('ERCC-00', expr$Geneid), seq(7,102)]))
erccCor=as.numeric(cor(t(expr[,seq(7,102)]), erccPC1$x[,1]))
names(erccCor)=expr$Geneid

  GPX4        COPRS  ZSCAN16-AS1        HSPB1         FIS1   SERF2      TMEM256        PPDPF   ERCC-00085        SURF1        RAB13         CIB1 
must highly correlated gens to ERCC controls

# plan is to use protein coding only genes with genes above with ERCC correciton. Not usefull to try protein coding only because already measured correlation and it is quite por with quant norm, as expected. if this fails, will need to consider repeating with ERCC controls. Will also need to derive QC metrics as spelled out in my presentation.

expr=read.csv('~/opentrons_testing_run/06_expression_matrix/counts_all_features.csv')
expr=expr[,-seq(2,5)]
zeros=apply(expr[,-1], 1, function(x) length(which(x==0)))
expr=expr[(which(zeros<6)),]
genes=expr[,1]
normG=expr[which(expr$Gene %in% c('HSPB1', 'PPDPF', 'RAB13')),]
rownames(normG)=normG$Gene
normC=cor(t(normG[,-1]))
library('RUVSeq')
setwd('/home/patrick/opentrons_testing_run/07_sample_networks')
pdf('ruvg_test.pdf')
for(i in c(1,2,3)){
	exprRUV=RUVg(round(as.matrix(expr[,-1])), which(expr$Gene %in% c('HSPB1', 'PPDPF', 'RAB13')), k=i)
	plotRLE(exprRUV$normalizedCounts, outline=T)
}
dev.off()
expr=RUVg(round(as.matrix(expr[,-1])), which(expr$Gene %in% c('HSPB1', 'PPDPF', 'RAB13')), k=3)
expr=data.frame(Gene=genes, expr$normalizedCounts)
i=1
infile=list()
for(file in list.files('/home/shared/genesets/fidelity_cell_types/hifi100')){
	infile[[i]]=read.csv(paste0('/home/shared/genesets/fidelity_cell_types/hifi100/',file), header=F)
	names(infile)[i]=gsub('_.*', '', file)
	i=i+1
}
infileDf=as.data.frame(infile)
colnames(infileDf)=names(infile)
outL=list()
for(colN in seq(1,ncol(infileDf))){
	outL[[colN]]=cor(t(expr[which(expr$Gene %in% infileDf[,colN]),-c(1,2,3)]))
}
summary(unlist(outL[[1]][upper.tri(outL[[1]])]))
#cSum=apply(expr[,-c(1)], 2, sum)/1E6
#expr[,-c(1)]=expr[,-c(1)]/cSum
expr[,-c(1)]=log2(expr[,-c(1)]+1)
outL=list()
for(colN in seq(1,ncol(infileDf))){
	outL[[colN]]=cor(t(expr[which(expr$Gene %in% infileDf[,colN]),-1]))
}
summary(unlist(outL[[1]][upper.tri(outL[[1]])]))
expr=aggregate(expr[,-c(1)], by=list(expr$Gene), FUN=sum)
