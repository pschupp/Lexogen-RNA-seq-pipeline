# deduplicate
library('future')
library('data.table')
library('reshape2')
library('ggplot2')
WD = '/home/patrick/SF11949_run/SF11949_run_1_0.8x_t2'
exprDir=paste0(WD, '/05_feature_count/')
outDir=paste0(WD, '/07_reports/')
plan(multicore, workers=15)

test_cor=function(expr){
	if(!(exists('infileDf'))){
		i=1
		infile=list()
		for(file in list.files('/home/shared/genesets/fidelity_cell_types/hifi100')){
			infile[[i]]=read.csv(paste0('/home/shared/genesets/fidelity_cell_types/hifi100/',file), header=F)
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
gtf=fread('/home/shared/hg_align_db/GRCh38_gencode_primary/gencode.v38.primary_assembly.annotation.ercc.phix.gtf')
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
test_cor(expr)
#          Astrocyte  Cholinergic      Choroid Dopaminergic Endothelial  Ependymal  Microglia
# Min.    -0.8872317 -0.824993562 -0.853883545   -0.8758288 -0.86576165 -0.7626769 -0.7799344
# 1st Qu.  0.3392508  0.007228488  0.001196181    0.1013966  0.05627019  0.1105487 -0.1463454
# Median   0.6294498  0.280882313  0.266191419    0.4163428  0.33588900  0.3738504  0.1623001
# Mean     0.5193329  0.262050933  0.236902411    0.3480505  0.28960311  0.3311612  0.1535857
# 3rd Qu.  0.7999434  0.541069145  0.496611701    0.6496634  0.57150635  0.5774853  0.4887060
# Max.     0.9858948  0.960718155  0.936418197    0.9796844  0.97700503  0.9327160  0.9800072
#              Mural     Neuron Oligodendrocyte        OPC    Purkinje
# Min.    -0.8629047 -0.8299676      -0.8613016 -0.6366964 -0.86710090
# 1st Qu. -0.1727524  0.3631321       0.1001923  0.5417798 -0.02923694
# Median   0.1565457  0.6383878       0.4124925  0.7517258  0.27450478
# Mean     0.1333825  0.5196932       0.3510526  0.6678712  0.22207821
# 3rd Qu.  0.4636497  0.7956905       0.6454395  0.8663823  0.50837766
# Max.     0.9568189  0.9899450       0.9808465  0.9951426  0.97086631
write.csv(expr, 'counts_ruvg_all_features.csv', row.names=F, quote=F)
zeros=apply(expr[,-c(1,2)], 1, function(x) length(which(x==0)))
expr=expr[(which(zeros<6)),]
expr=expr[!is.na(expr$Gene),]
expr=aggregate(expr[,-c(1,2,3)], by=list(expr$Gene), FUN=sum)
expr=data.frame(Gene=expr[,1], genes[match(expr[,1], genes$Gene), c(2,3)], expr[,-1])
test_cor(expr)
#          Astrocyte Cholinergic     Choroid Dopaminergic Endothelial  Ependymal  Microglia
# Min.    -0.8872317  -0.8249936 -0.85388354   -0.8758288 -0.86576165 -0.7626769 -0.7799344
# 1st Qu.  0.3392508   0.0262490  0.07818136    0.1013966  0.05877334  0.1402865 -0.1498205
# Median   0.6294498   0.3087639  0.34306015    0.4163428  0.33569567  0.3965144  0.1710808
# Mean     0.5193329   0.2788825  0.28353350    0.3480505  0.28981830  0.3490776  0.1561468
# 3rd Qu.  0.7999434   0.5568977  0.53623054    0.6496634  0.57079584  0.5962381  0.4932113
# Max.     0.9858948   0.9607182  0.93641820    0.9796844  0.97700503  0.9327160  0.9800072
#              Mural     Neuron Oligodendrocyte        OPC    Purkinje
# Min.    -0.8629047 -0.8299676      -0.8613016 -0.6366964 -0.86710090
# 1st Qu. -0.1789355  0.3631321       0.1001923  0.5417798 -0.01512147
# Median   0.1659637  0.6383878       0.4124925  0.7517258  0.29816689
# Mean     0.1374084  0.5196932       0.3510526  0.6678712  0.23430409
# 3rd Qu.  0.4749576  0.7956905       0.6454395  0.8663823  0.52241725
# Max.     0.9568189  0.9899450       0.9808465  0.9951426  0.97086631
write.csv(expr, 'counts_ruvg_expressed_features.csv', row.names=F, quote=F)
expr=expr[expr$type=='protein_coding',]
write.csv(expr, 'counts_ruvg_protein_coding.csv', row.names=F, quote=F)
expr[,-c(1,2,3)]=normalizeQuantiles(expr[,-c(1,2,3,4,5,6)])
test_cor(expr)
#          Astrocyte Cholinergic     Choroid Dopaminergic  Endothelial  Ependymal  Microglia
# Min.    -0.9240909  -0.9294336 -0.90213049   -0.9358621 -0.913529344 -0.8413732 -0.8761605
# 1st Qu. -0.2192502  -0.1615079 -0.15474737   -0.2341216 -0.296813456 -0.1319557 -0.2248822
# Median   0.1342160   0.1337200  0.08394982    0.1459606 -0.003822604  0.1153300  0.1007053
# Mean     0.1038116   0.1274052  0.08267510    0.1131662  0.007727227  0.1118117  0.1048464
# 3rd Qu.  0.4449895   0.4511859  0.32604906    0.4939164  0.302177134  0.3557703  0.4273394
# Max.     0.9715759   0.9687023  0.94939713    0.9759425  0.965846471  0.9365889  0.9908282
#               Mural      Neuron Oligodendrocyte        OPC     Purkinje
# Min.    -0.94091470 -0.91676389     -0.95427164 -0.9542691 -0.902278876
# 1st Qu. -0.26801134 -0.29246250     -0.29926553 -0.2765174 -0.276388236
# Median   0.09465093  0.05598889      0.04925125  0.2010396  0.005479464
# Mean     0.07742318  0.06067621      0.03979461  0.1390255  0.004738940
# 3rd Qu.  0.41077784  0.41351307      0.37380033  0.5893545  0.290709425
# Max.     0.97463066  0.97567666      0.98128494  0.9854445  0.961643359
expr[,-c(1,2,3)]=log2(expr[,-c(1,2,3,4,5,6)]+1)
test_cor(expr)
#          Astrocyte Cholinergic     Choroid Dopaminergic  Endothelial  Ependymal   Microglia
# Min.    -0.9682336  -0.9514949 -0.94271290   -0.9637136 -0.903225956 -0.8706444 -0.90006572
# 1st Qu. -0.2221610  -0.1704805 -0.15609639   -0.2712998 -0.303124928 -0.1093326 -0.25590633
# Median   0.1323329   0.1630221  0.09954864    0.1778037  0.003090642  0.1475579  0.10151599
# Mean     0.1035754   0.1400852  0.08699569    0.1185280  0.005190999  0.1368731  0.09276797
# 3rd Qu.  0.4322817   0.4832233  0.34346311    0.5186143  0.295216972  0.3908997  0.42305773
# Max.     0.9771554   0.9588152  0.93326015    0.9639850  0.950142584  0.9281977  0.98421666
#               Mural      Neuron Oligodendrocyte        OPC     Purkinje
# Min.    -0.96961578 -0.92546545     -0.97273784 -0.9832560 -0.951858408
# 1st Qu. -0.29940356 -0.29466764     -0.30143751 -0.2514677 -0.294539688
# Median   0.09683370  0.05957541      0.03802917  0.2347910 -0.004393351
# Mean     0.06350507  0.05701961      0.03283576  0.1638876  0.001203100
# 3rd Qu.  0.41261903  0.40998771      0.36032110  0.6373170  0.305562739
# Max.     0.96964553  0.97547892      0.97910682  0.9845487  0.894029506
expr=read.csv('counts_ruvg_protein_coding.csv')
exprPC=prcomp(t(expr[,-c(1,2,3,4,5,6)]))
m1=apply(expr[,-c(1,2,3,4,5,6)], 1, function(x) lm(x ~ exprPC$x[,1]+exprPC$x[,2]))
m2=t(data.frame(lapply(m1, function(x) resid(x))))
rownames(m2)=seq(1,nrow(m2))
m2=data.frame(expr[,c(1,2,3)], m2)
test_cor(m2)
#         Astrocyte Cholinergic     Choroid Dopaminergic Endothelial  Ependymal   Microglia
# Min.    -0.9085823 -0.89259421 -0.91111781  -0.91079778 -0.93135273 -0.8131257 -0.88428309
# 1st Qu. -0.1100306  0.07303689  0.06446005   0.02861538 -0.22416337  0.1125072 -0.25655666
# Median   0.1992836  0.35584973  0.31181810   0.38881210  0.12726910  0.3886364  0.08329688
# Mean     0.1450541  0.30925670  0.25569410   0.28440749  0.09575702  0.3211956  0.06964994
# 3rd Qu.  0.4327655  0.60328520  0.51656554   0.63949944  0.42191304  0.5833105  0.40566727
# Max.     0.9570737  0.96168852  0.93456579   0.96803800  0.97415231  0.9311089  0.98922978
#              Mural      Neuron Oligodendrocyte         OPC   Purkinje
# Min.    -0.9658270 -0.91719344     -0.94647187 -0.92177922 -0.8842796
# 1st Qu. -0.1634284 -0.03273419     -0.20877086 -0.30437409 -0.1505503
# Median   0.2389955  0.33933717      0.12389736  0.07422443  0.1009127
# Mean     0.1656373  0.26887517      0.09142659  0.06770101  0.0800408
# 3rd Qu.  0.5209450  0.61581775      0.41047198  0.44908315  0.3229567
# Max.     0.9548142  0.95067062      0.94692557  0.98443886  0.8818593
m1=apply(expr[,-c(1,2,3,4,5,6)], 1, function(x) lm(x ~ exprPC$x[,1]+exprPC$x[,2]))
m2=t(data.frame(lapply(m1, function(x) resid(x))))
rownames(m2)=seq(1,nrow(m2))
m2=data.frame(expr[,c(1,2,3)], m2)
test_cor(m2)
#          Astrocyte Cholinergic     Choroid Dopaminergic Endothelial   Ependymal   Microglia
# Min.    -0.82800125 -0.73155404 -0.63763214  -0.71990420 -0.79055188 -0.71326611 -0.71684974
# 1st Qu. -0.08573238 -0.06816807 -0.09654815  -0.03809689 -0.09832522 -0.08626172 -0.07287826
# Median   0.13590861  0.12928620  0.08856739   0.18126858  0.08820872  0.09117333  0.12497769
# Mean     0.13056512  0.12125605  0.08533541   0.16899540  0.07958547  0.08905282  0.11932365
# 3rd Qu.  0.34938956  0.31310673  0.27321832   0.38148887  0.26691738  0.26640357  0.31252167
# Max.     0.91539375  0.87532810  0.83096202   0.90066688  0.87808016  0.80301195  0.94876681
#               Mural       Neuron Oligodendrocyte         OPC    Purkinje
# Min.    -0.79952890 -0.760131308     -0.75619271 -0.76612944 -0.74929702
# 1st Qu. -0.06362493 -0.004363875     -0.12088893 -0.16713258 -0.13480605
# Median   0.13501474  0.254921996      0.06618108  0.08784230  0.05092864
# Mean     0.13353940  0.226512822      0.06288144  0.08062644  0.05090579
# 3rd Qu.  0.34239367  0.475911354      0.25732740  0.33114569  0.23844125
# Max.     0.90402615  0.925436036      0.82914838  0.89940498  0.81835081
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

setwd('/home/patrick/opentrons_testing_run/06_expression_matrix')
expr=read.csv('counts_ruvg_expressed_features.csv')
expr=read.csv('counts_ruvg_protein_coding.csv')
source('~/code/git/FindModules/FindModules.R')
setwd('~/opentrons_testing_run/07_sample_networks')
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
