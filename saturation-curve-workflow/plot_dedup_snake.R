library('data.table', quietly=T)
library('openblasctl')
openblas_set_num_threads(15)

fileL=unlist(snakemake@input[grep('*metrics.txt', snakemake@input)])
input=data.frame(matrix(nrow=9, ncol=length(fileL)))
for(i in seq(1, length(fileL))){
    input[,i]=t(fread(fileL[i], skip='LIBRARY')[,-1])
    if(i==length(fileL)){
        rownames(input)=colnames(fread(fileL[i], skip='LIBRARY')[,-1])
    }
}
colnames(input)=gsub('.*\\/|\\.deduplicated\\.metrics.txt', '', fileL)
input=t(input[c(1,5,8), order(colnames(input))])

file2=unlist(snakemake@input[grep('*featureCounts', snakemake@input)])
feats=data.frame(matrix(ncol=2, nrow=length(file2)))
feats[,1]=file2
for(i in seq(1, length(file2))){
    temp=data.frame(fread(file2[i]))
    temp=temp[-grep('^ERCC-', temp$Geneid),]
    temp=temp[which(temp[,7]>10),]
    feats[i,2]=nrow(temp)
}
feats[,1]=gsub('.*\\/|\\.deduplicated.*', '', feats[,1])
colnames(feats)=c('sample', 'Features')

library('ggplot2')
plotIn=data.frame(samples=gsub('_.*', '', rownames(input)),  input)
plotIn=data.frame(plotIn, UNIQUE_READS=plotIn$UNPAIRED_READS_EXAMINED-plotIn$UNPAIRED_READ_DUPLICATES, Features=feats$Features)
colnames(plotIn)=c('samples', 'totalReads', 'duplicateReads', 'percentDuplication', 'uniqueReads', 'features')
write.csv(plotIn, file=snakemake@output[[1]], quote=F)

pdf(snakemake@output[[2]])
print(ggplot(plotIn, aes(x=totalReads, y=uniqueReads, color=samples)) + 
    geom_point() + 
    scale_y_continuous(labels = function(x) format(x, scientific = TRUE)) +
    scale_x_continuous(labels = function(x) format(x, scientific = TRUE)) +
    geom_abline(mapping=aes(slope=1, intercept=lm(plotIn$uniqueReads~plotIn$totalReads)$coefficients[1]), color='red')+
    geom_abline(mapping=aes(slope=lm(plotIn$uniqueReads~plotIn$totalReads)$coefficients[2], intercept=lm(plotIn$uniqueReads~plotIn$totalReads)$coefficients[1]), color='black')+
    labs(title='Unique reads v. Total reads', subtitle='Red line with slope of 1, black line of best fit', x='Total reads', y='Unique reads', color='Samples')+
    theme_linedraw()
)

print(ggplot(plotIn, aes(x=totalReads, y=features, color=samples)) + 
    geom_point() + 
    scale_y_continuous(labels = function(x) format(x, scientific = TRUE)) +
    scale_x_continuous(labels = function(x) format(x, scientific = TRUE)) +
    geom_abline(mapping=aes(slope=lm(plotIn$features~plotIn$totalReads)$coefficients[2], intercept=lm(plotIn$features~plotIn$totalReads)$coefficients[1]))+
    labs(title='Features v. Total reads', subtitle='Black line of best fit', x='Total reads', y='Features', color='Samples')+
    theme_linedraw()
)

print(ggplot(plotIn, aes(x=totalReads, y=percentDuplication, color=samples)) +
    geom_point() +
    labs(title='Total reads v. percent duplication', x='Total reads', y='Percent duplication', color='Samples')+
    theme_linedraw()
)

plotIn[,-1]=log2(plotIn[,-1])
print(ggplot(plotIn, aes(x=totalReads, y=uniqueReads, color=samples)) +
    geom_point() +
    scale_y_continuous(labels = function(x) format(x, scientific = TRUE)) +
    scale_x_continuous(labels = function(x) format(x, scientific = TRUE)) +
    geom_abline(mapping=aes(slope=lm(plotIn$uniqueReads~plotIn$totalReads)$coefficients[2], intercept=lm(plotIn$uniqueReads~plotIn$totalReads)$coefficients[1]))+
    labs(title='Log(2) Unique reads v. Total reads', subtitle='Line of best fit', x='Total reads (log2)', y='Unique reads (log2)', color='Samples')+
    theme_linedraw()
)

garbage=dev.off()

