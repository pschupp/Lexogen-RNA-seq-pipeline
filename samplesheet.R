# goal: automate production of sample shwet for sequencing and demultiplexing with minimal user input

# input:
# 1. csv of sample name, sample barcode plate, and sample barcode well
# 2. lexogen barcode sheet DONE 

# output: 
# 1. Sample sheet with following columns:
#       Sample_ID	Sample_Name	Sample_Plate	Sample_Well	I7_Index_ID	index	I5_Index_ID	index2	Sample_Project	Description
# 2. Optional: full sample sheet with appropriate header, etc., but this is fairly trivial and be a later TODO
########################
# USER DEFINED PORTION HERE, MUST NAME FILE "samples_associated_barcodes.csv"
WD = '~/@patrick/SF11949_SF11055_run/00_raw_seq_data'
########################

library('xlsx')
setwd('~/code/git/workflows')
barcodes = lapply(seq(1,5), function(i) 
    {read.xlsx('~/code/git/workflows/lexogen_barcodes/107UI264V0104_UDI-12-nt-Index-Sequences-for-Illumina_2021-03-26.xlsx', sheetIndex=i)[1:96,seq(1,7)]}
)
names(barcodes)=c(paste0('Set_A', seq(1,4)), 'Set_B1')
input=read.csv(paste0(WD, '/samples_associated_barcodes.csv'))
outA = lapply(unique(input$index_plate), function(plate) {
    inputT = input[which(input$index_plate == plate),]
    barcI = grep(plate, names(barcodes))
    joinedBarc=gsub(' ' , '', apply(barcodes[[barcI]][,c(1,2)], 1, function(x) paste0(x[1], x[2])))
    out=cbind(inputT[,1], inputT, 
        barcodes[[barcI]][match(inputT$index_well,  joinedBarc),-seq(1,3)], 
        gsub('(.*)_.*', '\\1', inputT$name), 
        gsub('(.*)_.*', '\\1', inputT$name))
    colnames(out) = c('Sample_ID','Sample_Name','Sample_Plate','Sample_Well','I7_Index_ID','index','I5_Index_ID','index2','Sample_Project','Description')
    return(out)
})
outA = do.call(rbind, outA)
write.csv(outA, file=paste0(WD, '/sample_sheet_bottom.csv'), row.names=F, quote=F)
