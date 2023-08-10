#!/usr/bin/env Rscript
suppressPackageStartupMessages(library("argparse"))

# goal: automate production of sample shwet for sequencing and demultiplexing with minimal user input

# input:
# 1. sample-associated-barcode.csv - csv of sample name, sample barcode plate, and sample barcode well
# 2. DONE - lexogen barcode sheet 

# output: 
# 1. Sample sheet with following columns:
#       Sample_ID	Sample_Name	Sample_Plate	Sample_Well	I7_Index_ID	index	I5_Index_ID	index2	Sample_Project	Description
# 2. Optional: full sample sheet with appropriate header, etc., but this is fairly trivial and be a later TODO
########################
parser <- ArgumentParser(description="Generate Sample Sheets for sequencing of Lexogen FWD 3' Libraries. Supports both basic i5/i7 6nt indexing and the UDIs.")
parser$add_argument("-v", "--verbose", action="store_true", default=TRUE, help="Print extra output [default]")
parser$add_argument("--working-dir", type='character',  help="current working directory")
parser$add_argument("--mean", default=0, type="double", help="Mean if generator == \"rnorm\" [default %(default)s]")
args <- parser$parse_args()
print(args)
print('wd here')
WD <- args$working_dir
print(WD)
########################

library('xlsx')
barcodes = lapply(seq(1,5), function(i) 
    {read.xlsx('../lexogen-barcodes/107UI264V0104_UDI-12-nt-Index-Sequences-for-Illumina_2021-03-26.xlsx', sheetIndex=i)[1:96,seq(1,7)]}
)
names(barcodes)=c(paste0('Set_A', seq(1,4)), 'Set_B1')
input=read.csv(paste0(WD, '/sample-associated-barcodes.csv'))
outA = lapply(seq_along(input$Sample_Plate), function(i) {
    inputT = input[i,]
    plate=input$Sample_Plate[i]
    barcI = grep(plate, names(barcodes))
    joinedBarc=gsub(' ' , '', apply(barcodes[[barcI]][,c(1,2)], 1, function(x) paste0(x[1], x[2])))
    out=cbind(i, inputT, 
        barcodes[[barcI]][match(inputT$Sample_Well,  joinedBarc),-seq(1,3)], 
        gsub('(.*)_.*', '\\1', inputT$Sample_Name), 
        gsub('(.*)_.*', '\\1', inputT$Sample_Name))
    colnames(out) = c('Sample_ID','Sample_Name','Sample_Plate','Sample_Well','I7_Index_ID','index','I5_Index_ID','index2','Sample_Project','Description')
    return(out)
})
outA = do.call(rbind, outA)
system('cp sample_top.txt new_sample_sheet.txt')
lapply(seq_len(nrow(outA)), function(x) write(paste0(outA[x,], collapse=','), file='new_sample_sheet.txt', append = TRUE))
