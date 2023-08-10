#!/usr/bin/env Rscript
suppressPackageStartupMessages(library("argparse"))
suppressPackageStartupMessages(library('xlsx'))

# goal: automate production of sample shwet for sequencing and demultiplexing with minimal user input

# input:
# 1. sample-associated-barcode.csv - csv of sample name, sample barcode plate, and sample barcode well
# 2. DONE - lexogen barcode sheet 

# output: 
# 1. Sample sheet with following columns:
#       Sample_ID	Sample_Name	Sample_Plate	Sample_Well	I7_Index_ID	index	I5_Index_ID	index2	Sample_Project	Description
# 2. Optional: full sample sheet with appropriate header, etc., but this is fairly trivial and be a later TODO
########################
parser = ArgumentParser(description = "Generate Sample Sheets for sequencing of Lexogen FWD 3' Libraries. Supports both basic i5/i7 6nt indexing and the UDIs.")

parser$add_argument('--type', '-t', 
                    type='character',
                    metavar = '',
                    help="barcode type from Lexogen, one of 'UDI' or 'standard'")

parser$add_argument('--rc', '-c', 
                    type='logical',
                    default = 'FALSE',
                    metavar = '',
                    help="reverse complement the i5 barcode")

parser$add_argument('--samplebarcs', '-s',
                    type='character',
                    metavar = '',
                    help="location of csv containing columns: sample name, sample barcode plate, sample barcode well")

parser$add_argument('--workingdir', '-w',
                    type='character',
                    metavar = '',
                    help="current working directory")

parser$add_argument('--investigator', '-i',
                    type='character',
                    default = 'anonymous',
                    metavar = '',
                    help="name of the investigator, default is: '%(default)s'")

parser$add_argument('--experiment', '-e',
                    type='character',
                    default = 'UCSF experiment',
                    metavar = '',
                    help="experiment name, default is: '%(default)s'")

parser$add_argument('--date', '-d',
                    type='character',
                    default = format(Sys.time(), "%x"),
                    metavar = '',
                    help="date of the run, default is today: '%(default)s'")

parser$add_argument('--read1', '-r',
                    type = 'integer',
                    default = 66,
                    metavar = '',
                    help="length of read 1, default is: '%(default)s'")

parser$add_argument('--read2', '-u',
                    type = 'integer',
                    default = 0,
                    metavar = '',
                    help="length of read 2, default is: '%(default)s'")

parser$add_argument('--i7r', '-v',
                    type = 'integer',
                    metavar = '',
                    help = "length of the first index read, i7")

parser$add_argument('--i5r', '-x',
                    type = 'integer',
                    metavar = '',
                    help = "length of the second index read, i5")

args <- parser$parse_args()
WD <- args$working_dir
noArgs = unlist(lapply(args, is.null))
if(any(noArgs)){
    stop(paste0('need arguments for options: ', paste0(names(args)[noArgs], collapse=', ')))
}
if(args$type == 'UDI'){
    barcodes = lapply(seq(1,5), function(i) 
        {read.xlsx('../lexogen-barcodes/107UI264V0104_UDI-12-nt-Index-Sequences-for-Illumina_2021-03-26.xlsx', sheetIndex=i)[1:96,seq(1,7)]}
    )
    names(barcodes)=c(paste0('Set_A', seq(1,4)), 'Set_B1')
}
if(args$type == 'standard'){
    barcodes = lapply(seq(1,3), function(i) 
        {read.xlsx('../lexogen-barcodes/044UI263V0100_i7-and-i5-nt-Index-Sequences-for-Illumina.xlsx', sheetIndex=i)[1:96,seq(1,4)]}
    )
    names(barcodes)=c('i7', 'i5', 'i5rc')
}

sampleSheet = read.csv(args$samplebarcs)
i7Cycles = args$i7r
i5Cycles = args$i5r
fileName = paste0('samplesheet_', args$experiment, '_', gsub('/', '-', format(Sys.time(), "%x")), '.csv')
cat('[Header]\n', file = fileName)
cat('FileFormatVersion,2\n', file = fileName, append = TRUE)
cat(paste0('Investigator Name,', args$investigator, '\n'), file = fileName, append = TRUE)
cat(paste0('RunName,', args$experiment, '\n'), file = fileName, append = TRUE)
cat(paste0('Date,', args$date, '\n'), file = fileName, append = TRUE)
cat('InstrumentPlatform,NextSeq1k2k\nInstrumentType,NextSeq2000\n\n', file = fileName, append = TRUE)
cat('[Reads]\n', file = fileName, append = TRUE)
cat(paste0('Read1Cycles,', args$read1, '\n'), file = fileName, append = TRUE)
cat(paste0('Read2Cycles,', args$read2, '\n'), file = fileName, append = TRUE)
cat(paste0('Index1Cycles,', i7Cycles, '\n'), file = fileName, append = TRUE)
cat(paste0('Index2Cycles,', i5Cycles, '\n\n'), file = fileName, append = TRUE)
cat('[BCLConvert_Data]\n', file = fileName, append = TRUE)

input=read.csv(args$samplebarcs)

outA = lapply(seq_along(input$Sample_Plate), function(i) {
    inputT = input[i,]
    plate=input$Sample_Plate[i]
    barcI = grep(plate, names(barcodes))
    joinedBarc=gsub(' ' , '', apply(barcodes[[barcI]][,c(1,2)], 1, function(x) paste0(x[1], x[2])))
    out=cbind(i, inputT, 
        barcodes[[barcI]][match(inputT$Sample_Well,  joinedBarc),-seq(1,3)]) 
    colnames(out) = c('Sample_ID','Sample_Name','Sample_Plate','Sample_Well','I7_Index_ID','index','I5_Index_ID','index2')
    return(out)
})

write.table(do.call(rbind, outA), file = fileName, append = TRUE, quote = FALSE, row.names = FALSE, sep = ',')

#outA = do.call(rbind, outA)
#system('cp sample_top.txt new_sample_sheet.txt')
#lapply(seq_len(nrow(outA)), function(x) write(paste0(outA[x,], collapse=','), file='new_sample_sheet.txt', append = TRUE))
