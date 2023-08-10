#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("argparse"))
suppressPackageStartupMessages(library('xlsx'))

# take arguments from command line
parser = ArgumentParser(description = "Generate Sample Sheets for sequencing of Lexogen FWD 3' Libraries. Supports both basic i5/i7 6nt indexing and the UDIs.")

parser$add_argument('samplebarcs',
                    type='character',
                    metavar = 'samplebarcs',
                    help="location of the csv file matching sample names to barcode plates and wells")

parser$add_argument('--rc',
                    type='logical',
                    default = 'FALSE',
                    metavar = '',
                    help="reverse complement the i5 barcode, default is '%(default)s'")

parser$add_argument('--investigator',
                    type='character',
                    default = 'anonymous',
                    metavar = '',
                    help="name of the investigator, default is: '%(default)s'")

parser$add_argument('--experiment',
                    type='character',
                    default = 'UCSF experiment',
                    metavar = '',
                    help="experiment name, default is: '%(default)s'")

parser$add_argument('--date',
                    type='character',
                    default = format(Sys.time(), "%x"),
                    metavar = '',
                    help="date of the run, default is today: '%(default)s'")

parser$add_argument('--read1',
                    type = 'integer',
                    default = 66,
                    metavar = '',
                    help="length of read 1, default is: '%(default)s'")

parser$add_argument('--read2',
                    type = 'integer',
                    default = 0,
                    metavar = '',
                    help="length of read 2, default is: '%(default)s'")

parser$add_argument('--i7r',
                    type = 'integer',
                    metavar = '',
                    default = 6,
                    help = "length of the first index read, i7, default is: '%(default)s'")

parser$add_argument('--i5r',
                    type = 'integer',
                    metavar = '',
                    default = 6,
                    help = "length of the second index read, i5, default is: '%(default)s'")

# process arguments and raise error if there are any with no values
args = parser$parse_args()
noArgs = unlist(lapply(args, is.null))
if(any(noArgs)){
    stop(paste0('need arguments for options: ', paste0(names(args)[noArgs], collapse=', ')))
}

# reassign arguments to more intuitive variable names
i7Cycles = args$i7r
i5Cycles = args$i5r
sampleSheet = read.csv(args$samplebarcs)
input=read.csv(args$samplebarcs)
if (any(colnames(input)== 'Sample_Well')){
    typeLex = 'UDI'
} else{
    typeLex = 'standard'
}

# get the relevant barcode excel sheets for the UDI barcodes...
if(typeLex == 'UDI'){
    barcodes = lapply(seq(1,5), function(i) 
        {read.xlsx('../lexogen-barcodes/107UI264V0104_UDI-12-nt-Index-Sequences-for-Illumina_2021-03-26.xlsx', sheetIndex=i)[1:96,seq(1,7)]}
    )
    names(barcodes)=c(paste0('Set_A', seq(1,4)), 'Set_B1')
    outA = lapply(seq_along(input$Sample_Plate), function(i) {
        inputT = input[i,]
        plate=input$Sample_Plate[i]
        barcI = grep(plate, names(barcodes))
        joinedBarc=gsub(' ' , '', apply(barcodes[[barcI]][,c(1,2)], 1, function(x) paste0(x[1], x[2])))
        out=cbind(i, inputT,
            barcodes[[barcI]][match(inputT$Sample_Well,  joinedBarc),-seq(1,3)])
        colnames(out) = c('Sample_ID','Sample_Name','Sample_Plate','Sample_Well','I7_Index_ID','index','I5_Index_ID','index2')
        # process barcodes based on barcode length
        if(i7Cycles != 12){
            out$index = substr(out$index, 1, i7Cycles)
        }
        if(i5Cycles != 12){
            out$index = substr(out$index, 1, i5Cycles)
        }
        return(out)
    })
}

# and the standard 6nt barcodes...
if(typeLex == 'standard'){
    barcodes = lapply(seq(1,3), function(i) 
        {read.xlsx('../lexogen-barcodes/044UI263V0100_i7-and-i5-nt-Index-Sequences-for-Illumina.xlsx', sheetIndex=i)[1:96,seq(1,4)]}
    )
    names(barcodes)=c('i7', 'i5', 'i5rc')
    outA = lapply(seq_along(input$Sample_Name), function(i) {
        inputT = input[i,]
        joinedBarc=gsub(' ' , '', apply(barcodes[[1]][,c(1,2)], 1, function(x) paste0(x[1], x[2])))
        if(args$rc == FALSE){
            out=cbind(i, inputT,
                barcodes[[1]][match(inputT$Sample_Well_i7,  joinedBarc),-seq(1,3)], 
                barcodes[[2]][match(inputT$Sample_Well_i5,  joinedBarc),-seq(1,3)])
            colnames(out) = c('Sample_ID','Sample_Name', 'Sample_Well_i7', 'Sample_Well_i5','index','index2')
        } else if(args$rc == TRUE){
            out=cbind(i, inputT,
                barcodes[[1]][match(inputT$Sample_Well_i7,  joinedBarc),-seq(1,3)], 
                barcodes[[3]][match(inputT$Sample_Well_i5,  joinedBarc),-seq(1,3)])
            colnames(out) = c('Sample_ID','Sample_Name', 'Sample_Well_i7', 'Sample_Well_i5', 'index','index2')
        }
        # process barcodes based on barcode length
        if(i7Cycles == 0){
            out = out[,-c(3,5)]
        }
        if(i5Cycles == 0){
            out = out[,-c(4,6)]
        }
        return(out)
    })
}

# write out the new sample sheet
fileName = paste0('samplesheet_', gsub(' ', '-', args$experiment), '_', gsub('/', '-', format(Sys.time(), "%x")), '.csv')
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
write.table(do.call(rbind, outA), file = fileName, append = TRUE, quote = FALSE, row.names = FALSE, sep = ',')
