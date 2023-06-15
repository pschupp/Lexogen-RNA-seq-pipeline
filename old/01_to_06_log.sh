# in R prepare sample sheet for all barcocdes
library('stringr')
barc=fread('barcodes.csv', select=seq(1,8), nrow=96)
out=data.frame(Sample_ID=paste0('Oldham-PS-3616_', str_pad(seq(1,96), width=2, side='left', pad='0')), Sample_Name=paste0(barc[[1]],barc[[2]]), Description=barc[[3]], index=barc[[5]],I7_Index_ID=barc[[4]], index2=barc[[7]],I5_Index_ID=barc[[6]], Sample_Project=rep('Oldham-PS-3616',96))
write.table(out, file='SampleSheet.csv', row.names=F, quote=F, sep=',')

# add the following text
[Header]
Local Run Manager Analysis Id,76076
Experiment Name,Oldham-PS-3616
Date,2021-09-17
Module,GenerateFASTQ - 2.0.1
Workflow,GenerateFASTQ
Library Prep Kit,Custom
Chemistry,Amplicon
 
[Reads]
68
 
[Settings]
Read1EndWithCycle,68 
Read1StartFromCycle,10 
Read1UMILength,6 
Read1UMIStartFromCycle,1 
TrimUMI,1 

[Data]

# Settings explained
Read1EndWithCycle,68 # The last cycle to use for Read 1.
Read1StartFromCycle,10 # The first cycle to use for Read 1.
Read1UMILength,6 # The length of the UMI used for Read 1.
Read1UMIStartFromCycle,1 # The first cycle to use for UMI in Read 1. The cycle index is absolute and not affected by the Read1StartFromCycle setting. The software supports UMIs only at the beginning or end of reads. Unless paired with Read1UMILength, the software ignores this setting.
TrimUMI,1 # TrimUMI 0—False (default). 1—True. The software trims the UMI bases from Read 1 and Read 2.

mkdir /opt/illumina/bcl2fastq_source # source code is here
cd /opt/illumina/bcl2fastq_source
wget 'https://files.softwaredownloads.illumina.com/82660c4c-f46c-4743-8566-2437755e4329/bcl2fastq2-v2-20-0-tar.zip?Expires=1632771982&Policy=eyJTdGF0ZW1lbnQiOlt7IlJlc291cmNlIjoiaHR0cHM6Ly9maWxlcy5zb2Z0d2FyZWRvd25sb2Fkcy5pbGx1bWluYS5jb20vODI2NjBjNGMtZjQ2Yy00NzQzLTg1NjYtMjQzNzc1NWU0MzI5L2JjbDJmYXN0cTItdjItMjAtMC10YXIuemlwIiwiQ29uZGl0aW9uIjp7IkRhdGVMZXNzVGhhbiI6eyJBV1M6RXBvY2hUaW1lIjoxNjMyNzcxOTgyfX19XX0_&Signature=fS~4SIQjG~rQUh~Muj5XbMPkJWtur6ds-jWqrxCzy6G4EU-s-2UhZtVH7AF7yhd8FkHrsuhy7D25A3xaaX43SKi2j~0vMyOokznHObfLk2yJR8B4P8W6ttFxVoyc~ETfKVBcFKVglqMT9x0JeDIWVMZx7Zwgm9K-pMe0Z15H46tWHqP6nMQVb7jReZUMoyYX0dfrorXiROztA6nqA-CK2wADBl3cNtxdGvs1To-eNeMz2X2TFd~Uvu-Vj1NU-yHtL~5uPndEsFGtVrDSF4NZPdUQElV8Cg5tJxqUjkF-1KiDpMqUP8zWjzFoFj~CyecS7jmWWCK9Ltt~t8o-HfJYsg__&Key-Pair-Id=APKAJO3UYWPXK4A26FPQ' # temporary url, need to login
mkdir /opt/illumina/bcl2fastq # compiled binaries will be put here, according prefix
export C_INCLUDE_PATH=/usr/include/x86_64-linux-gnu # include essential libraries for compliation
./src/configure --prefix=../bcl2fastq
make
make install
# --minimum-trimmed-read-length arg (=35)         minimum read length after adapter trimming
# --with-failed-reads                             include non-PF clusters
# --no-lane-splitting                             do not split fastq files by lane.
# --barcode-mismatches arg (=1)                   number of allowed mismatches per index. Multiple, comma delimited, entries allowed. Each entry is applied to the corresponding index; last entry applies to all remaining indices. Accepted values: 0, 1, 2.
# --loading-threads Number of threads to load BCL data
# --processing-threads Number of threads to process demultiplexing data.
# --writing-threads Number of threads to write FASTQ data. This number must be lower than number of samples.

/opt/illumina/bcl2fastq/bin/bcl2fastq --runfolder ~/bdata/opentrons_testing_run/00_raw_output_from_NextSeq500 --output-dir ~/opentrons_testing_run/01_basecalled_fastq --minimum-trimmed-read-length 30 --with-failed-reads --sample-sheet ~/bdata/opentrons_testing_run/00_raw_output_from_NextSeq500/SampleSheet.csv --barcode-mismatches 2 --no-lane-splitting  --loading-threads 5 --processing-threads 15 --writing-threads 5

# say to use truseq
for ea in ~/opentrons_testing_run/01_basecalled_fastq/Oldham-PS-3616/*/*fastq.gz; do
	ea2=$( echo "$ea"  | sed 's/fastq/trim\.fastq/')
	ea2=$( echo "$ea2"  | sed 's/.*\///')
	TrimmomaticSE -threads 20 -trimlog ~/opentrons_testing_run/02_trimmed_reads/$ea2.log.trim  $ea $ea2 ILLUMINACLIP:/usr/share/trimmomatic/TruSeq3-SE.fa:2:30:10 LEADING:10 TRAILING:10 SLIDINGWINDOW:4:15 MINLEN:30 2>> ~/opentrons_testing_run/02_trimmed_reads/$ea2.log
done

#should consider how trimming effects output...maybe should not do...

STAR --runThreadN 15 --runMode genomeGenerate --genomeDir /home/shared/hg_align_db/GRCh38_gencode_primary --genomeFastaFiles GRCh38.primary_assembly.genome.fa --sjdbGTFfile /home/shared/hg_align_db/GRCh38_gencode_primary/gencode.v38.primary_assembly.annotation.gtf --sjdbOverhang 58

for ea in ~/opentrons_testing_run/02_trimmed_reads/*fastq.gz; do 
	ea2=$( echo "$ea"  | sed 's/fastq\.gz/align.bam/')
	ea2=$( echo "$ea2"  | sed 's/.*\///')
	STAR  --runThreadN 15 --genomeDir /home/shared/hg_align_db/GRCh38_gencode_primary --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate --readFilesIn $ea --outFilterMismatchNmax 999 --alignIntronMax 1000000 --outFilterMultimapNmax 20 --outSAMmultNmax 1 --outFilterScoreMinOverLread 0 --outFilterMatchNminOverLread 0 --genomeLoad LoadAndKeep --limitBAMsortRAM 322122382273 --outFileNamePrefix ~/opentrons_testing_run/03_alligned_bam/$ea2
done

for ea in ~/opentrons_testing_run/03_alligned_bam/*_001.align.bamAligned.sortedByCoord.out.bam; do
	ea2=$( echo "$ea"  | sed 's/001.align.bamAligned.sortedByCoord.out.bam/aligned.feature.tsv/')
	ea2=$( echo "$ea2"  | sed 's/.*\///')
	featureCounts -a  /home/shared/hg_align_db/GRCh38_gencode_primary/gencode.v38.primary_assembly.annotation.gtf -o ~/opentrons_testing_run/04_feature_count/$ea2 -R SAM $ea -T 15 -g gene_id -M --primary -O -t gene
done

# Saturation analysis

snakemake -n --forcerun $(snakemake --list-input-changes)

configfile: "config.yaml"

ea=E3_S21_R1_001.align.bamAligned.sortedByCoord.out.bam
samtools view -h $ea | sed -E 's/([0-9])(:)([A-Z]{6})(.*)$/\1:\3\4\tFL:Z:\3/g' | samtools view -O BAM -b -o $ea.fixed.bam -

for prob in {0.01,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9}; do
	sem -j 15 "picard-tools DownsampleSam \
	I=E3_S21_R1_001.align.bamAligned.sortedByCoord.out.bam.fixed.bam \
	O=E3_down_$prob.bam \
	STRATEGY=ConstantMemory \
	RANDOM_SEED=7 \
	PROBABILITY=$prob \
	METRICS_FILE=E3_metrics_$prob.txt"
done
sem --wait

for ea in E3_down_*fixed.bam; do 
	sem -j 15 "picard-tools MarkDuplicates INPUT=$ea \
	OUTPUT=$ea.marked.dup.bam \
	TAGGING_POLICY=All \
	BARCODE_TAG=FL \
	MOLECULAR_IDENTIFIER_TAG=FL \
	SORTING_COLLECTION_SIZE_RATIO=.8 \
	MAX_FILE_HANDLES=1000 \
	METRICS_FILE=$ea.marked.dup.metrics.txt"
done
sem --wait
