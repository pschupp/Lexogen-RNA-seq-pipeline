# How to use this repo
Each ```snakemake``` script can be run the usual way with Snakemake in its own virtual environment via:
```{bash}
# for dry run
#   -n : dry run
#   -r : reason for each excecuted rule
#   -s : snakemake file
snakemake -n -r -s [name of file.smk] 
# create DAG diagram
snakemake --dag -s [name of file.smk]| dot -Tsvg > dag.svg
# create DAG diagram
snakemake --rulegraph -s [name of file.smk]| dot -Tsvg > rules.svg
# for real run
snakemake --cores 15 -s [name of file.smk]
```
# Current issues
1. breaking on samples which are not actually there. need to supply the correct sample sheet with all samples included.
2. will break on barcodes that only have a few reads. so if it is poorly sequenced, will drop out...unlikely however, most likely just carry over from barcode via sequencing errors
3. Alignment appears to be all to all to phiX when this is demonstrably untrue
4. Are fasta and gtf file set up correctly. Some appear to be lacking the full assortment of phix and ercc
5. alignments might not be ideal because of read through, should limit maxium read length... but how to limit to cutoff when hitting a stretch of certain homopolymers?
6. Should trim out spacer sequence, already done or what is it?
# Fixes to pipeline
1. DONE Check on gtf and fasta files to make sure they are accurate
2. DONE If above isn't true, should regenerate STAR indeces with current version of START
3. DONE Fix sample sheet and include only the correct barcodes that were included
4. look up how to trim reads after poly A sequence
5. determine if spacer sequence (what is it?) has been cinluded as port of the sequenceds to edit out
6. optional: write quick tool in bash to extract read 1 and read 2 and list all barcodes and their frequency. probably a tool already exists to do this. Find it and run it.
# outputs to be analyzed
1. What percent is phiX? # excluded
2. What percent is mRNA v. rRNA? DONE
3. What percent are aligning internally? TODO, see 5
4. Sequence length distribution post-trimming DONE
5. Histogram of read alignment locations relatively to entire transcript (canonical)
# responses to problems
1. FASTA with combined was never generated. STAR index did use combined fasta. Alignments do exist to fasta and ercc, but not phix. 
REAL PROBLEM is whether there exists a gtf file with the gtf of ercc. i think the settings for feeature counts are such that it shouldn't care!
2. No PhiX is present. I conclude  that these reads must be eliminated prior to inclusion in fastq, i.e. filtered out by basecalling
3. There is a ERCC read but it is not included in feature counting...why?
4. DONE : add ERCC and phiX (just in case) to gtf file
5. DONE Sample sheet creation is involved and uses multiple spreadsheets. Hard to duplicate and there are some inevitable inconsistentices. For future work, it makes sense write a script to automate this process.
6. from trimmomatic add poly A sequence AND check to see if we can add adapter sequence, poly A trimming seems to have worked well. only possible caveat is that there are some short reads remaining (<30nt). These are likely to drop out as multimappers though.
7. cant remove multimappers via trimmomatic
8. alignments seem good
9. rerun featureCounts
10. all seems googd
