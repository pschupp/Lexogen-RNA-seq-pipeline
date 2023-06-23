threads = workflow.cores
import os
import numpy
import pandas
import re
cwd = os.getcwd()
sample=pandas.read_table('00_raw_seq_data/SampleSheet.csv', header=None, sep=',')
sample=sample[0][(int(numpy.where(sample[0]=='Sample_ID')[0])+1):]
rule all:
    input: 
        "09_multiqc/multiqc_report.html"
            
rule bcl_2_fastq:
    output:"00_raw_seq_data/00_basecalling.log.txt"
    shell:
        """
        /opt/illumina/bcl2fastq/bin/bcl2fastq --runfolder=00_raw_seq_data --output-dir=01_basecalled --sample-sheet=00_raw_seq_data/SampleSheet.csv --minimum-trimmed-read-length=30 --with-failed-reads --barcode-mismatches=1 --no-lane-splitting --loading-threads=5 --processing-threads=15 --writing-threads=5 2> 00_raw_seq_data/00_basecalling.log.txt;
        mv 01_basecalled/sf*/*fastq.gz 01_basecalled;
        mkdir 01_basecalled/Undetermined;
        mv 01_basecalled/Undetermined*fastq.gz 01_basecalled/Undetermined;
        rename 's/_S[0-9]+//' 01_basecalled/*fastq.gz
        """

rule umi_tools:
    input:"00_raw_seq_data/00_basecalling.log.txt"
    output:"01_basecalled/{file_names}.umi.fastq.gz"
    log:"01_basecalled/{file_names}.log.txt"
    threads: 1
    shell:
        """
        umi_tools extract \\
        --extract-method=string \\
        --bc-pattern=NNNNNN \\
        -L {log} \\
        --stdin {input} \\
        --stdout {output}
        """

rule trimmomatic:
    input:"01_basecalled/{file_names}.umi.fastq.gz"
    output:"02_trimmed/{file_names}.trimmed.fastq"
    params:
        headcrop=4,
        clipfile="/usr/share/trimmomatic/TruSeq3-SE.polyA.fa"
    threads: 1
    shell:
        """
        TrimmomaticSE {input} {output} \\
        -trimlog {output}.trim_full_log.txt \\
        -threads {threads} \\
        ILLUMINACLIP:{params.clipfile}:2:30:7 \\
        LEADING:10 \\
        TRAILING:10 \\
        SLIDINGWINDOW:4:15 \\
        MINLEN:30 \\
        HEADCROP:{params.headcrop} \\
        > 02_trimmed/{wildcards.file_names}.std_out.txt \\
        2> 02_trimmed/{wildcards.file_names}.std_err.txt
        """

rule fastqc_post:
    input: "02_trimmed/{file_names}.trimmed.fastq"
    output: "03_posttrim_fastqc/{file_names}.std_out.txt"
    threads: 3
    shell:
        """
        fastqc {input} \\
        --outdir 03_posttrim_fastqc/ \\
        --noextract \\
        --threads {threads} \\
        > {output} \\
        2> 03_posttrim_fastqc/{wildcards.file_names}.std_err.txt
        """

rule star_align:
    input: 
        dummy="03_posttrim_fastqc/{file_names}.std_out.txt",
        files="03_trimmed/{file_names}.trimmed.fastq"
    output: "04_aligned/{file_names}.Aligned.out.bam"
    threads: 15
    params:
        genomeDir=config['genomeDir']
    shell:
        """
        /opt/STAR_old/STAR-2.6.1e/bin/Linux_x86_64/STAR --genomeDir config['genomeDir'] \\
        --readFilesIn {input.files} --readFilesCommand cat \\
        --outFileNamePrefix 04_aligned/{wildcards.file_names}. \\
        --outSAMtype BAM Unsorted \\
        --outFilterMultimapNmax 4 \\
        --outFilterType BySJout \\
        --outSAMmultNmax 1  \\
        --alignSJoverhangMin 5 \\
        --alignSJDBoverhangMin 1 \\
        --outFilterMismatchNmax 999 \\
        --outFilterMismatchNoverReadLmax 0.3 \\
        --outFilterScoreMinOverLread 0.3 \\
        --outFilterMatchNminOverLread 0.3 \\
        --alignIntronMin 20 \\
        --alignIntronMax 1000000 \\
        --genomeLoad LoadAndKeep --limitBAMsortRAM 322122382273  \\
        --runThreadN {threads} \\
        > 04_aligned/{wildcards.file_names}.std_out.txt
        2> 04_aligned/{wildcards.file_names}.std_err.txt
        """
# already done by trimmomatic
#        --clip3pAdapterSeq AAAAAAAA --clip3pAdapterMMp 0.1 \\

rule create_field_OL:
    input: "04_aligned/{file_names}.Aligned.out.bam"
    output: "04_aligned/{file_names}.fixed_tag.bam"
    threads: 1
    shell:
        """
        set +o pipefail;
        samtools view -h -@ {threads} {input} \\
        | sed -E 's/([0-9]{{3,}})(:)([0-9]{{3,}})(:)([0-9]{{3,}})(_)([A-Z]{{6}})(.*)$/\\1:\\3:\\5:\\7\\8\tOL:Z:\\7/g' \\
        | samtools view -O BAM -b -@ {threads} -o {output} - 
        """

rule sort_query:
    input: "04_aligned/{file_names}.fixed_tag.bam"
    output: "05_sorted/{file_names}.bam"
    threads: 5
    shell:
        """
        PicardCommandLine SortSam \\
        I={input} \\
        O={output} \\
        SORT_ORDER=queryname \\
        > 05_sorted/{wildcards.file_names}.sort.std_out.txt \\
        2> 05_sorted/{wildcards.file_names}.sort.std_err.txt
        """

rule deduplicate:
    input: "05_sorted/{file_names}.bam"
    output: "06_deduplicated/{file_names}.deduplicated.bam"
    threads: 5
    shell:
        """
        PicardCommandLine MarkDuplicates \\
        I={input} \\
        O={output} \\
        TAGGING_POLICY=All \\
        ASSUME_SORT_ORDER=queryname \\
        BARCODE_TAG=OL \\
        MOLECULAR_IDENTIFIER_TAG=OL \\
        REMOVE_DUPLICATES=TRUE \\
        READ_NAME_REGEX=null \\
        SORTING_COLLECTION_SIZE_RATIO=.8 \\
        MAX_FILE_HANDLES=1000 \\
        METRICS_FILE={output}.metrics \\
        > 06_deduplicated/{wildcards.file_names}.markdup.std_out.txt \\
        2> 06_deduplicated/{wildcards.file_names}.markdup.std_err.txt
        """

rule feature_counts:
    input: "06_deduplicated/{file_names}.deduplicated.bam"
    output: "07_feature_count/{file_names}.deduplicated.bam.featureCounts.bam"
    threads: 15
    params:
        gtfFile=config['gtfFile']
    shell:
        """
        featureCounts \\
        {input} \\
        -a {params.gtfFile} \\
        -o 08_feature_count/{wildcards.file_names}.count_matrix.tsv \\
        -R BAM \\
        -g gene_id \\
        -M -O \\
        -t gene \\
        --fraction \\
        -T {threads} \\
        > 08_feature_count/{wildcards.file_names}.std_out.txt \\
        2> 08_feature_count/{wildcards.file_names}.std_err.txt
        """
# options:
# -a gtf file input
# -o output file count matrix
# -R input file format
# -g ninth column of gtf file. the attribute type by which to group features. look to gtf file for correct naming. could also be transcript_id or exon_id
# -t third column of gtf file. indicates what type of features should be considered. When using -g gene_id should use gene. If using transcript_id or exon_id, should use processed transcript.
# -M multi-mapping reads will counted. Uses NH tag for multimappers.
# -O assign reads to all overlaping features
# --fraction Assign fractional counts to features

rule bam_sort:
    input: "07_feature_count/{file_names}.deduplicated.bam.featureCounts.bam"
    output: "08_sort_index/{file_names}.sort.bam"
    shell:
        """
        samtools sort \\
        {input} \\
        -o {output}
        > 08_sort_index/{wildcards.file_names}.sort.std_out.txt \\
        2> 08_sort_index/{wildcards.file_names}.sort.std_err.txt
        """

rule bam_index:
    input: "08_sort_index/{file_names}.sort.bam"
    output: "08_sort_index/{file_names}.sort.bam.bai"
    shell:
        """
        samtools index \\
        {input} \\
        > 08_sort_index/{wildcards.file_names}.index.std_out.txt \\
        2> 08_sort_index/{wildcards.file_names}.index.std_err.txt
        """

rule genebody_cov:
    input: 
        real="08_sort_index/{file_names}.sort.bam",
        dummy="08_sort_index/{file_names}.sort.bam.bai"
    output: "08_sort_index/{file_names}.geneBodyCoverage.curves.pdf"
    params:
        bedFile=config['bedFile']
    shell:
        """
        geneBody_coverage.py \\
        -i {input.real} \\
        -r {params.bedFile} \\
        -o 11_sort_index/{wildcards.file_names} \\
        -f pdf
        """

rule multiqc:
    input: 
        inDir=".",
#        bodycov=expand("09_sort_index/{file_names}.sort.bam.bai", file_names=sample),
        extra=expand("08_sort_index/{file_names}.geneBodyCoverage.curves.pdf", file_names=sample)
#         dummy=expand("08_feature_count/{file_names}.deduplicated.bam.featureCounts.bam", file_names=sample)
    output:
        outFile="09_multiqc/multiqc_report.html"
    shell:
        """
        multiqc {input.inDir} \\
        -o 09_multiqc; \\
        multiqc 02_pretrim_fastqc -o 02_pretrim_fastqc
        """
