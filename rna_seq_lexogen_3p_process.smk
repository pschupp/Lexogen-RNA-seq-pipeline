threads = workflow.cores
import os
cwd = os.getcwd()
import numpy
import pandas
import re
sample=dirpath = '01_basecalled'
sample=[f for f in os.listdir(dirpath) if os.path.isfile(os.path.join(dirpath, f))]
sample=list(filter(lambda x:'.fastq' in x, sample))
sample=[re.sub(".fastq.*", "", x) for x in sample]

rule all:
    input: 
        "09_multiqc/multiqc_report.html"

rule umi_tools:
    input:"01_basecalled/{file_names}.fastq.gz"
    output:"01_basecalled/{file_names}.umi.fastq.gz"
    log:"01_basecalled/{file_names}.log.txt"
    threads: 1
    shell:
        """
        umi_tools extract \\
        --extract-method string \\
        --bc-pattern NNNNNN \\
        -L {log} \\
        --stdin {input} \\
        --stdout {output}
        """

rule fastqc_pre:
    input: "01_basecalled/{file_names}.umi.fastq.gz"
    output: "02_pretrim_fastqc/{file_names}.pretrim.std_out.txt"
    threads: 3
    shell:
        """
        fastqc {input} \\
        --outdir 02_pretrim_fastqc/ \\
        --noextract \\
        --threads {threads} \\
        > {output} \\
        2> 02_pretrim_fastqc/{wildcards.file_names}.pretrim.std_err.txt
        """

rule trimmomatic:
    input:
        real="01_basecalled/{file_names}.fastq.gz",
        dummy="02_pretrim_fastqc/{file_names}.pretrim.std_out.txt"
    output:"03_trimmed/{file_names}.trimmed.fastq"
    params:
        headcrop=4,
        clipfile="/usr/share/trimmomatic/TruSeq3-SE.polyA.fa"
    threads: 1
    shell:
        """
        TrimmomaticSE {input.real} {output} \\
        -trimlog {output}.trim_full_log.txt \\
        -threads {threads} \\
        ILLUMINACLIP:{params.clipfile}:2:30:7 \\
        LEADING:10 \\
        TRAILING:10 \\
        SLIDINGWINDOW:4:15 \\
        MINLEN:30 \\
        HEADCROP:{params.headcrop} \\
        > 03_trimmed/{wildcards.file_names}.std_out.txt \\
        2> 03_trimmed/{wildcards.file_names}.std_err.txt
        """

rule fastqc_post:
    input: "03_trimmed/{file_names}.trimmed.fastq"
    output: "04_posttrim_fastqc/{file_names}.std_out.txt"
    threads: 3
    shell:
        """
        fastqc {input} \\
        --outdir 04_posttrim_fastqc/ \\
        --noextract \\
        --threads {threads} \\
        > {output} \\
        2> 04_posttrim_fastqc/{wildcards.file_names}.std_err.txt
        """

rule star_align:
    input: 
        dummy="04_posttrim_fastqc/{file_names}.std_out.txt",
        files="03_trimmed/{file_names}.trimmed.fastq"
    output: "05_aligned/{file_names}.Aligned.out.bam"
    threads: 15
    params:
        genomeDir="/home/shared/hg_align_db/GRCh38_gencode_primary/star_index_ercc_phix"
    shell:
        """
        /opt/STAR_old/STAR-2.6.1e/bin/Linux_x86_64/STAR --genomeDir {params.genomeDir} \\
        --readFilesIn {input.files} --readFilesCommand cat \\
        --outFileNamePrefix 05_aligned/{wildcards.file_names}. \\
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
        > 05_aligned/{wildcards.file_names}.std_out.txt
        2> 05_aligned/{wildcards.file_names}.std_err.txt
        """
# already done by trimmomatic
#        --clip3pAdapterSeq AAAAAAAA --clip3pAdapterMMp 0.1 \\

rule create_field_OL:
    input: "05_aligned/{file_names}.Aligned.out.bam"
    output: "05_aligned/{file_names}.fixed_tag.bam"
    threads: 1
    shell:
        """
        set +o pipefail;
        samtools view -h -@ {threads} {input} \\
        | sed -E 's/([0-9]{{3,}})(:)([0-9]{{3,}})(:)([0-9]{{3,}})(_)([A-Z]{{6}})(.*)$/\\1:\\3:\\5:\\7\\8\tOL:Z:\\7/g' \\
        | samtools view -O BAM -b -@ {threads} -o {output} - 
        """

rule sort_query:
    input: "05_aligned/{file_names}.fixed_tag.bam"
    output: "06_sorted/{file_names}.bam"
    threads: 5
    shell:
        """
        PicardCommandLine SortSam \\
        I={input} \\
        O={output} \\
        SORT_ORDER=queryname \\
        > 06_sorted/{wildcards.file_names}.sort.std_out.txt \\
        2> 06_sorted/{wildcards.file_names}.sort.std_err.txt
        """

rule deduplicate:
    input: "06_sorted/{file_names}.bam"
    output: "07_deduplicated/{file_names}.deduplicated.bam"
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
        > 07_deduplicated/{wildcards.file_names}.markdup.std_out.txt \\
        2> 07_deduplicated/{wildcards.file_names}.markdup.std_err.txt
        """

rule feature_counts:
    input: "07_deduplicated/{file_names}.deduplicated.bam"
    output: "08_feature_count/{file_names}.deduplicated.bam.featureCounts.bam"
    threads: 15
    params:
        gtfFile="/home/shared/hg_align_db/GRCh38_gencode_primary/gencode.v38.primary_assembly.annotation.ercc.phix.gtf"
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
    input: "08_feature_count/{file_names}.deduplicated.bam.featureCounts.bam"
    output: "09_sort_index/{file_names}.sort.bam"
    shell:
        """
        samtools sort \\
        {input} \\
        -o {output}
        > 09_sort_index/{wildcards.file_names}.sort.std_out.txt \\
        2> 09_sort_index/{wildcards.file_names}.sort.std_err.txt
        """

rule bam_index:
    input: "09_sort_index/{file_names}.sort.bam"
    output: "09_sort_index/{file_names}.sort.bam.bai"
    shell:
        """
        samtools index \\
        {input} \\
        > 09_sort_index/{wildcards.file_names}.index.std_out.txt \\
        2> 09_sort_index/{wildcards.file_names}.index.std_err.txt
        """

rule genebody_cov:
    input: 
        real="09_sort_index/{file_names}.sort.bam",
        dummy="09_sort_index/{file_names}.sort.bam.bai"
    output: "09_sort_index/{file_names}.geneBodyCoverage.curves.pdf"
    shell:
        """
        geneBody_coverage.py \\
        -i {input.real} \\
        -r ~/@patrick/parsebio/hg38_GENCODE.v38.10000.bed \\
        -o 11_sort_index/{wildcards.file_names} \\
        -f pdf
        """

rule multiqc:
    input: 
        inDir=".",
#        bodycov=expand("09_sort_index/{file_names}.sort.bam.bai", file_names=sample),
        extra=expand("09_sort_index/{file_names}.geneBodyCoverage.curves.pdf", file_names=sample)
#         dummy=expand("08_feature_count/{file_names}.deduplicated.bam.featureCounts.bam", file_names=sample)
    output:
        outFile="09_multiqc/multiqc_report.html"
    shell:
        """
        multiqc {input.inDir} \\
        -o 09_multiqc; \\
        multiqc 02_pretrim_fastqc -o 02_pretrim_fastqc
        """
