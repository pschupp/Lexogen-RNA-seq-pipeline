# curl sftp://fastq.ucsf.edu/SSD/211118_K00153_0941_BHMJCHBBXY_PSMO --user hiseq_user:b1gdata -o nov_2021_run
threads = workflow.cores
import os
cwd = os.getcwd()
import numpy
import pandas
sample=pandas.read_table('00_raw_seq_data/SampleSheet.csv', sep=',', header=20)
IDS=list(sample["Sample_ID"])
NAMES=list(sample["Sample_Name"])
NAMES_T=list(zip(NAMES , range(1,len(NAMES)+1)))
NAMES_F=[]
for names_i in NAMES_T:
    NAMES_F.append(names_i[0]+'_S'+ str(names_i[1])+'_R1_001')

IDS=IDS[1]
NAMES_F=NAMES_F[1]

rule all:
    input: expand("05_feature_count/{file_name}.deduplicated.bam.featureCounts.bam", file_name=NAMES_F)

rule bcl_2_fastq:
    input: "00_raw_seq_data/RTAComplete.txt"
    output: "01_basecalled"
    shell:
        """
        /opt/illumina/bcl2fastq/bin/bcl2fastq \\
        --runfolder 00_raw_seq_data \\
        --output-dir {output} \\
        --sample-sheet 00_raw_seq_data/SampleSheet.csv \\
        --minimum-trimmed-read-length 30 --with-failed-reads  --barcode-mismatches 1 --no-lane-splitting \\
        --loading-threads 5 --processing-threads 15 --writing-threads 5 \\
        --tiles s_[2-5] 2> 00_basecalling.log.txt
        """

rule trimmomatic:
    input: expand("01_basecalled/PSOM/{name}/{file_name}.fastq.gz", zip, name=IDS, file_name=NAMES_F)
    output: expand("02_trimmed/{file_name}.trimmed.fastq", file_name=NAMES_F)
    threads: 15
    run:
        for i in range(len(input)):
            shell("TrimmomaticSE {input_i} {output_i} -trimlog {output_i}.trim_full_log.txt -threads {threads} ILLUMINACLIP:/usr/share/trimmomatic/TruSeq3-SE.fa:2:30:10 LEADING:10 TRAILING:10 SLIDINGWINDOW:4:15 MINLEN:30 2>> {output_i}.std_out_log.txt".format(input_i=input[i], output_i=output[i], threads=threads))

rule fastqc:
    input: "02_trimmed/{file_names}.trimmed.fastq"
    output: "02_trimmed/{file_names}.fastqc_log.txt"
    threads: 15
    shell:
        """
        fastqc {input} \\
        -o 02_trimmed/ \\
        --noextract \\
        --threads {threads} > {output}
        """

rule star_align:
    input: "02_trimmed/{file_names}.fastqc_log.txt" 
    output: "03_aligned/{file_names}.Aligned.out.bam"
    threads: 15
    params:
        genomeDir="/home/shared/hg_align_db/GRCm39_gencode_primary/star_index_ercc_phix"
    shell:
        """
        STAR --genomeDir {params.genomeDir} \\
        --readFilesIn 02_trimmed/{wildcards.file_names}.trimmed.fastq --readFilesCommand cat \\
        --outFileNamePrefix 03_aligned/{wildcards.file_names}. \\
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
        --runThreadN {threads}  2> {output}.log.txt
        """

rule create_field_OL:
    input: "03_aligned/{file_names}.Aligned.out.bam"
    output: "03_aligned/{file_names}.fixed_tag.bam"
    threads: 15
    shell:
        """
        set +o pipefail;
        samtools view -h -@ {threads} {input} \\
        | sed -E 's/([0-9]{{3,}})(:)([0-9]{{3,}})(:)([0-9]{{3,}})(:)([A-Z]{{6}})(.*)$/\\1:\\3:\\5:\\7\\8\tOL:Z:\\7/g' \\
        | samtools view -O BAM -b -@ {threads} -o {output} - 
        """

rule sort_query:
    input: "03_aligned/{file_names}.fixed_tag.bam"
    output: "04_deduplicated/{file_names}.sort_query.bam"
    threads: 15
    shell:
        """
        PicardCommandLine SortSam \\
        I={input} \\
        O={output} \\
        SORT_ORDER=queryname 
        """

rule deduplicate:
    input: "04_deduplicated/{file_names}.sort_query.bam"
    output: "04_deduplicated/{file_names}.deduplicated.bam"
    threads: 15
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
        METRICS_FILE={output}.metrics 2> {output}_log.txt
        """

rule feature_counts:
    input: "04_deduplicated/{file_names}.deduplicated.bam"
    output: "05_feature_count/{file_names}.deduplicated.bam.featureCounts.bam"
    threads: 15
    params:
        gtfFile="/home/shared/hg_align_db/GRCm39_gencode_primary/gencode.vM27.primary_assembly.annotation.ercc.phix.gtf"
    shell:
        """
        featureCounts -a {params.gtfFile} \\
        {input} -o 05_feature_count/{wildcards.file_names}.count_matrix.tsv \\
        -R BAM \\
        -g gene_id \\
        -M -O \\
        -t gene \\
        --fraction \\
        -T {threads}
        """

rule multimap_4mer_equal_dedup_umi:
    input: "04_feature_count/{file_names}.bam"
    output: 
    threads: 15
