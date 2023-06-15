threads = workflow.cores
import numpy
probs = list(numpy.round(numpy.arange(0.1, 1.1, 0.1), 2))[0:2]

rule all:
    input: "05_downsample/saturation_curve.csv",
           "05_downsample/saturation_curve.pdf"

ID, = glob_wildcards("04_deduplicated/{sample}.sort_query.bam")
ID=ID[0:3]

rule downsample_bam:
    input: 
            "04_deduplicated/{sample}.sort_query.bam"
    output: 
            bam=    "05_downsample/01_downsample/{sample}.downsampled_{pr}.bam",
            metrics="05_downsample/01_downsample/{sample}.downsampled_{pr}.metrics.txt"
    log: 
            "05_downsample/01_downsample/{sample}.downsampled_{pr}.log.txt"
    threads: 1
    shell:
        """
        PicardCommandLine DownsampleSam \\
        I={input} \\
        O={output.bam} \\
        STRATEGY=ConstantMemory \\
        RANDOM_SEED=7 \\
        PROBABILITY={wildcards.pr} \\
        METRICS_FILE={output.metrics} \\
        2> {log}
        """

rule deduplicate:
    input: 
            "05_downsample/01_downsample/{sample}.downsampled_{pr}.bam"
    output: 
            bam=    "05_downsample/02_deduplicate/{sample}.downsampled_{pr}.deduplicated.bam",
            metrics="05_downsample/02_deduplicate/{sample}.downsampled_{pr}.deduplicated.metrics.txt"
    log: 
            "05_downsample/02_deduplicate/{sample}.downsampled_{pr}.deduplicated.log.txt"
    threads: 5
    shell:
        """
        PicardCommandLine MarkDuplicates \\
        I={input} \\
        O={output.bam}\\
        TAGGING_POLICY=All \\
        ASSUME_SORT_ORDER=queryname \\
        BARCODE_TAG=OL \\
        MOLECULAR_IDENTIFIER_TAG=OL \\
        REMOVE_DUPLICATES=TRUE \\
        READ_NAME_REGEX=null \\
        SORTING_COLLECTION_SIZE_RATIO=.8 \\
        MAX_FILE_HANDLES=1000 \\
        METRICS_FILE={output.metrics} \\
        2> {log}
        """

rule feature_counts:
    input: 
            "05_downsample/02_deduplicate/{sample}.downsampled_{pr}.deduplicated.bam"
    output: 
            bam=    "05_downsample/03_features/{sample}.downsampled_{pr}.deduplicated.bam.featureCounts",
            summary="05_downsample/03_features/{sample}.downsampled_{pr}.deduplicated.bam.featureCounts.summary"
    log: 
            "05_downsample/03_features/{sample}.downsampled_{pr}.deduplicated.features.log.txt"
    threads: 15
    params:
        gtfFile="/home/shared/hg_align_db/GRCm39_gencode_primary/gencode.vM27.primary_assembly.annotation.ercc.phix.gtf"
    shell:
        """
        featureCounts \\
        {input} \\
        -a {params.gtfFile} \\
        -o 05_downsample/03_features/{wildcards.sample}.downsampled_{wildcards.pr}.deduplicated.bam.featureCounts \\
        -R CORE\\
        -g gene_id \\
        -M -O \\
        -t gene \\
        --fraction \\
        -T {threads} \\
        2> {log}
        """

rule plot_duplicates:
    input: expand("05_downsample/02_deduplicate/{sample}.downsampled_{pr}.deduplicated.metrics.txt", sample=ID, pr=probs),
           expand("05_downsample/03_features/{sample}.downsampled_{pr}.deduplicated.bam.featureCounts", sample=ID, pr=probs)
    output: "05_downsample/saturation_curve.csv",
            "05_downsample/saturation_curve.pdf"
    script: "plot_dedup_snake.R"
    #script: "R_snakemake.R"
