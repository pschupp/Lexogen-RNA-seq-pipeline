threads = workflow.cores
import numpy
probs = list(numpy.round(numpy.arange(0.1, 1.1, 0.1), 2))[0:2]

rule all:
    input: "10_downsample/saturation_curve.csv",
           "10_downsample/saturation_curve.pdf"

ID, = glob_wildcards("04_deduplicated/{sample}.sort_query.bam")
ID=ID[0:3]

rule downsample_bam:
    input: 
            "config['star-dir']/{sample}.sort_query.bam"
    output: 
            bam=    "10_downsample/01_downsample/{sample}.downsampled_{pr}.bam",
            metrics="10_downsample/01_downsample/{sample}.downsampled_{pr}.metrics.txt"
    log: 
            "10_downsample/01_downsample/{sample}.downsampled_{pr}.log.txt"
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


rule feature_counts:
    input: 
            "10_downsample/01_downsample/{sample}.downsampled_{pr}.deduplicated.bam"
    output: 
            bam=    "10_downsample/02_features/{sample}.downsampled_{pr}.deduplicated.bam.featureCounts",
            summary="10_downsample/02_features/{sample}.downsampled_{pr}.deduplicated.bam.featureCounts.summary"
    log: 
            "10_downsample/02_features/{sample}.downsampled_{pr}.deduplicated.features.log.txt"
    threads: 15
    params:
        gtfFile=config['gtf-file']
    shell:
        """
        featureCounts \\
        {input} \\
        -a {params.gtfFile} \\
        -o 10_downsample/02_features/{wildcards.sample}.downsampled_{wildcards.pr}.deduplicated.bam.featureCounts \\
        -R CORE\\
        -g gene_id \\
        -M -O \\
        -t gene \\
        --fraction \\
        -T {threads} \\
        2> {log}
        """

rule plot_duplicates:
    input: 
           expand("10_downsample/02_features/{sample}.downsampled_{pr}.deduplicated.bam.featureCounts", sample=ID, pr=probs)
    output: "10_downsample/saturation_curve.csv",
            "10_downsample/saturation_curve.pdf"
    script: "plot_dedup_snake.R"
    #script: "R_snakemake.R"
