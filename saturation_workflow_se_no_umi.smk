IDS, = glob_wildcards("{sample}.bam")
threads = workflow.cores
import numpy
probs = numpy.round(numpy.arange(0.1, 1.1, 0.1), 2)

rule all:
    input: "saturation_curve.pdf"
    #    input: expand("{sample}.prob_{pr}_downsampled.sort.deduplicated.bam", sample=IDS, pr=probs)

# rule tag_adjust:
#     input: "{sample}.align.bamAligned.sortedByCoord.out.bam"
#     output: "{sample}.fixed_tag.bam"
#     shell:
#         """
#         set +o pipefail;
#         samtools view -h -@ {threads} {input} \\
#         | sed -E 's/([0-9])(:)([A-Z]{{6}})(.*)$/\\1:\\3\\4\tRX:Z:\\3/g' \\
#         | samtools view -O BAM -b -@ {threads} -o {output} - 
#         """
#             ), sample=IDS, pr=lambda ldcards: *probs),


import numpy
rule downsample_bam:
    input: "{sample}.bam"
    output: bam="{sample}.prob_{pr}_downsampled.bam",
            metrics="{sample}.prob_{pr}_downsampled.metrics",
            logs="{sample}.prob_{pr}_downsampled.log.txt"
    shell:
        """
        PicardCommandLine DownsampleSam \\
        I={input} \\
        O={output.bam} \\
        STRATEGY=ConstantMemory \\
        RANDOM_SEED=7 \\
        PROBABILITY={wildcards.pr} \\
        METRICS_FILE={output.metrics} 2> {output.logs}
        """

rule sort_bam:
    input: "{sample}.prob_{pr}_downsampled.bam"
    output: "{sample}.prob_{pr}_downsampled.sort.bam"
    shell:
        """
        samtools sort \\
        -o {output} \\
        -O BAM \\
        -l 9 \\
        -m 100G \\
        -@ {threads} \\
        {input}
        """

rule mark_duplicates:
    input: "{sample}.prob_{pr}_downsampled.sort.bam"
    output: bam="{sample}.prob_{pr}_downsampled.sort.deduplicated.bam",
            metrics="{sample}.prob_{pr}_downsampled.dedup.metrics",
            logs="{sample}.prob_{pr}_downsampled.dedup.log.txt"
    shell:
        """
        PicardCommandLine MarkDuplicates \\
        I={input} \\
        O={output.bam} \\
        TAGGING_POLICY=All \\
        ASSUME_SORT_ORDER=coordinate \\
        BARCODE_TAG=RX \\
        MOLECULAR_IDENTIFIER_TAG=RX \\
        SORTING_COLLECTION_SIZE_RATIO=.8 \\
        MAX_FILE_HANDLES=1000 \\
        METRICS_FILE= {output.metrics} 2> {output.logs}
        """

#IDS2, = glob_wildcards("{sampleM}.deduplicated.bam.metrics")

rule plot_duplicates:
    input: expand("{sample}.prob_{pr}_downsampled.dedup.metrics", sample=IDS, pr=probs)
    output: "saturation_curve.csv",
            "saturation_curve.pdf"
    script: "plot_dedup_snake.R"

# junkyard
# IDS2, = glob_wildcards("{sampleM}.deduplicated.bam.metrics")

# rule plot_duplicates:
#     input: "05_downsample/03_features/{file_names}.downsampled_{pr}.deduplicated.features.bam"
#     output: "saturation_curve.csv",
#             "saturation_curve.pdf"
#     script: "plot_dedup_snake.R"

#     rule all:
#    input: expand("{sample}.fixed_tag.prob_{pr}_downsampled.sort.deduplicated.bam", sample=IDS, pr=probs)



