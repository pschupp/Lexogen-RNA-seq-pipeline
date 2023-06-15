IDS, = glob_wildcards("{sample}.align.bamAligned.sortedByCoord.out.bam")
threads = workflow.cores
import numpy
probs = numpy.round(numpy.arange(0.1, 1.1, 0.1), 2)

rule all:
    input: expand("{sample}.fixed_tag.prob_{pr}_downsampled.sort.deduplicated.bam", sample=IDS, pr=probs)

rule tag_adjust:
    input: "{sample}.align.bamAligned.sortedByCoord.out.bam"
    output: "{sample}.fixed_tag.bam"
    shell:
        """
        set +o pipefail;
        samtools view -h -@ {threads} {input} \\
        | sed -E 's/([0-9])(:)([A-Z]{{6}})(.*)$/\\1:\\3\\4\tRX:Z:\\3/g' \\
        | samtools view -O BAM -b -@ {threads} -o {output} - 
        """

import numpy
rule downsample_bam:
    input: "{sample}.fixed_tag.bam"
    output: "{sample}.fixed_tag.prob_{pr}_downsampled.bam"
    shell:
        """
        PicardCommandLine DownsampleSam \\
        I={input} \\
        O={output} \\
        STRATEGY=ConstantMemory \\
        RANDOM_SEED=7 \\
        PROBABILITY={wildcards.pr} \\
        METRICS_FILE={output}.metric 2> {wildcards.sample}.downsample_{wildcards.pr}_log.txt
        """

rule sort_bam:
    input: "{sample}.fixed_tag.prob_{pr}_downsampled.bam"
    output: "{sample}.fixed_tag.prob_{pr}_downsampled.sort.bam"
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
    input: "{sample}.fixed_tag.prob_{pr}_downsampled.sort.bam"
    output: "{sample}.fixed_tag.prob_{pr}_downsampled.sort.deduplicated.bam"
    shell:
        """
        PicardCommandLine MarkDuplicates \\
        I={input} \\
        O={output} \\
        TAGGING_POLICY=All \\
        ASSUME_SORT_ORDER=coordinate \\
        BARCODE_TAG=RX \\
        MOLECULAR_IDENTIFIER_TAG=RX \\
        SORTING_COLLECTION_SIZE_RATIO=.8 \\
        MAX_FILE_HANDLES=1000 \\
        METRICS_FILE={output}.metrics 2> {wildcards.sample}.downsample_{wildcards.pr}_log.txt
        """

IDS2, = glob_wildcards("{sampleM}.deduplicated.bam.metrics")

rule plot_duplicates:
    input: expand("{sampleM}.deduplicated.bam.metrics", sampleM=IDS2)
    output: "saturation_curve.csv",
            "saturation_curve.pdf"
    script: "plot_dedup_snake.R"
