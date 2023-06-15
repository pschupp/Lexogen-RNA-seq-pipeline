rule tag_adjust:
    input: 
        "{sample}.bam"
    output:
        "{sample}.fixed_tag.bam"
    params:
        config["parameter"]["tag_adjust"]
    shell:
        "samtools view -h {sample}.bam | sed -E {params} | samtools view -O BAM -b -o {sample}.bam -"
