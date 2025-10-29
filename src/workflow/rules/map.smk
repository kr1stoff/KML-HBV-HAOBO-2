rule bwa_mem:
    input:
        reads=rules.fastp.output.trimmed,
        ref=rules.copy_ref.output,
        idx=rules.bwa_index.output,
    output:
        "align/{sample}.bam",
    benchmark:
        ".log/align/{sample}.bwa_mem.bm"
    log:
        ".log/align/{sample}.bwa_mem.log",
    conda:
        config["conda"]["bwa"]
    params:
        extra=r"-M -Y -R '@RG\tID:{sample}\tSM:{sample}\tPL:ILLUMINA'",
    threads: config["threads"]["high"]
    shell:
        """
        (bwa mem -t {threads} {params.extra} {input.ref} {input.reads} \
            | samtools view -@ {threads} -hbS - \
            | samtools sort -@ {threads} -o {output} -) 2> {log}
        """


rule samtools_index:
    input:
        rules.bwa_mem.output,
    output:
        "align/{sample}.bam.bai",
    benchmark:
        ".log/align/{sample}.samtools_index.bm"
    log:
        ".log/align/{sample}.samtools_index.log",
    conda:
        config["conda"]["samtools"]
    shell:
        "samtools index {input} {output} 2> {log}"


rule samtools_stats:
    input:
        bam=rules.bwa_mem.output,
        bed=rules.make_target_bed.output,
    output:
        "align/{sample}.bam.target.stat",
    benchmark:
        ".log/align/{sample}.samtools_stats.bm"
    log:
        ".log/align/{sample}.samtools_stats.log",
    conda:
        config["conda"]["samtools"]
    threads: config["threads"]["low"]
    shell:
        "samtools stats --threads {threads} --target-regions {input.bed} {input.bam} > {output} 2> {log}"


rule samtools_stats_all:
    input:
        rules.bwa_mem.output,
    output:
        "align/{sample}.bam.stat",
    benchmark:
        ".log/align/{sample}.samtools_stats_all.bm"
    log:
        ".log/align/{sample}.samtools_stats_all.log",
    conda:
        config["conda"]["samtools"]
    threads: config["threads"]["low"]
    shell:
        "samtools stats --threads {threads} {input} > {output} 2> {log}"


rule samtools_depth:
    input:
        bam=rules.bwa_mem.output,
        bed=rules.make_target_bed.output,
    output:
        "align/{sample}.bam.target.depth",
    benchmark:
        ".log/align/{sample}.samtools_depth.bm"
    log:
        ".log/align/{sample}.samtools_depth.log",
    conda:
        config["conda"]["samtools"]
    shell:
        "samtools depth -a -b {input.bed} {input.bam} > {output} 2> {log}"


rule bam_stats:
    input:
        rules.samtools_stats_all.output,
        rules.samtools_stats.output,
        rules.samtools_depth.output,
    output:
        "align/{sample}.stats.csv",
    benchmark:
        ".log/align/{sample}.bam_stats.bm"
    log:
        ".log/align/{sample}.bam_stats.log",
    conda:
        config["conda"]["python"]
    script:
        "../scripts/bam_stats.py"


rule bam_stats_summary:
    input:
        expand("align/{sample}.stats.csv", sample=config["samples"]),
    output:
        "align/bam_summary.tsv",
    benchmark:
        ".log/align/bam_stats_summary.bm"
    log:
        ".log/align/bam_stats_summary.log",
    conda:
        config["conda"]["python"]
    script:
        "../scripts/bam_stats_summary.py"


rule samtools_bedcov:
    input:
        rules.bwa_mem.output,
        rules.make_target_bed.output,
        rules.samtools_index.output,
    output:
        "align/{sample}.bam.target.bedcov",
    benchmark:
        ".log/align/{sample}.samtools_bedcov.bm"
    log:
        ".log/align/{sample}.samtools_bedcov.log",
    conda:
        config["conda"]["samtools"]
    shell:
        "samtools bedcov -c {input[1]} {input[0]} > {output} 2> {log}"
