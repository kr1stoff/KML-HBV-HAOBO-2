# * 不同分型后续考虑
rule copy_ref:
    input:
        config["database"]["ref"]
    output:
        "prepare/ref.fasta"
    log:
        ".log/ref/copy_ref.log",
    benchmark:
        ".log/ref/copy_ref.bm"
    conda:
        config["conda"]["bwa"]
    shell:
        "cp {input} {output} 2> {log}"


rule bwa_index:
    input:
        rules.copy_ref.output
    output:
        multiext("prepare/ref.fasta", ".amb", ".ann", ".bwt", ".pac", ".sa")
    log:
        ".log/ref/bwa_index.log",
    benchmark:
        ".log/ref/bwa_index.bm"
    conda:
        config["conda"]["bwa"]
    shell:
        "bwa index {input} 2> {log}"


rule samtools_faidx:
    input:
        rules.copy_ref.output
    output:
        "prepare/ref.fasta.fai"
    log:
        ".log/ref/samtools_faidx.log",
    benchmark:
        ".log/ref/samtools_faidx.bm"
    conda:
        config["conda"]["samtools"]
    shell:
        "samtools faidx {input} 2> {log}"


rule tantan:
    input:
        config["database"]["ref"]
    output:
        "prepare/tantan.fasta"
    log:
        ".log/ref/tantan.log"
    benchmark:
        ".log/ref/tantan.bm"
    conda:
        config["conda"]["tantan"]
    shell:
        "tantan {input} > {output} 2> {log}"


rule extract_lcr_bed:
    input:
        rules.tantan.output
    output:
        "prepare/lcr_regions.bed"
    log:
        ".log/ref/extract_lcr_bed.log",
    benchmark:
        ".log/ref/extract_lcr_bed.bm"
    conda:
        config["conda"]["python"]
    script:
        "../scripts/tantan_lcr_extractor.py"
