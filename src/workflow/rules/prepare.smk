rule bwa_index:


rule bedtools_sort:
    input:
        in_file=config["database"]["bed"],
    output:
        ".temp/target.sorted.bed",
    log:
        ".log/prepare/bedtools_sort.log",
    benchmark:
        ".log/prepare/bedtools_sort.bm"
    conda:
        config["conda"]["bedtools"]
    shell:
        "bedtools sort -i {input.in_file} > {output} 2> {log}"


# * 不同分型后续考虑
rule tantan:
    input:
        config["database"]["ref"]
    output:
        ".temp/tantan.fasta"
    log:
        ".log/prepare/tantan.log"
    benchmark:
        ".log/prepare/tantan.bm"
    conda:
        config["conda"]["tantan"]
    shell:
        "tantan {input} > {output} 2> {log}"


rule extract_lcr_bed:
    input:
        rules.tantan.output
    output:
        ".temp/lcr_regions.bed"
    log:
        ".log/prepare/extract_lcr_bed.log",
    benchmark:
        ".log/prepare/extract_lcr_bed.bm"
    conda:
        config["conda"]["python"]
    script:
        "../scripts/tantan_lcr_extractor.py"
