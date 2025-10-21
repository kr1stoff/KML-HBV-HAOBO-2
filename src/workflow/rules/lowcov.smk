rule bedtools_genomecov:
    input:
        rules.bedtools_intersect.output.bam,
    output:
        temp("lowcov/{sample}.genomecov"),
    log:
        ".log/lowcov/{sample}.bedtools_genomecov.log",
    benchmark:
        ".log/lowcov/{sample}.bedtools_genomecov.bm"
    params:
        "-bga",
    conda:
        config["conda"]["bedtools"]
    shell:
        "bedtools genomecov {params} -ibam {input} > {output} 2> {log}"


rule bedtools_genomecov_filter:
    input:
        rules.bedtools_genomecov.output,
    output:
        temp("lowcov/{sample}.genomecov.filtered"),
    log:
        ".log/lowcov/{sample}.bedtools_genomecov_filter.log",
    benchmark:
        ".log/lowcov/{sample}.bedtools_genomecov_filter.bm"
    params:
        depth_threshold=1000,
    shell:
        "awk '$4<{params.depth_threshold}' {input} > {output} 2> {log}"


rule bedtools_merge:
    input:
        rules.bedtools_genomecov_filter.output,
    output:
        temp("lowcov/{sample}.lowcovmask.raw.bed"),
    log:
        ".log/lowcov/{sample}.bedtools_merge.log",
    benchmark:
        ".log/lowcov/{sample}.bedtools_merge.bm"
    conda:
        config["conda"]["bedtools"]
    shell:
        "bedtools merge -i {input} > {output} 2> {log}"


rule bedtools_merge_filter:
    input:
        rules.bedtools_merge.output,
    output:
        "lowcov/{sample}.lowcovmask.bed",
    log:
        ".log/lowcov/{sample}.bedtools_merge_filter.log",
    benchmark:
        ".log/lowcov/{sample}.bedtools_merge_filter.bm"
    params:
        length_threshold=20,
    shell:
        "awk '$3-$2>{params.length_threshold}' {input} > {output} 2> {log}"


rule bedtools_maskfasta:
    input:
        fa=config["database"]["ref"],
        bed=rules.bedtools_merge_filter.output,
    output:
        "lowcov/{sample}.masked.fa",
    log:
        ".log/lowcov/{sample}.bedtools_maskfasta.log",
    benchmark:
        ".log/lowcov/{sample}.bedtools_maskfasta.bm"
    conda:
        config["conda"]["bedtools"]
    shell:
        "bedtools maskfasta -fi {input.fa} -bed {input.bed} -fo {output} 2> {log}"


rule lowcov_samtools_faidx:
    input:
        rules.bedtools_maskfasta.output,
    output:
        "lowcov/{sample}.masked.fa.fai",
    log:
        ".log/lowcov/{sample}.lowcov_samtools_faidx.log",
    benchmark:
        ".log/lowcov/{sample}.lowcov_samtools_faidx.bm"
    conda:
        config["conda"]["samtools"]
    shell:
        "samtools faidx {input} > {output} 2> {log}"
