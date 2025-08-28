rule freebayes:
    input:
        alns=rules.bwa_mem.output,
        idxs=rules.samtools_index.output,
        ref=rules.bedtools_maskfasta.output,
        fai=rules.samtools_faidx.output,
    output:
        vcf="variant/{sample}.vcf",
    log:
        ".log/variant/{sample}.freebayes.log",
    benchmark:
        ".log/variant/{sample}.freebayes.bm"
    params:
        extra=(
            "--ploidy 1 --min-repeat-size 10 --read-indel-limit 15 --use-best-n-alleles 4 "
            "--theta 0.001 --haplotype-length 0 --min-alternate-fraction 0.001 --min-base-quality 30 "
            "--min-coverage 20 --min-alternate-count 2 --min-mapping-quality 30 --max-complex-gap 1 --trim-complex-tail "
            "--targets " + rules.bedtools_sort.output[0]
        ),
    threads: config["threads"]["low"]
    # resources:
    #     mem_mb=1024,
    conda:
        # ! 需要安装 freebayes 和 bcftools
        config["conda"]["basic2"]
    wrapper:
        f"file:{workflow.basedir}/wrappers/freebayes"


rule vcf2tab:
    input:
        rules.freebayes.output,
    output:
        "variant/{sample}.tsv",
    log:
        ".log/variant/{sample}.vcf2tab.log",
    benchmark:
        ".log/variant/{sample}.vcf2tab.bm"
    conda:
        config["conda"]["python"]
    script:
        "../scripts/vcf2tab.py"


rule csvtk_csv2xlsx:
    input:
        expand("variant/{sample}.tsv", sample=config["samples"]),
    output:
        "variant/vars.xlsx",
    log:
        ".log/variant/csvtk_csv2xlsx.log",
    benchmark:
        ".log/variant/csvtk_csv2xlsx.bm"
    conda:
        config["conda"]["basic"]
    params:
        extra="--tabs --format-numbers",
    shell:
        "csvtk csv2xlsx {params.extra} {input} -o {output} 2> {log}"
