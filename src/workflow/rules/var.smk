rule freebayes:
    input:
        alns=rules.bwa_mem.output,
        idxs=rules.samtools_index.output,
        ref=rules.bedtools_maskfasta.output,
        fai=rules.samtools_faidx.output,
        targets=rules.bedtools_sort.output,
    output:
        vcf="variant/{sample}.vcf",
    log:
        ".log/variant/{sample}.freebayes.log",
    benchmark:
        ".log/variant/{sample}.freebayes.bm"
    # * 用 thread 控制 freebayes 的并行数量，见效内存压力
    threads:
        config["custom"]["freebayes_threads"]
    params:
        (
            "--ploidy 1 --min-repeat-size 10 --read-indel-limit 15 --use-best-n-alleles 4 "
            "--theta 0.001 --haplotype-length 0 --min-alternate-fraction 0.001 --min-base-quality 30 "
            "--min-coverage 20 --min-alternate-count 2 --min-mapping-quality 30 --max-complex-gap 1 --trim-complex-tail "
        ),
    conda:
        config["conda"]["freebayes"]
    shell:
        "freebayes {params} --targets {input.targets} --fasta-reference {input.ref} {input.alns} > {output.vcf} 2> {log}"


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
        "variant/all-vars-sheets.xlsx",
    log:
        ".log/variant/csvtk_csv2xlsx.log",
    benchmark:
        ".log/variant/csvtk_csv2xlsx.bm"
    conda:
        config["conda"]["csvtk"]
    params:
        extra="--tabs --format-numbers",
    shell:
        "csvtk csv2xlsx {params.extra} {input} -o {output} 2> {log}"


rule all_vars_summary:
    input:
        rules.csvtk_csv2xlsx.output,
        config["database"]["known_sites"],
    output:
        "variant/all-vars-summary.tsv",
        "variant/all-vars-summary.xlsx",
    log:
        ".log/variant/all_vars_summary.log",
    benchmark:
        ".log/variant/all_vars_summary.bm"
    params:
        # [浩博过滤] 目标区域 1000 个reads
        core_depth_cutoff=1000,
        # [浩博过滤] 非目标区域 20000 个reads
        other_depth_cutoff=20000,
    conda:
        config["conda"]["python"]
    script:
        "../scripts/all_vars_summary.py"
