rule all_qc_summary:
    input:
        rules.fq_stats_summary.output,
        rules.bam_stats_summary.output,
    output:
        "upload/panel-qc-summary.tsv",
        "upload/panel-qc-summary.xlsx",
    benchmark:
        ".log/summary/all_qc_summary.bm"
    log:
        ".log/summary/all_qc_summary.log",
    conda:
        config["conda"]["python"]
    script:
        "../scripts/all_qc_summary.py"


rule csvtk_csv2xlsx:
    input:
        expand("variant/{sample}.vcf_test.tsv", sample=config["samples"]),
    output:
        "upload/all-vars-sheets.xlsx",
    log:
        ".log/summary/csvtk_csv2xlsx.log",
    benchmark:
        ".log/summary/csvtk_csv2xlsx.bm"
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
        "upload/all-vars-summary.tsv",
        "upload/all-vars-summary.xlsx",
    log:
        ".log/summary/all_vars_summary.log",
    benchmark:
        ".log/summary/all_vars_summary.bm"
    params:
        # [浩博过滤] 目标区域 1000 个reads
        core_depth_cutoff=1000,
        # [浩博过滤] 非目标区域 20000 个reads
        other_depth_cutoff=20000,
    conda:
        config["conda"]["python"]
    script:
        "../scripts/all_vars_summary.py"
