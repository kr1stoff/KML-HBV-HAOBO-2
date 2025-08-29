rule all_summary:
    input:
        rules.fq_stats_summary.output,
        rules.bam_stats_summary.output,
    output:
        "qc-summary/panel-qc-summary.tsv",
        "qc-summary/panel-qc-summary.xlsx",
    benchmark:
        ".log/qc/all_summary.bm"
    log:
        ".log/qc/all_summary.log",
    conda:
        config["conda"]["python"]
    script:
        "../scripts/all_qc_summary.py"
