rule fastp:
    input:
        sample=[
            ".rawdata/{sample}_1.fastq.gz",
            ".rawdata/{sample}_2.fastq.gz",
        ],
    output:
        trimmed=["fastp/{sample}.1.fastq.gz", "fastp/{sample}.2.fastq.gz"],
        html="fastp/{sample}.html",
        json="fastp/{sample}.json",
    log:
        ".log/fastp/{sample}.fastp.log",
    benchmark:
        ".log/fastp/{sample}.fastp.bm"
    conda:
        config["conda"]["basic2"]
    threads: config["threads"]["low"]
    params:
        # ! 用浩博的参数
        extra=(
            "--trim_poly_g --qualified_quality_phred 5 --unqualified_percent_limit 50 --n_base_limit 15 "
            "--length_required 150 --overlap_diff_limit 10 --overlap_diff_percent_limit 10"
        ),
    wrapper:
        f"file:{workflow.basedir}/wrappers/fastp"


rule fq_stats_summary:
    input:
        expand("fastp/{sample}.json", sample=config["samples"]),
    output:
        "fastp/fq_summary.tsv",
    benchmark:
        ".log/fastp/fq_all_samples_qc.bm"
    log:
        ".log/fastp/fq_all_samples_qc.log",
    conda:
        config["conda"]["python"]
    script:
        "../scripts/fq_all_samples_qc.py"
