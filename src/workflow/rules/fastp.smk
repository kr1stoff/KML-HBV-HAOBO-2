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
        config["conda"]["fastp"]
    threads: config["threads"]["low"]
    params:
        # * [20251013] 浩博 (ChatGPT) 建议的参数
        extra=(
            "--cut_front --cut_tail --cut_mean_quality 20 --qualified_quality_phred 20 "
            "--length_required 75 --detect_adapter_for_pe --trim_poly_g"
        ),
    wrapper:
        f"file:{workflow.basedir}/wrappers/fastp"


rule cleaned_fastp:
    input:
        sample=[
            "fastp/{sample}.1.fastq.gz",
            "fastp/{sample}.2.fastq.gz",
        ]
    output:
        trimmed=[temp("fastp/{sample}.trimmed.1.fastq.gz"), temp("fastp/{sample}.trimmed.2.fastq.gz")],
        html=temp("fastp/{sample}.trimmed.html"),
        json="fastp/{sample}.trimmed.json",
    log:
        ".log/fastp/{sample}.cleaned_fastp.log",
    benchmark:
        ".log/fastp/{sample}.cleaned_fastp.bm"
    conda:
        config["conda"]["fastp"]
    threads: config["threads"]["low"]
    params:
        extra=""
    wrapper:
        f"file:{workflow.basedir}/wrappers/fastp"


rule fq_stats_summary:
    input:
        fastp_json=expand("fastp/{sample}.json", sample=config["samples"]),
        cleaned_json=expand("fastp/{sample}.trimmed.json", sample=config["samples"]),
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
