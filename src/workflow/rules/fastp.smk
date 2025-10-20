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
    shell:
        """
        fastp -w {threads} {params.extra} \
            -i {input.sample[0]} -I {input.sample[1]} \
            -o {output.trimmed[0]} -O {output.trimmed[1]} \
            -h {output.html} -j {output.json} 2> {log}
        """


# Description:  对 cleaned fastq 进行 fastp 质控, 查看 adapter 残留比例.
#               若 cleaned fastq 为空, 则创建空的 trimmed fastq 文件, html 文件和 json 文件.
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
    run:
        import gzip
        import os
        if (os.stat(input[0]).st_size <= 100) or (os.stat(input[1]).st_size <= 100):
            # 创建空的 trimmed fastq 文件
            for out_file in output.trimmed:
                with gzip.open(out_file, "wb") as f:
                    f.write(b"")
            # 创建空的 html 文件
            with open(output.html, "w") as f:
                f.write('<html><body>No reads after first filtering</body></html>')
            # 创建空的 json 文件
            with open(output.json, "w") as f:
                f.write('{"adapter_cutting": {"adapter_trimmed_reads": 0}, "summary": {"before_filtering": {"total_reads": 0}}}')
        else:
            # 正常运行 fastp
            shell("""
            fastp -w {threads} {params.extra} \
                -i {input[0]} -I {input[1]} \
                -o {output.trimmed[0]} -O {output.trimmed[1]} \
                -h {output.html} -j {output.json} 2> {log}
            """)


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
