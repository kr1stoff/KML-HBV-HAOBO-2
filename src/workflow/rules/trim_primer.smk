rule ivar_trim:
    input:
        bam=rules.bwa_mem.output,
        bed=rules.make_primer_bed.output.bed,
    output:
        unsorted=temp("trim-primer/{sample}.trimmed.unsorted.bam"),
        bam="trim-primer/{sample}.trimmed.bam",
    log:
        ".log/trim_primer/{sample}.ivar_trim.log",
    benchmark:
        ".log/trim_primer/{sample}.ivar_trim.bm"
    threads: config["threads"]["low"]
    conda:
        config["conda"]["ivar"]
    shell:
        """
        prefix=$(dirname {output.unsorted})/$(basename {output.unsorted} .bam)
        ivar trim -i {input.bam} -b {input.bed} -p $prefix -e 2> {log}

        # ! ivar trim 如果输入没有条目会删掉信息头输出空文件
        # 检查unsorted文件是否有reads(除了信息头)
        if [[ $(samtools view -c {output.unsorted} 2>/dev/null || echo 0) -gt 0 ]]; then
            # 有reads时正常处理
            samtools sort -@ {threads} -o {output.bam} {output.unsorted} 2>> {log}
        else
            # 没有reads时创建空的已排序BAM文件
            echo "警告: {output.unsorted} 没有reads，创建空输出" >> {log}
            samtools view -hbS {input.bam} > {output.bam} 2>> {log}
        fi
        samtools index {output.bam} 2>> {log}
        """


rule bedtools_intersect:
    input:
        bam=rules.ivar_trim.output.bam,
        bed=rules.primer_mask.output.primer,
    output:
        unsorted=temp("trim-primer/{sample}.filtered.unsorted.bam"),
        bam="trim-primer/{sample}.filtered.bam",
    log:
        ".log/trim_primer/{sample}.bedtools_intersect.log",
    benchmark:
        ".log/trim_primer/{sample}.bedtools_intersect.bm"
    threads: config["threads"]["low"]
    conda:
        config["conda"]["bedtools"]
    shell:
        """
        bedtools intersect -v -a {input.bam} -b {input.bed} > {output.unsorted} 2> {log}
        samtools sort -@ {threads} -o {output.bam} {output.unsorted} 2>> {log}
        samtools index {output.bam} 2>> {log}
        """
