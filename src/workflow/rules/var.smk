rule freebayes:
    input:
        alns=rules.bedtools_intersect.output.bam,
        ref=rules.bedtools_maskfasta.output,
        fai=rules.lowcov_samtools_faidx.output,
        targets=rules.make_target_bed.output,
    output:
        vcf="variant/{sample}.vcf",
    log:
        ".log/variant/{sample}.freebayes.log",
    benchmark:
        ".log/variant/{sample}.freebayes.bm"
    # * 用 thread 控制 freebayes 的并行数量，减小内存压力
    threads:
        config["custom"]["freebayes_threads"]
    # * [20251014] 浩博建议的参数
    params:
        (
            "--pooled-continuous --min-repeat-size 10 --read-indel-limit 15 --use-best-n-alleles 4 "
            "--theta 0.005 --haplotype-length 0 --min-alternate-fraction 0.005 --min-base-quality 30 "
            "--min-coverage 1000 --min-alternate-count 10 --min-mapping-quality 30 --max-complex-gap 1 --trim-complex-tail"
        ),
    conda:
        config["conda"]["freebayes"]
    shell:
        "freebayes {params} --targets {input.targets} --fasta-reference {input.ref} {input.alns} > {output.vcf} 2> {log}"

# Description: 变异过滤(错误模型)第一步
#               1. 低质量: 深度 < 1000X, 等位基因深度 < 10X
#               2. 低于检测限: 最小等位基因频率 < 1%
#               3. 链偏倚: 正链/负链支持的 reads 数量 < 2, 正链/负链支持的 reads < 该变异总 reads 的 10%
#               4. 位置偏倚: 变异位置左侧/右侧支持的测序深度 < 2, 变异位置左侧/右侧支持的测序深度 < 该变异总深度的 10%
# Date: 20251015
rule filter_variant1:
    input:
        rules.freebayes.output,
    output:
        "variant/{sample}.filter1.vcf",
    log:
        ".log/variant/{sample}.filter_variant1.log",
    benchmark:
        ".log/variant/{sample}.filter_variant1.bm"
    conda:
        config["conda"]["bcftools"]
    shell:
        """
        bcftools filter --soft-filter LowQual --exclude 'INFO/DP < 1000 || INFO/AO < 10' {input} | \
        bcftools filter --mode + --soft-filter BelowLOD --exclude 'INFO/AO / INFO/DP < 0.01' | \
        bcftools filter --mode + --soft-filter StrandBias --exclude 'SAF < 2 || SAR < 2' | \
        bcftools filter --mode + --soft-filter StrandBias --exclude '(SAF + SAR) > 0 && SAF / (SAF + SAR) < 0.1' | \
        bcftools filter --mode + --soft-filter StrandBias --exclude '(SAF + SAR) > 0 && SAR / (SAF + SAR) < 0.1' | \
        bcftools filter --mode + --soft-filter PosBias --exclude 'RPL < 2 || RPR < 2' | \
        bcftools filter --mode + --soft-filter PosBias --exclude '(RPL + RPR) < 2 || RPL / (RPL + RPR) < 0.1' | \
        bcftools filter --mode + --soft-filter PosBias --exclude '(RPL + RPR) < 2 || RPR / (RPL + RPR) < 0.1' > {output}
        """


# Description: 变异过滤(错误模型)第二步
#              预先注释好基因组低复杂度 BED 区域, 然后用 soft-filter 方式添加上 LowComplexity 标签
# Date: 20251015
rule filter_variant2:
    input:
        rules.filter_variant1.output,
        rules.extract_lcr_bed.output,
    output:
        "variant/{sample}.filter.vcf",
    log:
        ".log/variant/{sample}.filter_variant2.log",
    benchmark:
        ".log/variant/{sample}.filter_variant2.bm"
    conda:
        config["conda"]["python"]
    script:
        "../scripts/vcf_low_lcr.py"


rule vcf2tab:
    input:
        rules.filter_variant2.output,
    output:
        "variant/{sample}.raw.tsv",
    log:
        ".log/variant/{sample}.vcf2tab.log",
    benchmark:
        ".log/variant/{sample}.vcf2tab.bm"
    conda:
        config["conda"]["python"]
    script:
        "../scripts/vcf2tab.py"

# Description: 变异检验
#               1. 二项分布, scipy.stats.binom_test, 使用 测序深度, 最小等位基因频率, 测序错误率
#               2. 泊松分布, scipy.stats.poisson_test, 使用 测序深度, 最小等位基因频率, 测序错误率
# Date: 20251015
rule variant_test:
    input:
        rules.vcf2tab.output,
    output:
        "variant/{sample}.tsv",
    log:
        ".log/variant/{sample}.variant_test.log",
    benchmark:
        ".log/variant/{sample}.variant_test.bm"
    conda:
        config["conda"]["python"]
    script:
        "../scripts/var_test.py"


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
