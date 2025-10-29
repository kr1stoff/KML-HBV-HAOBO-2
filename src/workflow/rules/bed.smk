rule bwa_mem_primer:
    input:
        primer=config["database"]["primer_fa"],
        ref=rules.copy_ref.output,
        idx=rules.bwa_index.output,
    output:
        temp("prepare/all.sam")
    log:
        ".log/bed/bwa_mem_primer.log"
    benchmark:
        ".log/bed/bwa_mem_primer.bm"
    params:
        "-k 8 -T 8"
    threads: config["threads"]["low"]
    conda:
        config["conda"]["bwa"]
    shell:
        "bwa mem -t {threads} {params} {input.ref} {input.primer} > {output} 2> {log}"


rule make_primer_bed:
    input:
        sam=rules.bwa_mem_primer.output[0],
        fai=rules.samtools_faidx.output[0],
    output:
        mapped=temp("prepare/mapped.sam"),
        bed="prepare/primer.bed"
    log:
        ".log/bed/make_primer_bed.log"
    benchmark:
        ".log/bed/make_primer_bed.bm"
    conda:
        config["conda"]["python"]
    script:
        "../scripts/make_primer_bed.py"


rule primer_mask:
    input:
        bed=rules.make_primer_bed.output.bed,
        fai=rules.samtools_faidx.output,
    output:
        primer="prepare/primer.mask.bed",
        # target="prepare/target.bed",
    log:
        ".log/bed/primer_mask.log"
    benchmark:
        ".log/bed/primer_mask.bm"
    params:
        flank=50,
    conda:
        config["conda"]["bedtools"]
    shell:
        "bedtools slop -b {params.flank} -i {input.bed} -g {input.fai} > {output.primer} 2> {log}"


rule make_target_bed:
    input:
        rules.primer_mask.output.primer
    output:
        "prepare/target.bed",
    log:
        ".log/bed/make_target_bed.log"
    benchmark:
        ".log/bed/make_target_bed.bm"
    conda:
        config["conda"]["python"]
    script:
        "../scripts/make_target_bed.py"
