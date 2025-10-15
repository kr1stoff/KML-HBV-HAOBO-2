# 使用低复杂度区域(Low Complexity Region) soft-filter 方式添加 LowComplexity 标签
import vcfpy
import sys

sys.stderr = open(snakemake.log[0], "w")

vcf_file = snakemake.input[0]
low_cpx_bed = snakemake.input[1]
out_file = snakemake.output[0]

# 获取低复杂度区域 BED 区域
low_cpx_regions = {}
with open(low_cpx_bed, "r") as f:
    for line in f:
        line = line.strip()
        fields = line.split("\t")
        chrom = fields[0]
        start = int(fields[1])
        end = int(fields[2])
        low_cpx_regions[chrom] = low_cpx_regions.get(chrom, []) + [(start, end)]

# 打开 VCF 文件, 添加 Header 信息
reader = vcfpy.Reader.from_path(vcf_file)
reader.header.add_filter_line({
    "ID": "LowComplexity",
    "Description": "Variant is in simple regions / low complexity / tandem repeats region"
})
writer = vcfpy.Writer.from_path(out_file, reader.header)

# 遍历 VCF 文件, 检查每个 variant 是否在低复杂度区域
for record in reader:
    if record.CHROM in low_cpx_regions:
        for region in low_cpx_regions[record.CHROM]:
            # BED 区域是 0-based, 而 VCF 是 1-based, 所以需要减 1
            if region[0] <= (record.POS - 1) <= region[1]:
                record.FILTER.append("LowComplexity")
                break
    writer.write_record(record)
