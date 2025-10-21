# 巢氏 PCR, 内引物扩增产物为文库实际区域
import sys

sys.stderr = open(snakemake.log[0], "w")


# * IO
in_bed = snakemake.input[0]
out_bed = snakemake.output[0]

chroms, plus_pos, minus_pos = [], [], []
with open(in_bed) as f:
    for line in f:
        chrom, start, end, _, _, strand = line.strip().split("\t")
        chroms.append(chrom)
        if strand == "+":
            plus_pos.append(int(end))
        else:
            minus_pos.append(int(start))

with open(out_bed, "w") as f:
    f.write(f"{chroms[0]}\t{max(plus_pos)}\t{min(minus_pos)}\n")
