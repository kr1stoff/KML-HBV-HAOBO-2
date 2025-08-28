import sys
import csv
import pandas as pd

sys.stderr = open(snakemake.log[0], "w")

vcf_file = snakemake.input[0]
out_file = snakemake.output[0]

with open(vcf_file) as f:
    lines = (line for line in f.readlines() if not line.startswith("##"))
reader = csv.DictReader(lines, delimiter="\t")
out_recs = []
for row in reader:
    chrom = row["#CHROM"]
    pos = int(row["POS"])
    ref = row["REF"]
    alt = row["ALT"]
    info = row["INFO"]
    # 如果过 info 中没有 AO 字段，则跳过该行
    if "AO" not in info:
        continue
    info_dict = dict(item.split("=") for item in info.split(";"))
    ao = info_dict["AO"]
    dp = int(info_dict["DP"])
    # 多 ALT 拆成多条
    if "," in ao:
        aos = [int(x) for x in ao.split(",")]
        alts = alt.split(",")
        for ia in range(len(aos)):
            curalt = alts[ia]
            if (len(ref) > 1) and (len(ref) == len(curalt)):
                for ir in range(len(ref)):
                    if ref[ir] != curalt[ir]:
                        out_recs.append([chrom, pos + ir, ref[ir], curalt[ir],
                                         aos[ia], dp, aos[ia] / dp])
            else:
                out_recs.append([chrom, pos, ref, alts[ia], aos[ia], dp, aos[ia] / dp])
    else:
        freq = int(ao) / dp
        # MNP 拆成单个 SNP
        if (len(ref) > 1) and (len(ref) == len(alt)):
            for i in range(len(ref)):
                if ref[i] != alt[i]:
                    out_recs.append([chrom, pos + i, ref[i], alt[i], ao, dp, freq])
        # 正常情况
        else:
            out_recs.append([chrom, pos, ref, alt, ao, dp, freq])
# 使用 Dataframe 处理复杂情况
df = pd.DataFrame(out_recs,
                  columns=["Chrom", "Pos", "Ref", "Alt", "Alt_Depth", "Total_Depth", "Alt_Freq"])
# 去重, sum 相同位置的变异的 Depth 和 Freq
dfgrp = df.groupby(['Chrom', 'Pos', 'Ref', 'Alt']).sum().reset_index()
dfgrp['Alt_Freq'] = dfgrp['Alt_Freq'].apply(lambda x: round(x, 4))
dfgrp.to_csv(out_file, index=False, sep='\t')
