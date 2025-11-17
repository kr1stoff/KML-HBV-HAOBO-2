import csv
import pandas as pd
import sys
from typing import Any

sys.stderr = open(snakemake.log[0], "w")    # type: ignore

vcf_file = snakemake.input[0]   # type: ignore
out_file = snakemake.output[0]  # type: ignore

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
    fltr = row["FILTER"]
    # 如果过 info 中没有 AO 字段，则跳过该行
    if "AO" not in info:
        continue
    info_dict = dict(item.split("=") for item in info.split(";"))
    ao = info_dict["AO"]
    dp = int(info_dict["DP"])

    # 多个 ALT
    if "," in ao:
        aos = [int(x) for x in ao.split(",")]
        alts = alt.split(",")
        for ia in range(len(aos)):
            curalt = alts[ia]
            # MNP 且 ALT 与 REF 长度相同
            if (len(ref) > 1) and (len(ref) == len(curalt)):
                for ir in range(len(ref)):
                    if ref[ir] != curalt[ir]:
                        out_recs.append(
                            [chrom, pos + ir, ref[ir], curalt[ir],
                                aos[ia], dp, aos[ia] / dp, fltr]
                        )
            # 直接输出
            else:
                out_recs.append([chrom, pos, ref, alts[ia],
                                aos[ia], dp, aos[ia] / dp, fltr])
    # 只有一个 ALT
    else:
        freq = int(ao) / dp
        # MNP 拆成单个 SNP
        if (len(ref) > 1) and (len(ref) == len(alt)):
            for i in range(len(ref)):
                if ref[i] != alt[i]:
                    out_recs.append(
                        [chrom, pos + i, ref[i], alt[i], ao, dp, freq, fltr])
        # 正常情况
        else:
            out_recs.append([chrom, pos, ref, alt, ao, dp, freq, fltr])
# 使用 Dataframe 处理复杂情况
column_names: Any = ["Chrom", "Pos", "Ref", "Alt", "Alt_Depth",
                     "Total_Depth", "Alt_Freq", "VCF_Filter"]
df = pd.DataFrame(
    out_recs,
    columns=column_names
)
# 去重, sum 相同位置的变异的 Depth 和 Freq
dfgrp = df.groupby(["Chrom", "Pos", "Ref", "Alt"]).sum().reset_index()
dfgrp["Alt_Freq"] = dfgrp["Alt_Freq"].apply(lambda x: round(x, 4))
dfgrp.to_csv(out_file, index=False, sep="\t")
