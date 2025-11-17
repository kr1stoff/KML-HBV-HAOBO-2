import pandas as pd
from functools import reduce
import sys

sys.stderr = open(snakemake.log[0], "w")


# 参考的位点
known_sites = pd.read_csv(snakemake.input[1])
known_sites['key'] = known_sites['Chrom'] + '_' + \
    known_sites['Pos'].astype(str) + '_' + known_sites['Ref'] + '_' + known_sites['Alt']
# 样本总表
rawtab = snakemake.input[0]
dfs = []
dfs.append(known_sites.drop(columns=['key']))
alldf = pd.ExcelFile(rawtab)
for sn in alldf.sheet_names:
    df = pd.read_excel(rawtab, sheet_name=sn,
                       usecols=['Chrom', 'Pos', 'Ref', 'Alt', 'AltDepth', 'AltFreq'])
    df['key'] = df['Chrom'] + '_' + df['Pos'].astype(str) + '_' + df['Ref'] + '_' + df['Alt']
    target_df = df[df['key'].isin(known_sites['key'])]
    target_df = target_df[target_df['AltDepth'] > int(snakemake.params.core_depth_cutoff)]
    non_target_df = df[~df['key'].isin(known_sites['key'])]
    non_target_df = df[df['AltDepth'] > int(snakemake.params.other_depth_cutoff)]
    # 合并
    df = pd.concat([target_df, non_target_df], axis=0)
    df.rename(columns={'AltDepth': f'Depth-{sn}', 'AltFreq': f'Freq-{sn}'}, inplace=True)
    df.drop(columns=['key'], inplace=True)
    dfs.append(df)
merged = reduce(lambda x, y: pd.merge(x, y, on=['Chrom', 'Pos', 'Ref', 'Alt'], how='outer'), dfs)
merged.drop_duplicates(inplace=True)
merged.to_csv(snakemake.output[0], index=False, sep='\t')
merged.to_excel(snakemake.output[1], index=False)
