import sys
from pathlib import Path
import json
import pandas as pd
import numpy as np

sys.stderr = open(snakemake.log[0], "w")

sys.stderr.write(f"Input: {snakemake.input['fastp_json'], snakemake.input['cleaned_json']}\n")

# def fastp_all_samples_qc(fastp_json, out_tsv):
title = ["Sample", "RawReads", "RawBases", "CleanReads", "CleanBases", "RawQ20",
         "RawQ30", "CleanQ20", "CleanQ30", "CleanAverageLength", "GC", "DuplicationRate", "ResidualAdapterRate"]
df = pd.DataFrame(columns=title)

# * 浩博新增接头残留计算
residual_dict = {}
for cj_path in snakemake.input['cleaned_json']:
    sample = Path(cj_path).stem.split('.')[0]
    cj_data = json.loads(open(cj_path, "r").read())
    residual = int(cj_data['adapter_cutting']['adapter_trimmed_reads'])
    all_reads = int(cj_data["summary"]["before_filtering"]["total_reads"])
    residual_rate = round(residual / all_reads, 6) if all_reads > 0 else 0
    residual_dict[sample] = residual_rate

for js_path in snakemake.input['fastp_json']:
    sample = Path(js_path).stem
    js_data = json.loads(open(js_path, "r").read())
    mean_lengths = [v for k, v in js_data["summary"]["after_filtering"].items()
                    if k.endswith("mean_length")]
    out = [
        sample,
        js_data["summary"]["before_filtering"]["total_reads"],
        js_data["summary"]["before_filtering"]["total_bases"],
        js_data["summary"]["after_filtering"]["total_reads"],
        js_data["summary"]["after_filtering"]["total_bases"],
        js_data["summary"]["before_filtering"]["q20_rate"],
        js_data["summary"]["before_filtering"]["q30_rate"],
        js_data["summary"]["after_filtering"]["q20_rate"],
        js_data["summary"]["after_filtering"]["q30_rate"],
        int(np.mean(mean_lengths)),
        js_data["summary"]["after_filtering"]["gc_content"],
        js_data["duplication"]["rate"],
        residual_dict[sample],
    ]
    df.loc[len(df)] = out
df.to_csv(snakemake.output[0], index=False, sep="\t")
