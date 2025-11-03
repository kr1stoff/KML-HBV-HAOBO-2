import click
import glob
import subprocess
import pandas as pd
from pathlib import Path


@click.command()
@click.option("--resdir", type=click.Path(exists=True), help="批次分析结果目录")
@click.option("--outfile", default="memory_stats.tsv", help="输出统计文件名")
@click.help_option("--help", help="显示帮助信息")
def main(resdir, outfile):

    resdir = Path(resdir)

    # 数据量
    qc_df = pd.read_csv(resdir / 'qc-summary/panel-qc-summary.tsv', sep='\t',
                        index_col=0, usecols=['Sample', 'RawBases'])
    qc_df['RawBases(GB)'] = round(qc_df['RawBases'] / 1e9, 2)
    qc_df.drop('RawBases', axis=1, inplace=True)

    # 变异数量（MNP和SNP）
    vardict = {}
    for vcf in glob.glob(f'{resdir}/variant/*.vcf'):
        name = Path(vcf).stem
        vcfnum = int(subprocess.run(
            f'grep -c -v "#" {vcf}', shell=True, capture_output=True, text=True).stdout.strip())
        vardict[name] = vcfnum

    # 内存 & 时间
    bmdict = {}
    for bm in glob.glob(f'{resdir}/.log/variant/*.freebayes.bm'):
        name = Path(bm).stem.replace('.freebayes', '')
        bm_df = pd.read_csv(bm, sep='\t')
        mem = bm_df['max_uss'][0]
        bmdict.setdefault(name, {})
        bmdict[name]['mem'] = round(mem/1000, 2)
        _time = bm_df['h:m:s'][0]
        bmdict[name]['time'] = _time

    # 输出
    qc_df['VariantNum'] = [vardict[idx] for idx in qc_df.index]
    qc_df['Memory(GB)'] = [bmdict[idx]['mem'] for idx in qc_df.index]
    qc_df['Time(h:m:s)'] = [bmdict[idx]['time'] for idx in qc_df.index]
    qc_df.to_csv(outfile, sep='\t', index=True)


if __name__ == '__main__':
    main()
