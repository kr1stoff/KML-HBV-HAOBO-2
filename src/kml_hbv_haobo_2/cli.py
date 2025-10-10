import click
import logging
from pathlib import Path

from src.kml_hbv_haobo_2.fastq import prepare_fastq_by_samptab
from src.kml_hbv_haobo_2.snakemake import run_snakemake
from src.kml_hbv_haobo_2 import __version__

logging.basicConfig(level=logging.DEBUG, format='%(asctime)s - %(levelname)s - %(message)s')


@click.command()
@click.option('--input-tab', required=True, help='输入样本表, 包含样本名, read1 和 read2 的路径')
@click.option('--output-dir', default='kml-hbv-haobo-result', show_default=True, help='输出文件夹')
@click.option('--threads', default=8, type=int, show_default=True, help='线程数')
@click.option('--freebayes-para-num', default=0, type=int, show_default=True, help='freebayes 线程数参数, 0 表示和线程数相同')
@click.version_option(version=__version__, message='KML HBV HAOBO Version: %(version)s', help='显示版本信息')
@click.help_option(help='获取帮助信息')
def main(input_tab, output_dir, threads, freebayes_para_num):
    """KML 浩博 HBV X 基因变异分析流程"""
    logging.info("开始 KML 浩博 HBV X 基因变异分析流程")
    output_dir = Path(output_dir).resolve()
    # fastq
    prepare_fastq_by_samptab(output_dir, input_tab, threads)
    # snakemake
    run_snakemake(input_tab, output_dir, threads, freebayes_para_num)
