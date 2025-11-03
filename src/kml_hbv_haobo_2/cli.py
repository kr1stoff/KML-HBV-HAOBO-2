import click
import logging

from src.kml_hbv_haobo_2.args import Argument
from src.kml_hbv_haobo_2.fastq import prepare_fastq_by_samptab
from src.kml_hbv_haobo_2.snakemake import run_snakemake
from src.kml_hbv_haobo_2.system import get_default_threads
from src.kml_hbv_haobo_2 import __version__

logging.basicConfig(level=logging.DEBUG, format='%(asctime)s - %(levelname)s - %(message)s')


@click.command()
@click.option('--input-tab', required=True, help='输入样本表, 包含样本名, read1 和 read2 的路径')
@click.option('--output-dir', default='kml-hbv-haobo-result', show_default=True, help='输出文件夹')
@click.option('--genotype', default='B', show_default=True, type=click.Choice(['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H']), help='参考基因型')
@click.option('--freebayes-para-num', default=0, type=int, show_default=True, help='freebayes 线程数参数, 0 表示和线程数相同')
@click.option('--threads', default=get_default_threads(), type=int, show_default=True, help='线程数')
@click.version_option(version=__version__, message='KML HBV HAOBO Version: %(version)s', help='显示版本信息')
@click.help_option(help='获取帮助信息')
def main(input_tab, output_dir, threads, freebayes_para_num, genotype):
    """KML 浩博 HBV 基因变异分析流程"""
    logging.info("开始 KML 浩博 HBV 基因变异分析流程")

    # 解析参数
    args = Argument(
        input_tab=input_tab,
        output_dir_str=output_dir,
        threads=threads,
        fb_para_num=freebayes_para_num,
        genotype=genotype,
    )
    # fastq
    prepare_fastq_by_samptab(args)
    # snakemake
    run_snakemake(args)

    logging.info("KML 浩博 HBV 基因变异分析流程完成")
