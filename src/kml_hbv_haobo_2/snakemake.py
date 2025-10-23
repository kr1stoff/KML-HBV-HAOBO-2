import yaml
import logging
from pathlib import Path
from subprocess import run

from src.kml_hbv_haobo_2.fastq import get_sample_names_by_samptab
from src.kml_hbv_haobo_2.config import get_thread_dict, get_custom_params
from src.config.env import CONDA_ENV_DICT
from src.config.software import SNAKEMAKE
from src.config.database import DATABASE


def create_snakemake_configfile(input_tab: str, workdir: Path, threads: int, freebayes_para_num: int) -> str:
    """创建 snakemake 配置文件"""
    logging.info('创建 snakemake 配置文件')
    samples = get_sample_names_by_samptab(input_tab)
    tempdir = workdir.joinpath('.temp')
    tempdir.mkdir(exist_ok=True, parents=True)
    dict_smk = {
        'workdir': str(workdir),
        'samples': samples,
        'threads': get_thread_dict(threads),
        'conda': CONDA_ENV_DICT,
        'database': DATABASE,
        'custom': get_custom_params(threads, freebayes_para_num),
    }
    configfile = f'{workdir}/.temp/snakemake.yaml'
    with open(configfile, 'w') as f:
        yaml.dump(dict_smk, f)
    return configfile


def run_snakemake(input_tab: str, workdir: Path, threads: int, freebayes_para_num: int) -> None:
    """
    运行 snakemake 工作流
    :param input_tab:       样本信息表, 样本名/read1/read2
    :param workdir:         分析结果目录
    :param threads:         最大线程数
    :param reference:       参考基因组路径
    :param bed:             BED 文件路径
    :return:
    """
    logging.info('运行 snakemake')
    # 创建 snakemake 配置文件
    configfile = create_snakemake_configfile(input_tab, workdir, threads, freebayes_para_num)
    # 运行 snakemake 流程
    snakefile = Path(__file__).resolve().parents[1].joinpath('workflow/Snakefile')
    logfile = f'{workdir}/.temp/snakemake.log'
    cml = f"{SNAKEMAKE} -c {threads} --use-conda -s {snakefile} --configfile {configfile} --ignore-incomplete --scheduler greedy"
    logging.debug(cml)
    proc = run(cml, shell=True, executable='/bin/bash', capture_output=True, encoding='utf-8')
    # 输出出来这段日志
    with open(logfile, 'w') as f:
        f.write(f'[STDOUT]\n{proc.stdout}\n[STDERR]\n{proc.stderr}')
