import yaml
import logging
from pathlib import Path
from subprocess import run

from src.kml_hbv_haobo_2.args import Argument
from src.kml_hbv_haobo_2.fastq import get_sample_names_by_samptab
from src.kml_hbv_haobo_2.config import get_thread_dict, get_custom_params, database_update_ref
from src.config.env import CONDA_ENV_DICT
from src.config.software import SNAKEMAKE


def create_snakemake_configfile(args: Argument) -> str:
    """创建 snakemake 配置文件"""
    logging.info('创建 snakemake 配置文件')
    samples = get_sample_names_by_samptab(args.input_tab)
    tempdir = args.output_dir.joinpath('.temp')
    tempdir.mkdir(exist_ok=True, parents=True)
    dict_smk = {
        'workdir': str(args.output_dir),
        'samples': samples,
        'threads': get_thread_dict(args.threads),
        'conda': CONDA_ENV_DICT,
        'database': database_update_ref(args.genotype),
        'custom': get_custom_params(args.threads, args.fb_para_num),
    }
    configfile = f'{args.output_dir}/.temp/snakemake.yaml'
    with open(configfile, 'w') as f:
        yaml.dump(dict_smk, f)
    return configfile


def run_snakemake(args: Argument) -> None:
    """运行 snakemake 工作流"""
    logging.info('运行 snakemake')
    # 创建 snakemake 配置文件
    configfile = create_snakemake_configfile(args)
    # 运行 snakemake 流程
    snakefile = Path(__file__).resolve().parents[1].joinpath('workflow/Snakefile')
    logfile = f'{args.output_dir}/.temp/snakemake.log'
    cml = f"{SNAKEMAKE} -c {args.threads} --use-conda -s {snakefile} --configfile {configfile} --ignore-incomplete --scheduler greedy"
    logging.debug(cml)
    proc = run(cml, shell=True, executable='/bin/bash', capture_output=True, encoding='utf-8')
    # 输出出来这段日志
    with open(logfile, 'w') as f:
        f.write(f'[STDOUT]\n{proc.stdout}\n[STDERR]\n{proc.stderr}')
