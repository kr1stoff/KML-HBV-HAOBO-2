from pathlib import Path
from src.kml_hbv_haobo_2.fastq import prepare_fastq_by_samptab
from src.kml_hbv_haobo_2.snakemake import create_snakemake_configfile


input_tab = '/data/mengxf/GitHub/KML-HBV-HAOBO-2/tests/input-1800bp.tsv'
output_dir = '/data/mengxf/Project/KML251013-HAOBOHBV-PIPE-UPDATE/results/251020-4'
threads = 32
fb_para_num = 0


output_dir = Path(output_dir).resolve()
# fastq
prepare_fastq_by_samptab(output_dir, input_tab, threads)
# snakemake
create_snakemake_configfile(input_tab, output_dir, threads, fb_para_num)
