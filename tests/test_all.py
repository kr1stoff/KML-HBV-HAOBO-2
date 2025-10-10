from pathlib import Path
from src.kml_hbv_haobo_2.fastq import prepare_fastq_by_samptab
from src.kml_hbv_haobo_2.snakemake import create_snakemake_configfile


input_tab = '/data/mengxf/Project/KML250929-HAOBO-SERUM-STD/input.2G.tsv'
output_dir = '/data/mengxf/Project/KML250929-HAOBO-SERUM-STD/results/251010'
threads = 32
freebayes_para_num = 4


output_dir = Path(output_dir).resolve()
# fastq
# prepare_fastq_by_samptab(output_dir, input_tab, threads)
# snakemake
create_snakemake_configfile(input_tab, output_dir, threads, freebayes_para_num)
