# 使用 tantan 获取低复杂度区域(Low Complexity Region)
import click
import logging
import subprocess
from pathlib import Path
from Bio import SeqIO

from src.config.software import TANTAN

logging.basicConfig(level=logging.DEBUG, format='%(asctime)s - %(levelname)s - %(message)s')


@click.command()
@click.option('--input-fasta', required=True, help='输入 FASTA 文件路径')
@click.option('--output-bed', default='low_cpx_region.bed', show_default=True, help='输出 BED 文件路径')
@click.option('--min-len', default=6, type=int, show_default=True, help='低复杂区最小长度(>=)')
@click.help_option(help='获取帮助信息')
def extract_tantan_lcr(input_fasta, output_bed, min_len):
    """使用 tantan 获取低复杂度区域(Low Complexity Region)"""
    logging.info("开始使用 tantan 获取低复杂度区域(Low Complexity Region)")
    output_bed = Path(output_bed).resolve()

    # 临时文件
    temp_dir = output_bed.parent / '.tmp'
    temp_dir.mkdir(parents=True, exist_ok=True)

    # tantan 预测低复杂度区域
    cmd = f"{TANTAN} {input_fasta} > {temp_dir / 'tantan.fasta'}"
    logging.debug(f"运行命令: {cmd}")
    subprocess.run(cmd, shell=True, check=True)

    # 解析 tantan 输出, 小写字母表示低复杂度区域
    record = list(SeqIO.parse(temp_dir / 'tantan.fasta', 'fasta'))[0]
    chrom = record.id
    lower_bases = []
    pos = 0
    for base in record.seq:
        if base.islower():
            lower_bases.append(pos)
        pos += 1

    # 转换为 BED 格式
    ranges = to_ranges(lower_bases)
    with open(output_bed, 'w') as f:
        for start, end in ranges:
            # * 限制低复杂区最小长度
            if end - start >= min_len:
                f.write(f"{chrom}\t{start}\t{end}\n")

    # 清理临时文件
    (temp_dir / 'tantan.fasta').unlink()

    logging.info("完成使用 tantan 获取低复杂度区域(Low Complexity Region)")


def to_ranges(nums):
    """把升序整数列表合并成连续区段字符串"""
    if not nums:
        return []
    nums = sorted(set(nums))          # 去重并排序
    ranges, start = [], nums[0]

    for pre, cur in zip(nums, nums[1:] + [None]):
        if cur is None or cur != pre + 1:      # 连续中断
            if start != pre:                   # 区间
                ranges.append([start, pre])
            start = cur
    return ranges


if __name__ == '__main__':
    extract_tantan_lcr()
