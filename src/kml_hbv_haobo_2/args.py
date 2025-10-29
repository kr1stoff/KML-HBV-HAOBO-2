from dataclasses import dataclass
from pathlib import Path


@dataclass
class Argument:
    input_tab: str
    output_dir_str: str = 'kml-hbv-haobo-result'
    threads: int = 8
    fb_para_num: int = 0
    genotype: str = 'B'

    def __post_init__(self):
        self.output_dir = Path(self.output_dir_str).resolve()
