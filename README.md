# KML-HBV-HAOBO-2

## 分析

  ```bash
  /home/mengxf/miniforge3/envs/python3.12/bin/python -m src.kml_hbv_haobo_2 \
    --input-tab /data/mengxf/Project/KML250829-HBVHAOBO-HWWKCAFX7/work/250829-input/input.hbv.tsv \
    --output-dir /data/mengxf/Project/KML250829-HBVHAOBO-HWWKCAFX7/results/250829 \
    --threads 32
  ```

## 同步结果

- 分析结果

  ```bash
  rsync -auvP --delete --exclude '**.bam' --exclude '**.gz' \
    /data/mengxf/Project/KML250829-HBVHAOBO-HWWKCAFX7/results/ \
    /data/share/samba/public/bioinfo/KML250829-HBVHAOBO-HWWKCAFX7/results/
  ```

- SAV

  ```bash
  rsync -auvP --delete \
    --include 'Images/***' --include 'InterOp/***' --include 'Thumbnail_Images/***' \
    --include 'RunInfo.xml' --include 'RunParameters.xml' \
    --exclude '*' \
    /data/rawdata/illumina/NEXTseq500/250827_NB501947_0947_AHWWKCAFX7/ \
    /data/share/samba/public/bioinfo/KML250829-HBVHAOBO-HWWKCAFX7/250827_NB501947_0947_AHWWKCAFX7/
  ```

## 注意

1. 不再使用 poetry, 在 snakemake 中总是影响环境, 报 numpy 和 pandas 版本冲突, 直接使用 python 解释器
