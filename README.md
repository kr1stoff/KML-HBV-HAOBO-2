# KML-HBV-HAOBO-2

## 分析

  ```bash
  /home/mengxf/miniforge3/envs/python3.12/bin/python -m src.kml_hbv_haobo_2 --input-tab /data/mengxf/Project/KML250829-HBVHAOBO-HWWKCAFX7/work/250829-input/input.hbv.tsv --output-dir /data/mengxf/Project/KML250829-HBVHAOBO-HWWKCAFX7/results/250829 --threads 32
  ```

## 软件环境部署

1. 安装 miniforge3, 并复制环境  

    - 安装 miniconda3

        ```bash
        # 一路默认确定
        bash Miniforge3-Linux-x86_64.sh
        ```

    - 复制环境

        ```bash
        # python3.12
        mkdir -p ~/miniforge3/envs/python3.12
        tar xvf python3.12.tar.gz -C ~/miniforge3/envs/python3.12
        conda activate python3.12
        conda-unpack

        # basic 同样
        # basic2 同样
        ```

2. 放置软件及数据库  
假设用户目录为 `/data/kml`

    - 软件

      ```bash
      mkdir -p /data/kml/release/
      cp KML-HBV-HAOBO-2.tar /data/kml/release/
      cd /data/kml/release/
      tar xvf KML-HBV-HAOBO-2.tar
      ```

    - 数据库

      ```bash
      mkdir -p /data/kml/database/
      cp D00330.fa /data/kml/database/
      ```

3. 对应修改 `KML-HBV-HAOBO-2\src\config` 中配置文件的对应路径
    - database.py  
    ref(需要提前创建 bwa/samtools 索引), bed, known_sites
    - env.py 无需更改
    - software.py  
    ACTIVATE

## 同步结果

- 分析结果

  ```bash
  rsync -auvP --delete --exclude '**.bam' --exclude '**.gz' /data/mengxf/Project/KML250829-HBVHAOBO-HWWKCAFX7/results/ /data/share/samba/public/bioinfo/KML250829-HBVHAOBO-HWWKCAFX7/results/
  ```

- SAV

  ```bash
  rsync -auvP --delete --include 'Images/***' --include 'InterOp/***' --include 'Thumbnail_Images/***' --include 'RunInfo.xml' --include 'RunParameters.xml' --exclude '*' /data/rawdata/illumina/NEXTseq500/250827_NB501947_0947_AHWWKCAFX7/ /data/share/samba/public/bioinfo/KML250829-HBVHAOBO-HWWKCAFX7/250827_NB501947_0947_AHWWKCAFX7/
  ```

## 注意

1. 不再使用 poetry, 在 snakemake 中总是影响环境, 报 numpy 和 pandas 版本冲突, 直接使用 python 解释器
