# Singularty 安装流程

1. docker 拉取 Ubuntu 22.04 镜像, 并打包成 `.tar` 文件

    ```bash
    docker pull ubuntu:22.04
    docker save -o ubuntu22.04.tar ubuntu:22.04
    ```

2. 创建沙盒模式容器

    ```bash
    singularity build --sandbox ubuntu22.04_sandbox docker-archive://ubuntu22.04.tar
    ```

3. 开发模式进入容器
    - 避免宿主机 `$HOME` 绑定到容器
    - 容器内提升到 `root` 权限
    - 写入模式，允许安装软件、修改配置、保存文件
    - 挂载本地目录

    ```bash
    singularity shell --fakeroot --no-home --writable --bind /data/mengxf:/mnt ubuntu22.04_sandbox
    ```

4. 容器内安装必要软件  

    ```bash
    apt update
    apt install -y vim less ca-certificates
    apt install -y tzdata
    ```

    注:
    - tzdata 需要交互式配置时区 ASIA/Shanghai (6->70)
    - ca-certificates 用于 HTTPS 下载, 否则 mamba 安装包时会报错

    *可以激活 /root/.bashrc 使用 root 环境变量*

    ```bash
    source /root/.bashrc
    ```

5. 安装 Conda

    ```bash
    cd /mnt/Download
    bash Miniforge3-Linux-x86_64.sh
    ```

    更新 conda, mamba

    ```bash
    mamba update -y conda mamba
    ```

6. 创建环境

    修改环境变量后安装环境

    ```bash
    export MAMBA_ROOT_PREFIX=/root/miniforge3
    export MAMBA_EXE=/root/miniforge3/bin/mamba
    ```

    - python3.12 环境

        ```bash
        mamba create -n python3.12 -y python=3.12
        mamba install -n python3.12 -c bioconda -c conda-forge -y biopython click numpy pandas pyyaml scipy vcfpy statsmodels openpyxl
        ```

    - basic 环境

        ```bash
        mamba create -n basic -y python=3.8
        mamba install -n basic -c bioconda -c conda-forge -y fastqc multiqc bwa samtools bedtools fastp freebayes csvtk bcftools tantan ivar
        ```

    - snakemake 环境

        ```bash
        mamba create -n snakemake -y bioconda::snakemake=9.6.0
        ```

7. 复制流程代码库

   ```bash
   cp -r /mnt/GitHub/KML-HBV-HAOBO-2 /opt/
   ```

   修改流程配置文件

   `/opt/KML-HBV-HAOBO-2/src/config` 目录

   - database.py

     ```python
     DATABASE = {
         'ref': '/opt/KML-HBV-HAOBO-2/assets/D00330/D00330.fasta',
         'known_sites': '/opt/KML-HBV-HAOBO-2/assets/known_sites.csv',
         'primer_fa': '/opt/KML-HBV-HAOBO-2/assets/primer.fasta',
     }
     ```

   - env.py

     ```python
     CONDA_ENV_DICT = {
         'python': 'python3.12',
         'fastqc': 'basic',
         'multiqc': 'basic',
         'bwa': 'basic',
         'samtools': 'basic',
         'bedtools': 'basic',
         'fastp': 'basic',
         'freebayes': 'basic',
         'csvtk': 'basic',
         'bcftools': 'basic',
         'tantan': 'basic',
         'ivar': 'basic',
     }
     ```

   - software.py

     ```python
     SNAKEMAKE = '/opt/miniconda3/envs/snakemake/bin/snakemake'
     ```

8. 包装成一个程序
    - 创建文件 `/bin/KML-HBV-HAOBO-2`

    ```bash
    #!/bin/bash
    
    cd /opt/KML-HBV-HAOBO-2
    /opt/miniconda3/envs/python3.12/bin/python -m src.kml_hbv_haobo_2 "$@"
    ```

    - 赋予可执行权限

    ```bash
    chmod +x /bin/KML-HBV-HAOBO-2
    ```

9. 运行程序

   todo: sandbox 打包回 SIF

   ```bash
   singularity exec --bind /data/mengxf:/data/mengxf ubuntu22.04_sandbox KML-HBV-HAOBO-2 --input-tab /data/mengxf/GitHub/KML-HBV-HAOBO-2/tests/input-1800bp.tsv --output-dir /data/mengxf/Project/KML251013-HAOBOHBV-PIPE-UPDATE/results/251022-2 --threads 32
   ```
