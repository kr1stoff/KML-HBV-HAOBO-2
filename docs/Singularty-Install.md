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
    - 容器内提升到 root 权限
    - 写入模式，允许安装软件、修改配置、保存文件
    - 挂载本地目录

    ```bash
    singularity shell --fakeroot --no-home --writable --bind /data/mengxf:/mnt ubuntu22.04_sandbox
    ```

4. 容器内安装必要软件  

    ```bash
    apt update
    apt install -y vim less tzdata ca-certificates
    ```

    注:
    - tzdata 需要交互式配置时区 ASIA/Shanghai
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

    - python3.12 环境

        ```bash
        mamba create -n python3.12 -y python=3.12
        mamba install -n python3.12 -c bioconda -c conda-forge -y biopython click numpy pandas pyyaml scipy vcfpy statsmodels
        ```

    - basic 环境

        ```bash
        mamba create -n basic -y python=3.8
        mamba install -n basic -c bioconda -c conda-forge -y fastqc multiqc bwa samtools bedtools fastp freebayes csvtk bcftools tantan ivar
        ```

    - snakemake 环境

        ```bash
        mamba create -n snakemake -y snakemake=9.6.0
        ```
