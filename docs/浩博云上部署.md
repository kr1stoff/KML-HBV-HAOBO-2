# 浩博云上部署

1. 安装 `singularity` 容器

    ```bash
    # 用不了就 scp 过去
    wget https://mirrors.bfsu.edu.cn/github-release/conda-forge/miniforge/LatestRelease/Miniforge3-Linux-x86_64.sh
    bash Miniforge3-Linux-x86_64.sh
    source ~/.bashrc
    # singularity
    mamba create -n singularity conda-forge::singularity=3.8.6 -y
    # 环境变量
    sudo ln -sf /home/ubuntu/miniforge3/envs/singularity/bin/singularity /usr/bin/
    # 重启终端
    ```

2. 配置工作目录

    ```bash
    # 创建工作目录
    mkdir -p /home/ubuntu/KML
    # sif 文件放在该目录
    cd /home/ubuntu/KML
    mkdir -p rawdata results
    ```

3. 运行分析流程

    ```bash
    # 原始数据放在 rawdata 目录下创建子目录, 生成 input.tsv 文件
    # 每次在 /home/ubuntu/KML 目录运行
    sudo singularity exec --containall --bind /home/ubuntu/KML:/home/ubuntu/KML kml-haobo-hbv-ubuntu22.04.sif bash -c "KML-HBV-HAOBO-2 --input-tab /home/ubuntu/KML/rawdata/251107/input.tsv --output-dir /home/ubuntu/KML/results/251107"
    ```

    软件参数:

    ```text
    Options:
      --input-tab TEXT              输入样本表, 包含样本名, read1 和 read2 的路径  [required]
      --output-dir TEXT             输出文件夹  [default: kml-hbv-haobo-result]
      --genotype [A|B|C|D|E|F|G|H]  参考基因型  [default: B]
      --freebayes-para-num INTEGER  freebayes 线程数参数, 0 表示和线程数相同  [default: 0]
      --threads INTEGER             线程数  [default: 32]
      --version                     显示版本信息
      --help                        获取帮助信息
    ```

注意：

- 架构要一致, 用 `x86` 架构
