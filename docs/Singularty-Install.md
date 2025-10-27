# Singularty 部署 KML-HBV-HAOBO-2

1. 创建沙盒容器

    ```bash
    singularity build --fakeroot --sandbox kml-haobo-hbv-ubuntu22.04-sandbox /data/mengxf/GitHub/KML-HBV-HAOBO-2/singularity.def
    ```

2. 进入容器

    ```bash
    singularity shell --fakeroot --containall --writable --bind /data/mengxf:/mnt kml-haobo-hbv-ubuntu22.04-sandbox
    ```

3. 安装 Conda 环境

    ```bash
    # 激活 conda 环境
    source /opt/miniconda3/etc/profile.d/conda.sh
    # snakemake 环境
    mamba create -n snakemake -y bioconda::snakemake=9.6.0
    # python3.12 环境
    mamba create -n python3.12 -y python=3.12
    mamba install -n python3.12 -c bioconda -c conda-forge -y biopython click numpy pandas pyyaml scipy vcfpy statsmodels openpyxl
    # basic 环境
    mamba create -n basic -y python=3.8
    mamba install -n basic -c bioconda -c conda-forge -y fastqc multiqc bwa samtools bedtools fastp freebayes csvtk bcftools tantan ivar
    ```

4. 复制流程代码库

   ```bash
   cp -r /mnt/GitHub/KML-HBV-HAOBO-2 /opt/
   ```

   修改流程配置文件. `/opt/KML-HBV-HAOBO-2/src/config` 目录

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

5. 包装成一个程序

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

6. 运行测试

   ```bash
   singularity exec --containall --bind /data:/data kml-haobo-hbv-ubuntu22.04-sandbox bash -c "KML-HBV-HAOBO-2 --input-tab /data/mengxf/GitHub/KML-HBV-HAOBO-2/tests/input-1800bp.tsv --output-dir /data/mengxf/Project/KML251013-HAOBOHBV-PIPE-UPDATE/results/251023 --threads 32"
   ```

7. sandbox 转 SIF

    ```bash
    singularity build --fakeroot kml-haobo-hbv-ubuntu22.04.sif kml-haobo-hbv-ubuntu22.04-sandbox
    ```

8. 运行项目数据
  
   ```bash
   singularity exec --containall --bind /data:/data kml-haobo-hbv-ubuntu22.04.sif bash -c "KML-HBV-HAOBO-2 --input-tab /data/mengxf/GitHub/KML-HBV-HAOBO-2/tests/input-1800bp.tsv --output-dir /data/mengxf/Project/KML251013-HAOBOHBV-PIPE-UPDATE/results/251023 --threads 32"
   ```
