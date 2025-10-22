# KML-HBV-HAOBO-2 软件环境部署

## 创建账户（可选）

使用 *root* 执行

```bash
# 创建
useradd -m test
# 修改密码 test123
passwd test
# 修改 shell 类型
chsh -s /bin/bash test
```

## 创建目录结构

使用 *root* 执行

（可选）数据盘如果不挂在 `/data`，就软链接过去

```bash
ln -sf /dev/shm /data
```

创建目录

````bash
mkdir -p /data/lab # 主目录
mkdir -p /data/lab/rawdata # 原始数据
mkdir -p /data/lab/database # 数据库
mkdir -p /data/lab/release # 软件
mkdir -p /data/lab/download # 其他下载文件
mkdir -p /data/lab/project # 任务目录
chown -R test:test /data/lab # 调整拥有者权限
````

- 上传到 `/data/lab/download`
  - python3.12.tar.gz
  - basic.tar.gz
  - basic2.tar.gz
  - Miniforge3-Linux-x86_64.sh
- 上传到 `/data/lab/database`
  - D00330.fa
- 上传到 `/data/lab/release`
  - KML-HBV-HAOBO-2.tar
- 原始 FQ 数据上传到 `/data/lab/rawdata`

## 安装 miniforge3, 并复制环境

使用 *test* 账户

- 安装 miniconda3

    ```bash
    # 一路默认确定
    cd /data/lab/download
    bash Miniforge3-Linux-x86_64.sh
    ```

- 复制环境

    ```bash
    cd ~/miniforge3/envs
    cp /data/lab/download/*tar.gz .
    # conda pack 包
    ls *gz | sed 's/.tar.gz//g' | while read tgz;do mkdir $tgz && tar zxvf ${tgz}.tar.gz -C $tgz && conda activate $tgz && conda-unpack;done
    rm *tar.gz
    # 在线安装 snakemake 环境，版本 2.3.0
    mamba create -n snakemake bioconda::snakemake
    ```

## 放置软件及数据库

- 软件

   ```bash
   cd /data/lab/release
   tar xvf KML-HBV-HAOBO-2.tar
   ```

- 数据库
   *需要提前使用 bwa, samtools 建立索引*

   ```bash
   cd /data/lab/database
   mamba run -n basic bwa index D00330.fa
   mamba run -n basic samtools faidx D00330.fa
   ```

## 修改软件配置文件

对应修改 `KML-HBV-HAOBO-2\src\config` 中配置文件的对应路径

```bash
cd /data/lab/release
# 按实际路径修改
vi src/config/database.py
# 按实际路径修改
vi src/config/software.py
```

database.py

```python
DATABASE = {
    'ref': '/data/lab/database/D00330.fa',
    'bed': '/data/lab/release/KML-HBV-HAOBO-2/assets/target.bed',
    'known_sites': '/data/lab/release/KML-HBV-HAOBO-2/assets/known_sites.csv',
}
```

env.py 无需更改

software.py

```python
ACTIVATE = '/home/test/miniforge3/bin/activate'
```
