# KML-HBV-HAOBO-2

## 分析

- 主流程

  ```bash
  /home/mengxf/miniforge3/envs/python3.12/bin/python -m src.kml_hbv_haobo_2 \
    --input-tab /data/mengxf/Project/KML250829-HBVHAOBO-HWWKCAFX7/work/250829-input/input.hbv.tsv \
    --output-dir /data/mengxf/Project/KML250829-HBVHAOBO-HWWKCAFX7/results/250829
  ```

- snakemake 运行

  ```bash
  snakemake -c 32 --use-conda -s /data/mengxf/GitHub/KML-HBV-HAOBO-2/src/workflow/Snakefile --configfile .temp/snakemake.yaml --scheduler greedy
  ```

- Singularity 容器运行

  ```bash
  singularity exec --containall \
    --bind /data:/data /data/mengxf/Software/Singularity/kml-haobo-hbv-ubuntu22.04.sif bash \
    -c "KML-HBV-HAOBO-2 \
     --input-tab /data/mengxf/GitHub/KML-HBV-HAOBO-2/tests/input-1800bp.tsv \
     --output-dir /data/mengxf/Project/KML251013-HAOBOHBV-PIPE-UPDATE/results/251023-2"
  ```

## 辅助工具

- 内存使用统计

    ```bash
    /home/mengxf/miniforge3/envs/python3.12/bin/python -m src.util.memory_stats \
      --resdir /data/mengxf/Project/KML251030-HAOBOHBV-HWWV3AFX7/results/251031-bds-4 \
      --outfile /data/mengxf/Project/KML251030-HAOBOHBV-HWWV3AFX7/results/251031-bds-4/memory_stats.tsv
    ```

## 同步结果

- 分析结果

  ```bash
  rsync -auvP --delete --exclude '**.bam' --exclude '**.gz' /data/mengxf/Project/KML250829-HBVHAOBO-HWWKCAFX7/results/ /data/share/samba/public/bioinfo/KML250829-HBVHAOBO-HWWKCAFX7/results/
  ```

- SAV

  ```bash
  rsync -auvP --delete --include 'Images/***' --include 'InterOp/***' --include 'Thumbnail_Images/***' --include 'RunInfo.xml' --include 'RunParameters.xml' --exclude '*' /data/rawdata/illumina/NEXTseq500/250827_NB501947_0947_AHWWKCAFX7/ /data/share/samba/public/bioinfo/KML250829-HBVHAOBO-HWWKCAFX7/250827_NB501947_0947_AHWWKCAFX7/
  ```

## 更新

- [20251114]
  - 修改 vcf 过滤部分  
    - `bcftools norm` + `vt decompose_blocksub` 把多等位基因 MNP 转换为 单等位基因 SNP
    - 重复的 SNP 合并 AO, SAF/SAR, RPL/RLR，reads 是加和关系
    - vcf filter 标签从 bcftools filter 转为自建脚本，LowComplexity 标签脚本也整合进来
  - 新增软件 vt 对应修改 singularity 中 conda 环境
  
- [20251029] 新增支持输入分型
  - 修复 target 区域错误的 BUG
  - 新增支持输入分型, 并根据分型调整参考序列

- [20251023] 新增 Singularity 容器支持
  - 容器所在地址为 /data/mengxf/Software/Singularity/kml-haobo-hbv-ubuntu22.04.sif

- [20251021] 根据浩博脚本更新分析流程
  - 使用 bwa 比对引物位置, 输出引物位置, 引物flank 50bp屏蔽区域, 靶区域
  - ivar 去引物
  - bedtools intersect 去除引物flank 50bp屏蔽区域reads

- [20251013] 根据浩博的修改意见修改流程和参数
  - 修改 fastp 参数 --cut_front --cut_tail --cut_mean_quality 20 --qualified_quality_phred 20 --length_required 75 --detect_adapter_for_pe --trim_poly_g
  - 新增 fastp 接头残留计算
  - 低深度区域屏蔽深度阈值20调整至1000
  - 新增比对质控均一性 P90/P10 数值
  - 靶区域从 1300-1800bp 调整至 70-1730bp (排除引物+flank区域70bp)
  - 修改 freebayes 参数 --pooled-continuous --min-repeat-size 10 --read-indel-limit 15 --use-best-n-alleles 4 --theta 0.005 --haplotype-length 0 --min-alternate-fraction 0.005 --min-base-quality 30 --min-coverage 1000 --min-alternate-count 10 --min-mapping-quality 30 --max-complex-gap 1 --trim-complex-tail
  - 加入错误模型过滤 深度,最小等位基因深度,变异频率,链偏倚,位置偏倚,低复杂区域变异

- [20251010]
  - 新增 freebayes 线程数参数, 默认为 0, 表示和线程数相同

## 注意

- 不再使用 poetry, 在 snakemake 中总是影响环境, 报 numpy 和 pandas 版本冲突, 直接使用 python 解释器
