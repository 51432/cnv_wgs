# SCCC Tumor-Normal WGS CNV/SV Pipeline (hg38, SLURM)

用于**宫颈小细胞癌（SCCC）配对 tumor-normal WGS**数据的学术型、模块化分析流程。核心分析为：**ASCAT（CNV）+ Manta/GRIDSS（SV）+ HPV integration 联动**，适配 **SLURM 集群**，支持 array 并行与依赖驱动 DAG。

---

## 目录

- [1. 项目背景与目标](#1-项目背景与目标)
- [2. 项目特点 / 核心功能](#2-项目特点--核心功能)
- [3. 目录结构说明](#3-目录结构说明)
- [4. 输入文件格式](#4-输入文件格式)
- [5. 参考资源文件](#5-参考资源文件)
- [6. 流程总体架构（Phase + Stage）](#6-流程总体架构phase--stage)
- [7. 各模块详细说明](#7-各模块详细说明)
- [8. SLURM 运行方式](#8-slurm-运行方式)
- [9. 配置文件说明](#9-配置文件说明)
- [10. 输出结果说明](#10-输出结果说明)
- [11. 可选输入缺失时的行为](#11-可选输入缺失时的行为)
- [12. 重跑与断点续跑机制](#12-重跑与断点续跑机制)
- [13. 常见问题 / Troubleshooting](#13-常见问题--troubleshooting)
- [14. 开发状态与后续扩展](#14-开发状态与后续扩展)
- [15. 示例命令（可直接改路径运行）](#15-示例命令可直接改路径运行)

---

## 1. 项目背景与目标

本项目面向 SCCC 配对 tumor-normal WGS 数据，目标是提供一套可发表级扩展的分析框架：

1. 以 `markdup` 后 BAM 为起点（不从 FASTQ 重比对）。
2. 使用 ASCAT 输出 purity / ploidy / allele-specific CN / LOH 等 CNV 结果。
3. 使用 Manta + GRIDSS 双工具进行 SV 检测，并进行后处理、合并与注释。
4. 支持导入 HPV integration 断点并与宿主 CN/SV 做联动统计。
5. 在 SLURM 环境实现模块化、可重跑、可 array 并行、可阶段化控制。

---

## 2. 项目特点 / 核心功能

- **WGS 原生设计**（默认不提供 `wes` 模式）。
- **配对 tumor-normal 输入**（`pairbam.tsv`）。
- **样本级 array 并行**（soft_qc、ASCAT、SV calling、SV 注释、HPV 联动）。
- **依赖驱动 DAG**（非严格串行）。
- **Phase + Stage 双层控制**（粗粒度 + 细粒度）。
- **自动跳过可选模块**（无 `group.tsv` / 无 `hpv_breakpoints.tsv` 时自动跳过对应步骤）。
- **可重跑与断点续跑**（`.done` + 结果目录快照）。

---

## 3. 目录结构说明

```text
.
├── 01_submit_slurm_array.sh         # 主提交入口（phase + stage + DAG）
├── config/
│   └── config.yaml.example          # 配置模板（参考、资源、模块参数）
├── env/
│   └── environment.yaml             # 环境依赖示例
├── examples/
│   ├── pairbam.example.tsv
│   ├── group.example.tsv
│   └── hpv_breakpoints.example.tsv
├── logs/                            # Slurm 日志根目录
├── results/                         # 输出根目录
├── scripts/
│   ├── lib/common.sh                # 公共函数（目录初始化、done、快照）
│   ├── modules/*.sh                 # 各模块 wrapper
│   ├── python/*.py                  # 占位/接口脚本
│   ├── r/ascat_run.R                # ASCAT R 占位
│   ├── run_all.sh                   # 快速全流程提交包装
│   └── submit_module.sh             # 预留模块提交脚本
├── slurm/*.sbatch                   # 每模块 sbatch 模板
├── slurm/templates/                 # array/single 通用模板
└── workflow/Snakefile               # 迁移 Snakemake 的入口占位
```

---

## 4. 输入文件格式

### 4.1 `pairbam.tsv`（必需）

TAB 分隔，前三列必须固定：

| 列名 | 含义 | 约束 |
|---|---|---|
| sample_id | 样本 ID | 必须唯一 |
| tumor_bam | 肿瘤 BAM 绝对路径 | 文件必须存在 |
| normal_bam | 对照 BAM 绝对路径 | 文件必须存在且不同于 tumor_bam |

示例：

```tsv
sample_id	tumor_bam	normal_bam
TSDX001	/data/.../TSDX001.markdup.bam	/data/.../DSDX001.markdup.bam
```

> 要求 BAM 对应 `.bai` 可用；input_check 会检查可读性与兼容性。

### 4.2 `group.tsv`（可选）

第一列 `sample_id`，后续为任意分组变量（subtype/hpv_state/response 等）。

```tsv
sample_id	subtype	hpv_state	response
TSDX001	SCCC-A	HPV+	PR
```

### 4.3 `hpv_breakpoints.tsv`（可选）

支持无表头 raw 行（兼容外部断点导出文本），脚本会尝试解析 host `chr:pos`。

示例：

```text
ID=0 chr15:-46006922 HPV39REF|lcl|Human:-6786 SUPPORTING_PAIRS=10 SPLIT_READS=0
```

---

## 5. 参考资源文件

默认 hg38 资源：

- fasta：`/data/person/wup/public/liusy_files/reference_genomes/hg38/reference/Homo_sapiens_assembly38.fasta`
- ASCAT loci/alleles/GC 资源：统一放置在 `ascat_resources.root_dir`（正式环境建议：`/data/person/wup/public/liusy_files/reference_genomes/hg38/resources`）。
- accessible BED：`/data/person/wup/public/liusy_files/reference_genomes/hg38/resources/access.hg38.bed`

用途：

- fasta：SV/CNV 坐标体系、染色体命名检查。
- ASCAT loci/alleles/GC：用于固定位点计数、normal 杂合筛选、tumor BAF 计算、GC 矫正 logR。
- accessible BED：可访问区域过滤，避免低可比区域引入噪音。

---

## 6. 流程总体架构（Phase + Stage）

### Phase 是粗粒度控制
- `precheck`
- `sample`
- `cohort`
- `all`

### Stage 是细粒度控制
支持 `--start-stage` / `--end-stage` 在 phase 内截取窗口。

### 阶段定义

#### 阶段0：前置检查
- `input_check`（必须最先；默认不跳过）

#### 阶段1：样本级分析（array）

1A 可并行分叉：
- `soft_qc`
- `ascat_prepare`
- `sv_call_manta`
- `sv_call_gridss`

1B 依赖上游：
- `ascat_run`（依赖 `ascat_prepare`）
- `sv_postfilter`（依赖 `sv_call_manta + sv_call_gridss`）

1C 样本级整合：
- `sv_merge`（依赖 `sv_postfilter`）

1D 样本级注释：
- `sv_annotation`（依赖 `sv_merge + ascat_run`）

1E HPV 联动（可选）：
- `hpv_link`（依赖 `ascat_run + sv_annotation`；无 `hpv_breakpoints.tsv` 自动跳过）

#### 阶段2：队列级汇总（单次）
- `cohort_summary`（依赖 `ascat_run + sv_annotation`）
- `group_compare`（依赖 `cohort_summary`；若比较 HPV linked 指标可附加 `hpv_link`；无 group 文件自动跳过）
- `final_report`（依赖 `cohort_summary + soft_qc`；若存在 group/hpv 结果则纳入）

> 说明：本流程不再提供 `--enable-annotation`。`sv_annotation` 是主流程必需阶段。

---

## 7. 各模块详细说明

> 下表为“当前版本”的接口契约与输出约定，便于后续替换为真实算法实现。

### 7.1 input_check
- 输入：`pairbam.tsv`、`config.yaml`
- 输出：`results/input_check/input_check.report.tsv`
- 依赖：无
- 功能：表头/唯一性/绝对路径/BAM+BAI/可读性检查 + ASCAT 正式资源检查
  - 自动探测或读取 `ascat_resources.{loci_path,alleles_path,gc_path}`
  - 检查 loci/alleles/GC 文件存在性
  - 检查 BAM 染色体命名与 ASCAT 资源 chr 风格一致性
  - 在报告中写入资源识别结果（INFO 行）
- array：否

### 7.2 soft_qc
- 输入：每样本 tumor/normal BAM
- 输出：`results/soft_qc/<sample>/soft_qc.summary.tsv`
- 依赖：`input_check`
- 功能：基础 QC 汇总（不硬停）
  - reads/mapping/duplicate 指标来自 `samtools flagstat`
  - 覆盖度使用**抽样窗口估计**（`samtools bedcov`），默认主染色体随机窗口；避免 `samtools depth -a` 全基因组逐位点扫描
  - `summary.tsv` 额外记录 `coverage_method / n_windows / window_size` 便于追踪估计策略
- array：是

### 7.3 ascat_prepare
- 输入：BAM + 固定 ASCAT 资源（loci/alleles/GC）
- 输出：`baf.tsv`、`logr.tsv`、`ascat_prepare.run_info.tsv`
- 依赖：`input_check`
- 功能：按 nf-core/sarek 的 ASCAT/HTS 思路组织 ASCAT 预处理
  - 不再从 `reference.snp_vcf` 抽样位点；直接使用固定 loci + alleles + GC 资源
  - BAF：normal 先筛选高质量杂合位点，再用 tumor 在同位点的等位计数计算 BAF
  - 计数优先 `alleleCounter`；若缺失或失败，自动 fallback 到 `samtools mpileup`
  - logR：基于 tumor/normal depth ratio，并结合固定 GC 资源做线性 GC 矫正（不足时 fallback 到 raw ratio）
  - `run_info` 明确记录 loci_path / alleles_path / gc_path / BAF 方法 / logR 方法 / 是否 fallback
- array：是

### 7.4 ascat_run
- 输入：ascat_prepare 结果
- 输出：`purity_ploidy.tsv`、`segments.tsv`、`allele_specific_cn.tsv`、`loh.tsv`、`arm_level_cn.tsv`、`gene_level_cn.tsv`、`run_info.tsv`
- 依赖：`ascat_prepare`
- 功能：ASCAT 拟合与结果导出
- array：是

### 7.5 sv_call_manta
- 输入：配对 BAM
- 输出：`somaticSV.vcf.gz`、`manta.summary.tsv`
- 依赖：`input_check`
- 功能：Manta 体细胞 SV calling
- array：是

### 7.6 sv_call_gridss
- 输入：配对 BAM（**normal 在前 tumor 在后**）
- 输出：`gridss.somatic.vcf.gz`、`gridss.summary.tsv`
- 依赖：`input_check`
- 功能：GRIDSS SV calling
- array：是

### 7.7 sv_postfilter
- 输入：Manta + GRIDSS VCF
- 输出：`manta.postfilter.tsv`、`gridss.postfilter.tsv`
- 依赖：`sv_call_manta + sv_call_gridss`
- 功能：支持证据抽取、标准化
- array：是

### 7.8 sv_merge
- 输入：postfilter TSV
- 输出：`merged.sv.tsv`
- 依赖：`sv_postfilter`
- 功能：样本级 union merge + source 标记（MANTA/GRIDSS/BOTH）
- array：是

### 7.9 sv_annotation
- 输入：`merged.sv.tsv` + ascat 结果（+可选 HPV 位点）
- 输出：`annotated.sv.tsv`
- 依赖：`sv_merge + ascat_run`
- 功能：当前为**占位实现**，仅输出固定表头 + `NA` 结果，用于验证流程接口与依赖，不用于真实生物学注释结论
- array：是

### 7.10 hpv_link
- 输入：`hpv_breakpoints.tsv` + ascat + annotated SV
- 输出：`hpv_cn_link.tsv`、`hpv_sv_link.tsv`、`hpv_link.summary.tsv`
- 依赖：`ascat_run + sv_annotation`
- 功能：按窗口（10kb/100kb/1Mb）联动统计
- array：是

### 7.11 cohort_summary
- 输入：全样本 ASCAT + SV annotation
- 输出：`cohort.purity_ploidy.tsv`、`cohort.arm_level_cn.matrix.tsv`、`cohort.gene_level_cn.matrix.tsv`、`cohort.loh.matrix.tsv`、`cohort.sv_burden.tsv`、`cohort.sv_event_table.tsv`
- 依赖：`ascat_run + sv_annotation`
- 功能：队列级汇总与矩阵化
- array：否

### 7.12 group_compare
- 输入：`group.tsv` + cohort summary
- 输出：按分组列导出的 summary / plot 占位结果
- 依赖：`cohort_summary`（可选附加 `hpv_link`）
- 功能：分组汇总比较
- array：否

### 7.13 final_report
- 输入：cohort/group/hpv/soft_qc 结果
- 输出：`results/final_report/final_report.md`
- 依赖：`cohort_summary + soft_qc`（可纳入 group/hpv）
- 功能：项目级总览报告
- array：否

---



### Group 驱动的 cohort 行为

- `group.tsv` 存在：允许 `cohort_summary`，若存在有效分组列（除 `sample_id` 外）则运行 `group_compare`，`final_report` 可生成完整 cohort 报告。
- `group.tsv` 不存在：
  - `--phase cohort`：直接报错退出。
  - `--phase all`：自动跳过 `cohort_summary/group_compare/final_report`，仅运行 precheck + sample，并给 warning。
  - `--phase sample`：不要求 `group.tsv`。
- 若 `group.tsv` 仅有 `sample_id` 而无有效分组列：跳过 `group_compare`，但 `cohort_summary` 与 `final_report` 仍执行。

## 8. SLURM 运行方式

### 8.1 提交全流程

```bash
bash 01_submit_slurm_array.sh \
  --pairs input/sample_pairs.tsv \
  --group input/group.tsv \
  --hpv-breakpoints input/hpv_breakpoints.tsv \
  --results-dir results_run_20260421 \
  --phase all \
  --max-parallel 2 \
  --config config/config.yaml.example
```

### 8.2 只跑某个 phase

```bash
bash 01_submit_slurm_array.sh \
  --pairs input/sample_pairs.tsv \
  --results-dir results_run_20260421 \
  --phase sample \
  --config config/config.yaml.example
```

### 8.3 只跑 stage 窗口

```bash
bash 01_submit_slurm_array.sh \
  --pairs input/sample_pairs.tsv \
  --results-dir results_run_20260421 \
  --phase all \
  --start-stage sv_postfilter \
  --end-stage final_report \
  --config config/config.yaml.example
```

### 8.4 控制 array 并发

```bash
# 全局并发
--max-parallel 2

# 单 stage 覆盖（可重复）
--max-parallel-stage sv_call_gridss=1 \
--max-parallel-stage ascat_run=1
```

### 8.5 跳过 input_check（可选，不推荐）

```bash
--skip-input-check 1
```

默认 `--skip-input-check 0`。

### 8.6 只重跑某样本 / 某阶段

- 某样本：删除对应样本目录 `.done` 后，重投该阶段 array（可用 `--start-stage`/`--end-stage` 限定）。
- 某阶段：删除阶段输出 `.done` 或阶段目录中的目标结果后重投。

---

## 9. 配置文件说明

`config/config.yaml.example` 关键内容：

- `reference`：hg38 fasta / access BED（`snp_vcf` 可保留给其他模块，不再作为 ASCAT 正式主输入）。
- `ascat_resources`：ASCAT 正式资源根目录与 loci/alleles/gc 路径、chr 风格、alleleCounter 可执行文件。
- `parallel.max_parallel_default`：默认并发。
- `modules.<stage>.slurm.*`：分模块资源策略（partition/cpu/mem/time/log）。
- `qc.*`：`soft_qc` 覆盖度估计策略参数（`coverage_method`、`n_windows`、`window_size`、`include_chroms`），用于平衡速度与稳健性。
- `ascat_prepare.*`：ASCAT 预处理参数（`max_sites`、`min_normal_depth`、`min_tumor_depth`、`include_chroms`、normal het BAF 区间）。
- `hpv_link.windows_bp`：联动窗口。

### CPU 分区策略建议

- `cpu1`：常规模块（soft_qc、ascat_prepare/run、manta、postfilter、summary）。
- `cpu2`：重计算模块（gridss）。

### 当前推荐资源配额（与示例配置一致）

- `input_check`: `cpu1`, 1 CPU, 4G
- `soft_qc`: `cpu1`, 4 CPU, 16G
- `ascat_prepare`: `cpu1`, 8 CPU, 32G
- `ascat_run`: `cpu1`, 4 CPU, 16G
- `sv_call_manta`: `cpu1`, 8 CPU, 32G
- `sv_call_gridss`: `cpu2`, 16 CPU, 64G
- `sv_postfilter`: `cpu1`, 2 CPU, 8G
- `sv_merge`: `cpu1`, 2 CPU, 8G
- `sv_annotation`: `cpu1`, 2 CPU, 8G
- `hpv_link`: `cpu1`, 2 CPU, 8G
- `cohort_summary`: `cpu1`, 4 CPU, 16G
- `group_compare`: `cpu1`, 2 CPU, 8G
- `final_report`: `cpu1`, 1 CPU, 4G, `02:00:00`

---

## 10. 输出结果说明

### 样本级 ASCAT
- `purity_ploidy.tsv`：样本纯度和倍性
- `segments.tsv`：CN 分段
- `allele_specific_cn.tsv`：等位基因 CN
- `loh.tsv`：LOH 区域
- `arm_level_cn.tsv` / `gene_level_cn.tsv`
- `run_info.tsv`：收敛/运行状态

### 样本级 SV
- `somaticSV.vcf.gz`（Manta 原始）
- `gridss.somatic.vcf.gz`（GRIDSS 原始）
- `*.postfilter.tsv`（标准化后）
- `merged.sv.tsv`（样本级并集）
- `annotated.sv.tsv`（样本级注释）

### 队列级
- `cohort.purity_ploidy.tsv`
- `cohort.arm_level_cn.matrix.tsv`
- `cohort.gene_level_cn.matrix.tsv`
- `cohort.loh.matrix.tsv`
- `cohort.sv_burden.tsv`
- `cohort.sv_event_table.tsv`

### HPV 联动
- `hpv_cn_link.tsv`
- `hpv_sv_link.tsv`
- `hpv_link.summary.tsv`

---

## 11. 可选输入缺失时的行为

- 缺少 `group.tsv`：
  - `--phase all`：自动跳过整个 cohort 阶段（`cohort_summary/group_compare/final_report`）。
  - `--phase cohort`：报错退出。
  - `--phase sample`：不受影响。
- 缺少 `hpv_breakpoints.tsv`：`hpv_link` 自动跳过（不影响主 CNV/SV）。

---

## 12. 重跑与断点续跑机制

- 每个模块输出目录写入 `.done`。
- `.done` 记录 `status/timestamp/stage/sample_id`。
- 每次提交会将当前 `config` 与 `pairbam` 快照至 `results_dir/00.config/`。
- 若用 `--start-stage/--end-stage` 从中间重跑，提交脚本会检查被截断上游的 `.done` 是否存在（样本级检查到每个 sample）；缺失时直接报错，避免“误跳过上游导致下游失败”。
- 若担心旧结果复用：
  - 推荐新建 `results_dir`；或
  - 删除目标模块 `.done` 强制重跑。

---

## 13. 常见问题 / Troubleshooting

1. **chr 命名不一致（chr1 vs 1）**  
   优先在 `input_check` 阶段修正参考与 BAM/VCF 一致性。

2. **BAM/BAI 缺失或不可读**  
   `input_check` 会失败并阻断后续；先修复路径与索引。

3. **input_check 失败后流程不启动**  
   这是预期行为；后续均依赖其 `afterok`。

4. **Slurm Pending 很久**  
   调低 `--max-parallel` 或对重模块使用 `--max-parallel-stage` 限流。

5. **array 日志怎么看**  
   查看 `logs/<module>/%A_%a.out` 和 `%A_%a.err`。

---

## 14. 开发状态与后续扩展

### 当前版本
- 已完成 phase+stage 控制、DAG 提交、array 框架、I/O 契约、可选模块自动跳过、断点续跑机制。
- `sv_annotation` 仍为占位模块，因此当前 pipeline 可用于**提交框架与可运行性测试**，但不能用于验证真实注释生物学结果。

### 计划扩展
- 将占位实现替换为真实生产命令（ASCAT/Manta/GRIDSS 参数细化）。
- 增强统计与可视化。
- 接入 `snv_indel_hook`（SNV/Indel + driver 注释）并与 CN/SV 联合解释。

---

## 15. 示例命令（可直接改路径运行）

### 示例1：全流程
```bash
bash 01_submit_slurm_array.sh \
  --pairs input/sample_pairs.tsv \
  --group input/group.tsv \
  --hpv-breakpoints input/hpv_breakpoints.tsv \
  --results-dir results/run_all \
  --phase all \
  --max-parallel 2 \
  --config config/config.yaml.example
```

### 示例2：只做 precheck
```bash
bash 01_submit_slurm_array.sh \
  --pairs input/sample_pairs.tsv \
  --results-dir results/precheck \
  --phase precheck \
  --config config/config.yaml.example
```

### 示例3：只跑样本级 phase
```bash
bash 01_submit_slurm_array.sh \
  --pairs input/sample_pairs.tsv \
  --results-dir results/sample_only \
  --phase sample \
  --max-parallel 2 \
  --config config/config.yaml.example
```

### 示例4：只跑 sv_postfilter 到 final_report
```bash
bash 01_submit_slurm_array.sh \
  --pairs input/sample_pairs.tsv \
  --group input/group.tsv \
  --results-dir results/rerun_tail \
  --phase all \
  --start-stage sv_postfilter \
  --end-stage final_report \
  --max-parallel 2 \
  --config config/config.yaml.example
```

### 示例5：限制 GRIDSS 并发
```bash
bash 01_submit_slurm_array.sh \
  --pairs input/sample_pairs.tsv \
  --results-dir results/gridss_limited \
  --phase all \
  --max-parallel 3 \
  --max-parallel-stage sv_call_gridss=1 \
  --config config/config.yaml.example
```

### 示例6：已有检查结果时跳过 input_check（谨慎）
```bash
bash 01_submit_slurm_array.sh \
  --pairs input/sample_pairs.tsv \
  --results-dir results/skip_precheck \
  --phase sample \
  --skip-input-check 1 \
  --config config/config.yaml.example
```
