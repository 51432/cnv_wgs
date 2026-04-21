# SCCC 配对 tumor-normal WGS（hg38）CNV/SV 学术型模块化 Pipeline（SLURM）

> 第一版实现重点：模块化、可 array 并行、可阶段控制提交、便于迁移到 Snakemake。

## 1) 项目结构

```text
.
├── 01_submit_slurm_array.sh
├── config/config.yaml.example
├── examples/
├── scripts/lib/common.sh
├── scripts/modules/*.sh
├── scripts/python/*.py
├── scripts/r/ascat_run.R
├── slurm/*.sbatch
└── workflow/Snakefile
```

## 2) 输入

### `pairbam.tsv`（必需）
前三列固定：`sample_id tumor_bam normal_bam`。

### `group.tsv`（可选）
第一列 `sample_id`，其余分组变量自定义。

### `hpv_breakpoints.tsv`（可选）
支持无表头 raw 行，例如：
`ID=0 chr15:-46006922 HPV39REF|lcl|Human:-6786 ...`

## 3) 输出目录与“避免覆盖”策略

- 提交脚本使用 `--results-dir` 指定结果根目录。
- 仅在目录不存在时创建；已存在目录会保留并提示 `[keep]`，不会清空历史结果。
- 模块依赖 `.done` 跳过机制（默认 `SKIP_DONE=true`），避免重复覆盖已完成样本。

## 4) 模块列表与执行顺序

`input_check -> soft_qc -> ascat_prepare -> ascat_run -> sv_call_manta -> sv_call_gridss -> sv_postfilter -> sv_merge -> sv_annotation -> cohort_summary -> group_compare -> hpv_link -> final_report`

其中 array 模块：
- soft_qc
- ascat_prepare
- ascat_run
- sv_call_manta
- sv_call_gridss
- sv_postfilter

## 5) 一键提交脚本（新增）

新增 `01_submit_slurm_array.sh`，支持参数化阶段控制。

### 示例（与你给出的风格一致）

```bash
bash 01_submit_slurm_array.sh \
  --pipeline phase2 \
  --pairs input/sample_pairs.tsv \
  --mode wes \
  --max-parallel 1 \
  --end-stage filter \
  --enable-contamination 1 \
  --enable-orientation 1 \
  --enable-annotation 0
```

> 说明：`filter` 会映射到 `sv_postfilter`。

### 常用参数

- `--pipeline all|phase1|phase2|phase3`
- `--pairs <pairbam.tsv>`
- `--group <group.tsv>`
- `--hpv-breakpoints <hpv_breakpoints.tsv>`
- `--results-dir <dir>`
- `--max-parallel <N>`
- `--start-stage <stage>`
- `--end-stage <stage>`
- `--enable-annotation 0|1`
- `--enable-contamination 0|1`（第一版先透传记录）
- `--enable-orientation 0|1`（第一版先透传记录）

## 6) SLURM 资源默认策略

- 轻量模块：默认 `cpu1`，2-8 CPU，8-32G
- 中等模块：默认 `cpu1`，4-16 CPU，32-64G
- 重计算（GRIDSS）：默认 `cpu2`，16-32 CPU，128-256G
- 所有模块参数可在 `config/config.yaml.example` 单独覆盖。

## 7) 重跑策略

- 重跑某一步：删除该模块 `.done` 后重新提交。
- 重跑某样本：删除样本目录 `.done` 并提交相应 array task。
- 全量重跑：使用新 `results-dir`（推荐）或清理旧 `.done`。

## 8) Snakemake 迁移

保留 `workflow/Snakefile` 和 `profiles/slurm/config.yaml` 作为下一步迁移入口；当前主生产方式建议继续 `sbatch + array`。
