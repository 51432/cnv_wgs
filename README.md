# SCCC Tumor-Normal WGS（hg38）CNV/SV Pipeline（SLURM）

## 核心提交入口

```bash
bash 01_submit_slurm_array.sh --help
```

支持阶段化、依赖驱动和 array 并行。

---

## 1) 执行拓扑（依赖关系驱动，不是严格串行）

### 样本级 array 步骤
- `soft_qc`
- `ascat_prepare`
- `ascat_run`（依赖 `ascat_prepare`）
- `sv_call_manta`
- `sv_call_gridss`（normal 在前，tumor 在后）
- `sv_postfilter`（依赖 `sv_call_manta` + `sv_call_gridss`）
- `sv_merge`（按样本 merge）
- `sv_annotation`（按样本注释，依赖 `sv_merge` + `ascat_run`）
- `hpv_link`（按样本联动，依赖 `ascat_run` + `sv_annotation`，需提供 hpv 文件）

### 队列级单次任务
- `input_check`
- `cohort_summary`（依赖样本级结果完成）
- `group_compare`（依赖 `cohort_summary`；如比较 HPV-linked 指标，也依赖 `hpv_link`）
- `final_report`

### 可并行分叉
- 在 `input_check` 后，以下分支可并行：
  - `soft_qc`
  - `ascat_prepare -> ascat_run`
  - `sv_call_manta`
  - `sv_call_gridss`

### 必须等待依赖
- `ascat_run` 必须等 `ascat_prepare`
- `sv_postfilter` 必须等两个 caller
- `sv_merge` 必须等 `sv_postfilter`
- `sv_annotation` 必须等 `sv_merge` 和 `ascat_run`
- `hpv_link` 必须等 `ascat_run` 和 `sv_annotation`
- `cohort_summary` 等样本级完成
- `group_compare` 等 `cohort_summary`（必要时再等 `hpv_link`）
- `final_report` 等汇总完成

---

## 2) 结果目录策略

`--results-dir` 指定输出根目录。目录仅在不存在时创建；已有目录不会被清空，避免覆盖历史结果。模块依赖 `.done` 跳过已完成样本。

---

## 3) 运行示例（你可直接改路径后执行）

### 全流程（推荐）
```bash
bash 01_submit_slurm_array.sh \
  --pipeline all \
  --pairs input/sample_pairs.tsv \
  --group input/group.tsv \
  --hpv-breakpoints input/hpv_breakpoints.tsv \
  --mode wgs \
  --results-dir results_run_20260421 \
  --max-parallel 2 \
  --start-stage input_check \
  --end-stage final_report \
  --enable-contamination 1 \
  --enable-orientation 1 \
  --enable-annotation 1 \
  --config config/config.yaml.example
```

### 仅跑 phase2（示例风格）
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

> `filter` 会映射到 `sv_postfilter`。

### 简化入口
```bash
bash scripts/run_all.sh config/config.yaml.example input/sample_pairs.tsv results_run_20260421
```

---

## 4) 关键文件
- 主提交脚本：`01_submit_slurm_array.sh`
- 公共函数：`scripts/lib/common.sh`
- 模块脚本：`scripts/modules/*.sh`
- sbatch 模板：`slurm/*.sbatch`
- 示例输入：`examples/*.tsv`
