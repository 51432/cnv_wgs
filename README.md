# SCCC Tumor-Normal WGS（hg38）CNV/SV Pipeline（SLURM）

> 默认面向 **WGS**，不提供 `mode=wes`。

## 核心入口

```bash
bash 01_submit_slurm_array.sh --help
```

## 依赖关系（DAG）

### 强制规则
1. `input_check` 永远先提交。  
2. 后续所有任务都带 `afterok:<input_check_jobid>` 依赖。  

### 样本级 array
- `soft_qc`
- `ascat_prepare`
- `ascat_run`（依赖 `ascat_prepare`）
- `sv_call_manta`
- `sv_call_gridss`（normal 在前 tumor 在后）
- `sv_postfilter`（依赖 `sv_call_manta` + `sv_call_gridss`）
- `sv_merge`（按样本）
- `sv_annotation`（依赖 `sv_merge` + `ascat_run`）
- `hpv_link`（可选，依赖 `ascat_run` + `sv_annotation`）

### 队列级单次任务
- `input_check`
- `cohort_summary`（依赖 `ascat_run` + `sv_annotation`）
- `group_compare`（依赖 `cohort_summary`；若比较 HPV linked 指标可附加 `hpv_link`）
- `final_report`（依赖 `cohort_summary` + `group_compare` + `soft_qc`）

## 结果目录与 `.done` 机制

- `--results-dir` 指定结果根目录，只在目录不存在时创建。
- 运行时会把当前配置快照到：`results_dir/00.config/`。
- `.done` 文件记录 `status/timestamp/stage/sample_id`，用于跳过机制并减少误复用风险。

## 运行示例

### 全流程（推荐）
```bash
bash 01_submit_slurm_array.sh \
  --pairs input/sample_pairs.tsv \
  --group input/group.tsv \
  --hpv-breakpoints input/hpv_breakpoints.tsv \
  --results-dir results_run_20260421 \
  --max-parallel 2 \
  --max-parallel-stage ascat_run=1 \
  --max-parallel-stage sv_call_gridss=1 \
  --start-stage input_check \
  --end-stage final_report \
  --enable-annotation 1 \
  --config config/config.yaml.example
```

### 只重跑部分阶段
```bash
bash 01_submit_slurm_array.sh \
  --pairs input/sample_pairs.tsv \
  --results-dir results_run_20260421 \
  --start-stage sv_postfilter \
  --end-stage final_report \
  --max-parallel 2 \
  --config config/config.yaml.example
```

## 关键文件
- `01_submit_slurm_array.sh`
- `scripts/lib/common.sh`
- `scripts/modules/*.sh`
- `slurm/*.sbatch`
- `config/config.yaml.example`
