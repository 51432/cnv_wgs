# SCCC 配对 tumor-normal WGS（hg38）CNV/SV 学术型模块化 Pipeline（SLURM）

> 第一版实现重点：**目录结构 + 配置模板 + 模块脚本骨架 + sbatch 模板 + 重跑机制**。  
> 设计目标：发表级 CNV/SV 分析，兼容后续 HPV integration 联动。

## 1. 项目结构

```text
.
├── config/
│   └── config.yaml.example
├── env/
│   └── environment.yaml
├── examples/
│   ├── pairbam.example.tsv
│   ├── group.example.tsv
│   └── hpv_breakpoints.example.tsv
├── logs/
├── profiles/
│   └── slurm/config.yaml
├── results/
├── scripts/
│   ├── lib/common.sh
│   ├── modules/
│   │   ├── input_check.sh
│   │   ├── soft_qc.sh
│   │   ├── ascat_prepare.sh
│   │   ├── ascat_run.sh
│   │   ├── sv_call_manta.sh
│   │   ├── sv_call_gridss.sh
│   │   ├── sv_postfilter.sh
│   │   ├── sv_merge.sh
│   │   ├── sv_annotation.sh
│   │   ├── cohort_summary.sh
│   │   ├── group_compare.sh
│   │   ├── hpv_link.sh
│   │   ├── snv_indel_hook.sh
│   │   └── final_report.sh
│   ├── python/*.py
│   ├── r/ascat_run.R
│   ├── run_all.sh
│   └── submit_module.sh
├── slurm/
│   ├── *.sbatch
│   └── templates/
└── workflow/Snakefile
```

---

## 2. 输入规范

### 2.1 `pairbam.tsv`（必需）
- TAB 分隔，前三列必须固定：
  - `sample_id`
  - `tumor_bam`
  - `normal_bam`
- 要求：
  - BAM 为绝对路径
  - 对应 `.bai` 必须存在
  - `sample_id` 唯一
  - tumor/normal 不能同文件

### 2.2 `group.tsv`（可选）
- 第一列必须是 `sample_id`
- 其余列可自由扩展（如 `subtype`, `hpv_state`, `response`）

### 2.3 `hpv_breakpoints.tsv`（可选）
- 当前支持 **无表头 raw 行格式**（示例已给）
- `scripts/python/hpv_link.py` 会解析 host 侧 `chr:pos`。

---

## 3. 固定 hg38 参考
在 `config/config.yaml.example` 已预置：
- fasta: `/data/person/wup/public/liusy_files/reference_genomes/hg38/reference/Homo_sapiens_assembly38.fasta`
- snp vcf: `/data/person/wup/public/liusy_files/reference_genomes/hg38/resources/1000G_phase1.snps.high_confidence.hg38.vcf.gz`
- accessible bed: `/data/person/wup/public/liusy_files/reference_genomes/hg38/resources/access.hg38.bed`

---

## 4. 模块与 I/O（第一版骨架）

| 模块 | 类型 | 关键输入 | 关键输出 | 成功判定 |
|---|---|---|---|---|
| input_check | 单任务 | pairbam, config | input_check.report.tsv | `.done` |
| soft_qc | array | pairbam, bam | soft_qc.summary.tsv | `sample/.done` |
| ascat_prepare | array | pairbam, snp vcf, access bed | baf.tsv, logr.tsv | `sample/.done` |
| ascat_run | array | ascat_prepare 输出 | purity_ploidy.tsv, segments.tsv, loh.tsv... | `sample/.done` |
| sv_call_manta | array | pairbam | somaticSV.vcf.gz, manta.summary.tsv | `sample/.done` |
| sv_call_gridss | array | pairbam | gridss.somatic.vcf.gz, gridss.summary.tsv | `sample/.done` |
| sv_postfilter | array | manta/gridss VCF | manta.postfilter.tsv, gridss.postfilter.tsv | `sample/.done` |
| sv_merge | 单任务 | postfilter | merged.sv.tsv | `.done` |
| sv_annotation | 单任务 | merged SV + ascat + hpv | annotated.sv.tsv | `.done` |
| cohort_summary | 单任务 | ascat + sv annotation | cohort.* 矩阵与表 | `.done` |
| group_compare | 单任务 | group.tsv + cohort 输出 | 分组 summary/图占位 | `.done` |
| hpv_link | 单任务 | hpv_breakpoints + ascat/sv | per-sample + cohort summary | `.done` |
| snv_indel_hook | 单任务 | 预留 | README.txt | `.done` |
| final_report | 单任务 | 项目结果目录 | final_report.md | `.done` |

> `skip_done=true` 时自动跳过已完成样本/模块。

---

## 5. SLURM 与资源策略

- 轻量模块默认 `cpu1`（2-8 CPU, 8-32G）
- 中等模块默认 `cpu1`（4-16 CPU, 32-64G）
- GRIDSS 默认 `cpu2`（16-32 CPU, 128-256G）
- 全部模块参数均可在 `config.yaml` 独立覆盖：
  - `partition`
  - `cpus_per_task`
  - `mem`
  - `time`
  - `job_name`
  - `stdout/stderr`

### 5.1 日志路径
自动使用并创建：
- `logs/input_check/`
- `logs/soft_qc/`
- `logs/ascat_prepare/`
- `logs/ascat/`
- `logs/manta/`
- `logs/gridss/`
- `logs/sv_merge/`
- `logs/report/`
- 及其余模块目录

array 日志格式默认：`%A_%a.out` / `%A_%a.err`

---

## 6. 执行顺序（推荐）

1. input_check
2. soft_qc（array）
3. ascat_prepare（array）
4. ascat_run（array）
5. sv_call_manta（array）
6. sv_call_gridss（array，normal 在前 tumor 在后）
7. sv_postfilter（array）
8. sv_merge
9. sv_annotation
10. cohort_summary
11. group_compare（若提供 group.tsv）
12. hpv_link（若提供 hpv_breakpoints.tsv）
13. final_report

---

## 7. 提交方式示例

### 7.1 单模块提交
```bash
sbatch slurm/input_check.sbatch examples/pairbam.example.tsv config/config.yaml.example results/input_check
```

### 7.2 Array 模块（以 soft_qc 为例）
```bash
N=$(awk 'NR>1{n++} END{print n+0}' examples/pairbam.example.tsv)
sbatch --array=1-${N}%2 slurm/soft_qc.sbatch examples/pairbam.example.tsv results/soft_qc
```

### 7.3 一键骨架流程提交（无 dependency 编排）
```bash
bash scripts/run_all.sh config/config.yaml.example examples/pairbam.example.tsv results
```

---

## 8. 重跑策略

### 8.1 只重跑某一步
删除模块级 `.done`（或模块输出目录），再重新 `sbatch`。

### 8.2 只重跑某个样本
删除样本输出目录下 `.done`：
```bash
rm results/ascat/TSDX001/.done
```
然后用 array 仅投单 task：
```bash
sbatch --array=1-1 slurm/ascat_run.sbatch examples/pairbam.example.tsv config/config.yaml.example results/ascat_prepare results/ascat
```

### 8.3 全量重跑
将 `runtime.skip_done` 设为 `false`，或清空 `results/`。

---

## 9. Snakemake 迁移位点
- 已提供 `workflow/Snakefile` 与 `profiles/slurm/config.yaml` 占位。
- 当前建议先用模块化 `sbatch` 稳定生产。
- 后续可将 `scripts/modules/*` 封装为 Snakemake rule shell 命令。

---

## 10. 第一版说明
当前为可运行骨架，重点保证：
- SLURM 友好（array + 独立 sbatch）
- 模块化与重跑
- 输入输出契约清晰
- compatibility check 独立且可阻断
- HPV 联动接口可扩展

后续可在各 `scripts/python/*.py` 与 `scripts/r/ascat_run.R` 内逐步替换为真实算法调用。
