#!/usr/bin/env python3
import argparse
import math
import re
import subprocess
from pathlib import Path


def run_cmd(cmd):
    return subprocess.check_output(cmd, text=True, stderr=subprocess.STDOUT)


def parse_flagstat(flagstat_text):
    total_reads = None
    mapped_rate = math.nan
    duplicate_rate = math.nan
    dup_reads = None

    for line in flagstat_text.splitlines():
        m_total = re.match(r"^(\d+)\s+\+\s+\d+\s+in total", line)
        if m_total:
            total_reads = int(m_total.group(1))
            continue

        m_mapped = re.match(r"^(\d+)\s+\+\s+\d+\s+mapped\s+\(([-0-9.]+)%", line)
        if m_mapped:
            mapped_rate = float(m_mapped.group(2)) / 100.0
            continue

        m_dup = re.match(r"^(\d+)\s+\+\s+\d+\s+duplicates", line)
        if m_dup:
            dup_reads = int(m_dup.group(1))

    if total_reads and dup_reads is not None:
        duplicate_rate = dup_reads / max(total_reads, 1)

    return total_reads, mapped_rate, duplicate_rate


def mean_depth_from_bam(bam_path):
    proc = subprocess.Popen(
        ["samtools", "depth", "-a", bam_path],
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )
    assert proc.stdout is not None

    n = 0
    depth_sum = 0
    for line in proc.stdout:
        parts = line.rstrip("\n").split("\t")
        if len(parts) < 3:
            continue
        try:
            depth_sum += int(parts[2])
            n += 1
        except ValueError:
            continue

    _, stderr = proc.communicate()
    if proc.returncode != 0:
        raise RuntimeError(f"samtools depth failed for {bam_path}: {stderr.strip()}")

    if n == 0:
        return math.nan
    return depth_sum / n


def fmt_num(x, digits=6):
    if x is None or (isinstance(x, float) and math.isnan(x)):
        return "NA"
    if isinstance(x, int):
        return str(x)
    return f"{x:.{digits}f}"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--sample-id', required=True)
    ap.add_argument('--tumor-bam', required=True)
    ap.add_argument('--normal-bam', required=True)
    ap.add_argument('--outdir', required=True)
    args = ap.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    t_total, t_mapped, t_dup = parse_flagstat(run_cmd(["samtools", "flagstat", args.tumor_bam]))
    n_total, n_mapped, n_dup = parse_flagstat(run_cmd(["samtools", "flagstat", args.normal_bam]))
    mean_cov_t = mean_depth_from_bam(args.tumor_bam)
    mean_cov_n = mean_depth_from_bam(args.normal_bam)

    tn_ratio = mean_cov_t / mean_cov_n if mean_cov_n and not math.isnan(mean_cov_n) and mean_cov_n > 0 else math.nan

    qc_flag = "PASS"
    if (
        math.isnan(t_mapped)
        or math.isnan(n_mapped)
        or t_mapped < 0.90
        or n_mapped < 0.90
        or (not math.isnan(t_dup) and t_dup > 0.80)
        or (not math.isnan(n_dup) and n_dup > 0.80)
        or (not math.isnan(tn_ratio) and (tn_ratio < 0.5 or tn_ratio > 2.0))
    ):
        qc_flag = "REVIEW"

    header = [
        "sample_id",
        "total_reads_tumor",
        "mapped_rate_tumor",
        "duplicate_rate_tumor",
        "total_reads_normal",
        "mapped_rate_normal",
        "duplicate_rate_normal",
        "mean_cov_tumor",
        "mean_cov_normal",
        "tn_cov_ratio",
        "qc_flag",
    ]
    row = [
        args.sample_id,
        fmt_num(t_total),
        fmt_num(t_mapped),
        fmt_num(t_dup),
        fmt_num(n_total),
        fmt_num(n_mapped),
        fmt_num(n_dup),
        fmt_num(mean_cov_t),
        fmt_num(mean_cov_n),
        fmt_num(tn_ratio),
        qc_flag,
    ]

    (outdir / 'soft_qc.summary.tsv').write_text("\t".join(header) + "\n" + "\t".join(row) + "\n")


if __name__ == '__main__':
    main()
