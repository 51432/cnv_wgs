#!/usr/bin/env python3
import argparse
import math
import random
import re
import subprocess
import tempfile
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


def parse_qc_config(config_path):
    defaults = {
        "coverage_method": "sampled_windows",
        "n_windows": 200,
        "window_size": 100000,
        "include_chroms": [
            *(f"chr{i}" for i in range(1, 23)),
            "chrX",
            "chrY",
        ],
    }

    lines = Path(config_path).read_text().splitlines()
    in_qc = False

    for raw in lines:
        line = raw.rstrip()
        stripped = line.strip()
        if not stripped or stripped.startswith("#"):
            continue

        if not in_qc:
            if re.match(r"^qc\s*:\s*$", stripped):
                in_qc = True
            continue

        if raw and not raw.startswith(" ") and not raw.startswith("\t"):
            break

        m = re.match(r"^\s*(coverage_method|n_windows|window_size)\s*:\s*(.+?)\s*$", raw)
        if m:
            key, val = m.group(1), m.group(2).strip().strip('"').strip("'")
            if key in {"n_windows", "window_size"}:
                try:
                    defaults[key] = int(val)
                except ValueError:
                    pass
            else:
                defaults[key] = val
            continue

        m_chroms = re.match(r"^\s*include_chroms\s*:\s*\[(.*)\]\s*$", raw)
        if m_chroms:
            body = m_chroms.group(1).strip()
            if body:
                chroms = [x.strip().strip('"').strip("'") for x in body.split(",") if x.strip()]
                if chroms:
                    defaults["include_chroms"] = chroms

    defaults["n_windows"] = max(1, int(defaults["n_windows"]))
    defaults["window_size"] = max(1000, int(defaults["window_size"]))
    return defaults


def norm_chrom(chrom):
    c = chrom.lower()
    if c.startswith("chr"):
        c = c[3:]
    return c


def read_idxstats(bam_path):
    out = run_cmd(["samtools", "idxstats", bam_path])
    contigs = []
    for line in out.splitlines():
        p = line.split("\t")
        if len(p) < 2:
            continue
        chrom = p[0]
        if chrom == "*":
            continue
        try:
            length = int(p[1])
        except ValueError:
            continue
        if length > 0:
            contigs.append((chrom, length))
    return contigs


def sample_windows_from_contigs(contigs, include_chroms, n_windows, window_size, seed):
    include_norm = {norm_chrom(c) for c in include_chroms}
    candidates = [(c, l) for c, l in contigs if norm_chrom(c) in include_norm and l >= window_size]
    if not candidates:
        raise RuntimeError("No eligible contigs for coverage sampling; check include_chroms/window_size")

    weights = [l for _, l in candidates]
    rng = random.Random(seed)

    windows = []
    for _ in range(n_windows):
        chrom, length = rng.choices(candidates, weights=weights, k=1)[0]
        max_start = length - window_size
        start0 = rng.randint(0, max_start) if max_start > 0 else 0
        end = start0 + window_size
        windows.append((chrom, start0, end))
    return windows


def mean_depth_sampled_windows(tumor_bam, normal_bam, include_chroms, n_windows, window_size, seed):
    contigs = read_idxstats(tumor_bam)
    windows = sample_windows_from_contigs(contigs, include_chroms, n_windows, window_size, seed)

    tumor_sum = 0.0
    normal_sum = 0.0
    with tempfile.NamedTemporaryFile("w", suffix=".bed", delete=True) as bedf:
        for chrom, s0, e1 in windows:
            bedf.write(f"{chrom}\t{s0}\t{e1}\n")
        bedf.flush()

        out = run_cmd(["samtools", "bedcov", bedf.name, tumor_bam, normal_bam])

    n_rows = 0
    for line in out.splitlines():
        p = line.split("\t")
        if len(p) < 5:
            continue
        try:
            tumor_sum += float(p[3])
            normal_sum += float(p[4])
            n_rows += 1
        except ValueError:
            continue

    if n_rows == 0:
        return math.nan, math.nan

    total_bases = n_rows * window_size
    return tumor_sum / total_bases, normal_sum / total_bases


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
    ap.add_argument('--config', required=True)
    ap.add_argument('--outdir', required=True)
    args = ap.parse_args()

    qc_cfg = parse_qc_config(args.config)
    coverage_method = qc_cfg["coverage_method"]
    n_windows = qc_cfg["n_windows"]
    window_size = qc_cfg["window_size"]
    include_chroms = qc_cfg["include_chroms"]

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    t_total, t_mapped, t_dup = parse_flagstat(run_cmd(["samtools", "flagstat", args.tumor_bam]))
    n_total, n_mapped, n_dup = parse_flagstat(run_cmd(["samtools", "flagstat", args.normal_bam]))

    if coverage_method != "sampled_windows":
        raise RuntimeError(f"Unsupported qc.coverage_method: {coverage_method}")

    seed = f"{args.sample_id}:{n_windows}:{window_size}"
    mean_cov_t, mean_cov_n = mean_depth_sampled_windows(
        args.tumor_bam,
        args.normal_bam,
        include_chroms,
        n_windows,
        window_size,
        seed,
    )

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
        "coverage_method",
        "n_windows",
        "window_size",
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
        coverage_method,
        str(n_windows),
        str(window_size),
    ]

    (outdir / 'soft_qc.summary.tsv').write_text("\t".join(header) + "\n" + "\t".join(row) + "\n")


if __name__ == '__main__':
    main()
