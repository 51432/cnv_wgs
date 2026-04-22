#!/usr/bin/env python3
import argparse
import csv
import gzip
from pathlib import Path
import subprocess
import sys

import yaml


def check_bam_readable(bam: Path) -> bool:
    cmd = ["samtools", "quickcheck", str(bam)]
    r = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    return r.returncode == 0


def detect_bam_chr_style(bam: Path) -> str:
    cmd = ["samtools", "idxstats", str(bam)]
    r = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    if r.returncode != 0:
        return "unknown"
    for line in r.stdout.splitlines():
        if not line.strip():
            continue
        chrom = line.split("\t", 1)[0]
        if chrom == "*":
            continue
        return "chr" if chrom.startswith("chr") else "no_chr"
    return "unknown"


def open_text(path: Path):
    if str(path).endswith(".gz"):
        return gzip.open(path, "rt")
    return path.open("r")


def normalize_chr(c: str) -> str:
    x = c.strip()
    if x.lower().startswith("chr"):
        return x[3:].upper()
    return x.upper()


def parse_config(path: Path) -> dict:
    with path.open("r") as fh:
        cfg = yaml.safe_load(fh) or {}
    if not isinstance(cfg, dict):
        raise RuntimeError("config.yaml 必须是 YAML 字典")
    return cfg


def probe_resource_file(root_dir: Path, explicit_path: str, keywords: list[str], label: str) -> Path:
    if explicit_path:
        p = Path(explicit_path)
        if not p.exists():
            raise RuntimeError(f"ascat_resources.{label} 不存在: {p}")
        return p

    if not root_dir.exists():
        raise RuntimeError(f"ascat_resources.root_dir 不存在: {root_dir}")

    candidates = []
    for p in root_dir.rglob("*"):
        if not p.is_file():
            continue
        low = p.name.lower()
        if low.endswith(".zip"):
            continue
        if all(k in low for k in keywords):
            candidates.append(p)

    if len(candidates) == 1:
        return candidates[0]
    if not candidates:
        raise RuntimeError(
            f"无法在 {root_dir} 自动识别 {label}，请在 config 里显式设置 ascat_resources.{label}"
        )

    # 倾向更短文件名（通常是主资源，而不是 index/aux）
    candidates = sorted(candidates, key=lambda p: (len(p.name), str(p)))
    return candidates[0]


def read_first_chr_style(path: Path) -> str:
    with open_text(path) as fh:
        for line in fh:
            s = line.strip()
            if not s or s.startswith("#"):
                continue
            cols = s.split("\t")
            if not cols:
                continue
            c0 = cols[0]
            if normalize_chr(c0) in {str(i) for i in range(1, 23)} | {"X", "Y", "MT", "M"}:
                return "chr" if c0.startswith("chr") else "no_chr"
    return "unknown"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--pairbam", required=True)
    ap.add_argument("--config", required=True)
    ap.add_argument("--outdir", required=True)
    args = ap.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    errors = []
    infos = []
    sample_ids = set()

    cfg = parse_config(Path(args.config))
    ascat_res = cfg.get("ascat_resources") or {}
    root_dir = Path(ascat_res.get("root_dir", "")) if ascat_res.get("root_dir") else None

    loci_path = None
    alleles_path = None
    gc_path = None
    try:
        if not root_dir:
            raise RuntimeError("缺少 ascat_resources.root_dir")
        loci_path = probe_resource_file(root_dir, str(ascat_res.get("loci_path", "") or ""), ["loci", "hg38"], "loci_path")
        alleles_path = probe_resource_file(root_dir, str(ascat_res.get("alleles_path", "") or ""), ["alleles", "hg38"], "alleles_path")
        gc_path = probe_resource_file(root_dir, str(ascat_res.get("gc_path", "") or ""), ["gc", "hg38"], "gc_path")
        infos.extend(
            [
                f"ASCAT资源识别: loci_path={loci_path}",
                f"ASCAT资源识别: alleles_path={alleles_path}",
                f"ASCAT资源识别: gc_path={gc_path}",
            ]
        )
    except Exception as e:
        errors.append(str(e))

    resource_chr_style = "unknown"
    if loci_path is not None:
        resource_chr_style = read_first_chr_style(loci_path)
        infos.append(f"ASCAT资源chr风格: {resource_chr_style}")

    with open(args.pairbam) as fh:
        reader = csv.reader(fh, delimiter="\t")
        header = next(reader)
        if header[:3] != ["sample_id", "tumor_bam", "normal_bam"]:
            errors.append("前三列必须是 sample_id tumor_bam normal_bam")
        for i, row in enumerate(reader, start=2):
            if len(row) < 3:
                errors.append(f"line {i}: 列数不足")
                continue
            sid, tb, nb = row[0], Path(row[1]), Path(row[2])
            if sid in sample_ids:
                errors.append(f"line {i}: sample_id 重复 {sid}")
            sample_ids.add(sid)
            if not str(tb).startswith("/") or not str(nb).startswith("/"):
                errors.append(f"line {i}: BAM 不是绝对路径")
            if tb == nb:
                errors.append(f"line {i}: tumor_bam 与 normal_bam 相同")
            for b in (tb, nb):
                if not b.exists():
                    errors.append(f"line {i}: BAM 不存在 {b}")
                bai = Path(str(b) + ".bai")
                if not bai.exists():
                    errors.append(f"line {i}: BAI 不存在 {bai}")
                if b.exists() and not check_bam_readable(b):
                    errors.append(f"line {i}: BAM 不可读 {b}")

            if tb.exists():
                bam_style = detect_bam_chr_style(tb)
                infos.append(f"sample={sid} tumor_chr_style={bam_style}")
                if resource_chr_style != "unknown" and bam_style != "unknown" and bam_style != resource_chr_style:
                    errors.append(
                        f"line {i}: BAM染色体命名({bam_style})与ASCAT资源({resource_chr_style})不一致; "
                        f"请调整 ascat_resources.chr_style 或替换对应资源"
                    )

    report = outdir / "input_check.report.tsv"
    with report.open("w") as fw:
        fw.write("status\tmessage\n")
        for m in infos:
            fw.write(f"INFO\t{m}\n")
        if errors:
            for e in errors:
                fw.write(f"ERROR\t{e}\n")
        else:
            fw.write("OK\tinput check passed\n")

    if errors:
        sys.exit(1)


if __name__ == "__main__":
    main()
