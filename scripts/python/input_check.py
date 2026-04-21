#!/usr/bin/env python3
import argparse
import csv
from pathlib import Path
import subprocess
import sys


def check_bam_readable(bam: Path) -> bool:
    cmd = ["samtools", "quickcheck", str(bam)]
    r = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    return r.returncode == 0


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--pairbam", required=True)
    ap.add_argument("--config", required=True)
    ap.add_argument("--outdir", required=True)
    args = ap.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    errors = []
    sample_ids = set()

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

    report = outdir / "input_check.report.tsv"
    with report.open("w") as fw:
        fw.write("status\tmessage\n")
        if errors:
            for e in errors:
                fw.write(f"ERROR\t{e}\n")
        else:
            fw.write("OK\tinput check passed\n")

    if errors:
        sys.exit(1)


if __name__ == "__main__":
    main()
