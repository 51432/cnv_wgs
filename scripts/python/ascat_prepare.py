#!/usr/bin/env python3
import argparse
import csv
import gzip
import math
import subprocess
import tempfile
from pathlib import Path


def get_config_value(config_path, key):
    prefix = f"{key}:"
    for line in Path(config_path).read_text().splitlines():
        stripped = line.strip()
        if stripped.startswith(prefix):
            return stripped.split(":", 1)[1].strip().strip('"').strip("'")
    return ""


def open_text_maybe_gz(path):
    if path.endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path, "rt")


def collect_snp_sites(vcf_path, max_sites=30000):
    sites = []
    with open_text_maybe_gz(vcf_path) as f:
        for line in f:
            if not line or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 2:
                continue
            chrom, pos = parts[0], parts[1]
            if chrom in {"chrM", "MT", "M"}:
                continue
            try:
                ipos = int(pos)
            except ValueError:
                continue
            sites.append((chrom, ipos))
            if len(sites) >= max_sites:
                break
    return sites


def run_depth_for_sites(tumor_bam, normal_bam, sites):
    rows = []
    with tempfile.NamedTemporaryFile("w", suffix=".bed", delete=True) as bedf:
        for chrom, pos in sites:
            bedf.write(f"{chrom}\t{pos - 1}\t{pos}\n")
        bedf.flush()

        cmd = ["samtools", "depth", "-a", "-b", bedf.name, tumor_bam, normal_bam]
        out = subprocess.check_output(cmd, text=True)

    for line in out.splitlines():
        p = line.split("\t")
        if len(p) < 4:
            continue
        chrom, pos_s, td_s, nd_s = p[0], p[1], p[2], p[3]
        try:
            pos = int(pos_s)
            td = int(td_s)
            nd = int(nd_s)
        except ValueError:
            continue

        if nd < 8:
            continue
        baf = td / (td + nd) if (td + nd) > 0 else math.nan
        logr = math.log2((td + 1) / (nd + 1))
        rows.append((chrom, pos, baf, logr, td, nd))

    return rows


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--sample-id', required=True)
    ap.add_argument('--tumor-bam', required=True)
    ap.add_argument('--normal-bam', required=True)
    ap.add_argument('--config', required=True)
    ap.add_argument('--outdir', required=True)
    args = ap.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    snp_vcf = get_config_value(args.config, "snp_vcf")
    if not snp_vcf:
        raise RuntimeError("Cannot find reference.snp_vcf in config")

    sites = collect_snp_sites(snp_vcf)
    rows = run_depth_for_sites(args.tumor_bam, args.normal_bam, sites)

    with (outdir / 'baf.tsv').open("w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["chr", "pos", "baf"])
        for chrom, pos, baf, _, _, _ in rows:
            if math.isnan(baf):
                continue
            w.writerow([chrom, pos, f"{baf:.6f}"])

    with (outdir / 'logr.tsv').open("w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["chr", "pos", "logr"])
        for chrom, pos, _, logr, _, _ in rows:
            w.writerow([chrom, pos, f"{logr:.6f}"])

    with (outdir / 'ascat_prepare.run_info.tsv').open("w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["sample_id", "status", "n_sites_requested", "n_sites_retained", "method"])
        w.writerow([
            args.sample_id,
            "ok" if rows else "empty",
            len(sites),
            len(rows),
            "samtools_depth_proxy_baf_logr",
        ])


if __name__ == '__main__':
    main()
