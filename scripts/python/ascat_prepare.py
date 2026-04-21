#!/usr/bin/env python3
import argparse
import csv
import gzip
import hashlib
import math
import random
import subprocess
import sys
import tempfile
from collections import defaultdict
from pathlib import Path

import yaml


def log(msg):
    print(f"[ascat_prepare] {msg}", flush=True)


def load_config(config_path):
    with open(config_path, "r") as f:
        cfg = yaml.safe_load(f)

    if not isinstance(cfg, dict):
        raise RuntimeError("Config must be a YAML mapping")

    reference = cfg.get("reference") or {}
    if "snp_vcf" not in reference or not reference["snp_vcf"]:
        raise RuntimeError("Missing required config: reference.snp_vcf")

    prep_cfg = cfg.get("ascat_prepare") or {}
    out = {
        "snp_vcf": str(reference["snp_vcf"]),
        "max_sites": int(prep_cfg.get("max_sites", 50000)),
        "site_bed": str(prep_cfg.get("site_bed", "") or ""),
        "min_normal_depth": int(prep_cfg.get("min_normal_depth", 12)),
        "min_tumor_depth": int(prep_cfg.get("min_tumor_depth", 8)),
        "include_chroms": prep_cfg.get(
            "include_chroms",
            [*(f"chr{i}" for i in range(1, 23)), "chrX", "chrY"],
        ),
        "normal_het_min_baf": float(prep_cfg.get("normal_het_min_baf", 0.30)),
        "normal_het_max_baf": float(prep_cfg.get("normal_het_max_baf", 0.70)),
    }

    if out["max_sites"] < 1:
        raise RuntimeError("ascat_prepare.max_sites must be >= 1")
    if out["min_normal_depth"] < 1 or out["min_tumor_depth"] < 1:
        raise RuntimeError("ascat_prepare.min_normal_depth/min_tumor_depth must be >= 1")
    if not isinstance(out["include_chroms"], list) or not out["include_chroms"]:
        raise RuntimeError("ascat_prepare.include_chroms must be a non-empty list")

    return out


def open_text(path):
    if str(path).endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path, "rt")


def norm_chrom(c):
    c2 = c.strip()
    if c2.lower().startswith("chr"):
        return c2[3:].upper()
    return c2.upper()


def load_sites_from_site_bed(site_bed, include_chroms, max_sites):
    include = {norm_chrom(c) for c in include_chroms}
    sites = []
    with open(site_bed, "r") as f:
        for line in f:
            if not line.strip() or line.startswith("#"):
                continue
            p = line.rstrip("\n").split("\t")
            if len(p) < 2:
                continue
            chrom = p[0]

            pos = None
            if len(p) >= 3:
                # BED style: chrom, start0, end
                try:
                    pos = int(p[1]) + 1
                except ValueError:
                    pos = None
            if pos is None:
                try:
                    pos = int(p[1])
                except ValueError:
                    continue

            if norm_chrom(chrom) not in include:
                continue
            sites.append((chrom, pos))
            if len(sites) >= max_sites:
                break
    return sites


def sample_sites_from_vcf(vcf_path, include_chroms, max_sites, cache_path):
    include = {norm_chrom(c) for c in include_chroms}

    if cache_path.exists():
        log(f"using cached sites: {cache_path}")
        out = []
        with cache_path.open() as f:
            for row in csv.DictReader(f, delimiter="\t"):
                out.append((row["chr"], int(row["pos"])))
        if out:
            return out

    # per-chromosome reservoirs to avoid taking only file head
    per_chr_target = max(1, max_sites // max(1, len(include)))
    reservoirs = defaultdict(list)
    seen = defaultdict(int)

    log(f"sampling sites from VCF: {vcf_path}")
    with open_text(vcf_path) as f:
        for line in f:
            if not line or line.startswith("#"):
                continue
            p = line.rstrip("\n").split("\t")
            if len(p) < 2:
                continue
            chrom = p[0]
            if norm_chrom(chrom) not in include:
                continue
            try:
                pos = int(p[1])
            except ValueError:
                continue

            k = norm_chrom(chrom)
            seen[k] += 1
            bucket = reservoirs[k]
            if len(bucket) < per_chr_target:
                bucket.append((chrom, pos))
            else:
                j = random.randint(1, seen[k])
                if j <= per_chr_target:
                    bucket[j - 1] = (chrom, pos)

    combined = []
    for k in sorted(reservoirs):
        combined.extend(reservoirs[k])

    if len(combined) > max_sites:
        random.shuffle(combined)
        combined = combined[:max_sites]

    cache_path.parent.mkdir(parents=True, exist_ok=True)
    with cache_path.open("w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["chr", "pos"])
        for chrom, pos in combined:
            w.writerow([chrom, pos])

    return combined


def parse_pileup_bases(ref_base, bases):
    ref = ref_base.upper()
    counts = {"A": 0, "C": 0, "G": 0, "T": 0}
    i = 0
    n = len(bases)
    while i < n:
        c = bases[i]
        if c == "^":
            i += 2
            continue
        if c == "$":
            i += 1
            continue
        if c in "+-":
            i += 1
            num = []
            while i < n and bases[i].isdigit():
                num.append(bases[i])
                i += 1
            indel_len = int("".join(num)) if num else 0
            i += indel_len
            continue
        if c in ".,":
            if ref in counts:
                counts[ref] += 1
            i += 1
            continue

        b = c.upper()
        if b in counts:
            counts[b] += 1
        i += 1

    return counts


def select_alt_base(normal_counts, ref_base):
    ref = ref_base.upper()
    candidates = [(b, c) for b, c in normal_counts.items() if b != ref]
    if not candidates:
        return None, 0
    alt, c = max(candidates, key=lambda x: x[1])
    return alt, c


def run_mpileup_and_collect(tumor_bam, normal_bam, sites, cfg):
    rows = []
    with tempfile.NamedTemporaryFile("w", suffix=".sites", delete=True) as sf:
        for chrom, pos in sites:
            sf.write(f"{chrom}\t{pos}\n")
        sf.flush()

        cmd = [
            "samtools",
            "mpileup",
            "-aa",
            "-A",
            "-B",
            "-Q",
            "20",
            "-q",
            "20",
            "-l",
            sf.name,
            normal_bam,
            tumor_bam,
        ]

        log("running samtools mpileup (normal + tumor)")
        proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        assert proc.stdout is not None

        n_in = 0
        for line in proc.stdout:
            p = line.rstrip("\n").split("\t")
            # chr pos ref n_dp n_bases n_qual t_dp t_bases t_qual
            if len(p) < 9:
                continue
            n_in += 1
            chrom, pos_s, ref = p[0], p[1], p[2]
            try:
                pos = int(pos_s)
                n_dp = int(p[3])
                t_dp = int(p[6])
            except ValueError:
                continue

            if n_dp < cfg["min_normal_depth"] or t_dp < cfg["min_tumor_depth"]:
                continue

            n_counts = parse_pileup_bases(ref, p[4])
            t_counts = parse_pileup_bases(ref, p[7])

            alt_base, alt_count_n = select_alt_base(n_counts, ref)
            if not alt_base:
                continue

            ref_u = ref.upper()
            ref_count_n = n_counts.get(ref_u, 0)
            n_ref_alt = ref_count_n + alt_count_n
            if n_ref_alt == 0:
                continue

            n_alt_frac = alt_count_n / n_ref_alt
            if n_alt_frac < cfg["normal_het_min_baf"] or n_alt_frac > cfg["normal_het_max_baf"]:
                continue

            # tumor BAF proxy: alt fraction at normal-selected het SNP
            t_alt = t_counts.get(alt_base, 0)
            t_ref = t_counts.get(ref_u, 0)
            t_ref_alt = t_ref + t_alt
            if t_ref_alt == 0:
                continue

            baf = t_alt / t_ref_alt
            logr = math.log2((t_dp + 1.0) / (n_dp + 1.0))
            rows.append((chrom, pos, baf, logr, n_dp, t_dp, ref_u, alt_base))

            if n_in % 5000 == 0:
                log(f"mpileup parsed sites: {n_in}, retained: {len(rows)}")

        stderr = proc.stderr.read() if proc.stderr else ""
        rc = proc.wait()
        if rc != 0:
            raise RuntimeError(f"samtools mpileup failed (rc={rc}): {stderr.strip()}")

    return rows


def write_outputs(outdir, rows):
    log("writing baf.tsv")
    with (outdir / "baf.tsv").open("w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["chr", "pos", "baf"])
        for chrom, pos, baf, *_ in rows:
            w.writerow([chrom, pos, f"{baf:.6f}"])

    log("writing logr.tsv")
    with (outdir / "logr.tsv").open("w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["chr", "pos", "logr"])
        for chrom, pos, _, logr, *_ in rows:
            w.writerow([chrom, pos, f"{logr:.6f}"])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--sample-id", required=True)
    ap.add_argument("--tumor-bam", required=True)
    ap.add_argument("--normal-bam", required=True)
    ap.add_argument("--config", required=True)
    ap.add_argument("--outdir", required=True)
    args = ap.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    log("loading config")
    cfg = load_config(args.config)

    snp_vcf = cfg["snp_vcf"]
    max_sites = cfg["max_sites"]
    include_chroms = cfg["include_chroms"]

    site_source = ""
    requested_sites = []
    if cfg["site_bed"]:
        site_bed = Path(cfg["site_bed"])
        if not site_bed.exists():
            raise RuntimeError(f"ascat_prepare.site_bed not found: {site_bed}")
        log(f"loading SNP sites from site_bed: {site_bed}")
        requested_sites = load_sites_from_site_bed(site_bed, include_chroms, max_sites)
        site_source = f"site_bed:{site_bed}"
    else:
        key = f"{snp_vcf}|{max_sites}|{','.join(include_chroms)}"
        digest = hashlib.md5(key.encode()).hexdigest()[:10]
        cache_path = Path(args.config).resolve().parent / f"ascat_prepare.sites.{digest}.tsv"
        requested_sites = sample_sites_from_vcf(snp_vcf, include_chroms, max_sites, cache_path)
        site_source = f"vcf_uniform_reservoir:{snp_vcf};cache={cache_path.name}"

    if not requested_sites:
        raise RuntimeError("No SNP sites prepared for ascat_prepare (check site_bed/vcf and include_chroms)")

    log(f"prepared SNP sites: {len(requested_sites)}")
    rows = run_mpileup_and_collect(args.tumor_bam, args.normal_bam, requested_sites, cfg)
    log(f"retained sites after filters: {len(rows)}")

    write_outputs(outdir, rows)

    log("writing ascat_prepare.run_info.tsv")
    with (outdir / "ascat_prepare.run_info.tsv").open("w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow([
            "sample_id",
            "status",
            "site_source",
            "n_sites_requested",
            "n_sites_retained",
            "baf_method",
            "logr_method",
            "min_normal_depth",
            "min_tumor_depth",
            "normal_het_min_baf",
            "normal_het_max_baf",
            "notes",
        ])
        w.writerow([
            args.sample_id,
            "ok" if rows else "empty",
            site_source,
            len(requested_sites),
            len(rows),
            "proxy_baf_from_tumor_alt_fraction_at_normal_het_sites(mpileup_counts)",
            "proxy_logr_log2((tumor_depth+1)/(normal_depth+1))",
            cfg["min_normal_depth"],
            cfg["min_tumor_depth"],
            cfg["normal_het_min_baf"],
            cfg["normal_het_max_baf"],
            "Proxy implementation for robust production-oriented preprocessing; not full ASCAT allele counting pipeline",
        ])


if __name__ == "__main__":
    try:
        main()
    except Exception as e:
        print(f"[ascat_prepare][ERROR] {e}", file=sys.stderr, flush=True)
        raise
