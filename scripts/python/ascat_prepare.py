#!/usr/bin/env python3
import argparse
import csv
import gzip
import math
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

import yaml


VALID_AUTOSOMES = {str(i) for i in range(1, 23)} | {"X", "Y", "M", "MT"}


def log(msg):
    print(f"[ascat_prepare] {msg}", flush=True)


def open_text(path: Path):
    if str(path).endswith(".gz"):
        return gzip.open(path, "rt")
    return path.open("r")


def norm_chr(c: str) -> str:
    x = c.strip()
    if x.lower().startswith("chr"):
        return x[3:].upper()
    return x.upper()


def infer_chr_style(chrom: str) -> str:
    return "chr" if chrom.startswith("chr") else "no_chr"


def convert_chr_style(chrom: str, style: str) -> str:
    core = norm_chr(chrom)
    if style == "chr":
        return f"chr{core}"
    return core


def load_yaml(path: Path) -> dict:
    with path.open("r") as fh:
        cfg = yaml.safe_load(fh) or {}
    if not isinstance(cfg, dict):
        raise RuntimeError("Config must be a YAML mapping")
    return cfg


def probe_resource_file(root_dir: Path, explicit_path: str, keywords: list[str], label: str) -> Path:
    if explicit_path:
        p = Path(explicit_path)
        if not p.exists():
            raise RuntimeError(f"ascat_resources.{label} not found: {p}")
        return p

    if not root_dir.exists():
        raise RuntimeError(f"ascat_resources.root_dir not found: {root_dir}")

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
            f"unable to auto-detect ascat_resources.{label} from root_dir={root_dir}; set explicit path in config"
        )

    candidates = sorted(candidates, key=lambda p: (len(p.name), str(p)))
    log(f"multiple {label} candidates found, choosing shortest name: {candidates[0]}")
    return candidates[0]


def parse_config(config_path: Path) -> dict:
    cfg = load_yaml(config_path)

    ascat_resources = cfg.get("ascat_resources") or {}
    if not ascat_resources:
        raise RuntimeError("Missing required config section: ascat_resources")

    root_dir_v = ascat_resources.get("root_dir")
    if not root_dir_v:
        raise RuntimeError("Missing required config: ascat_resources.root_dir")
    root_dir = Path(str(root_dir_v))

    loci_path = probe_resource_file(root_dir, str(ascat_resources.get("loci_path", "") or ""), ["loci", "hg38"], "loci_path")
    alleles_path = probe_resource_file(
        root_dir, str(ascat_resources.get("alleles_path", "") or ""), ["alleles", "hg38"], "alleles_path"
    )
    gc_path = probe_resource_file(root_dir, str(ascat_resources.get("gc_path", "") or ""), ["gc", "hg38"], "gc_path")

    prep_cfg = cfg.get("ascat_prepare") or {}

    return {
        "root_dir": root_dir,
        "loci_path": loci_path,
        "alleles_path": alleles_path,
        "gc_path": gc_path,
        "chr_style": str(ascat_resources.get("chr_style", "chr")),
        "allelecounter_exe": str(ascat_resources.get("allelecounter_exe", "alleleCounter") or "alleleCounter"),
        "min_normal_depth": int(prep_cfg.get("min_normal_depth", 12)),
        "min_tumor_depth": int(prep_cfg.get("min_tumor_depth", 8)),
        "normal_het_min_baf": float(prep_cfg.get("normal_het_min_baf", 0.30)),
        "normal_het_max_baf": float(prep_cfg.get("normal_het_max_baf", 0.70)),
        "max_sites": int(prep_cfg.get("max_sites", 50000)),
        "include_chroms": prep_cfg.get(
            "include_chroms", [*(f"chr{i}" for i in range(1, 23)), "chrX", "chrY"]
        ),
    }


def detect_bam_chr_style(bam_path: str) -> str:
    cmd = ["samtools", "idxstats", bam_path]
    r = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    if r.returncode != 0:
        raise RuntimeError(f"samtools idxstats failed for {bam_path}: {r.stderr.strip()}")
    for line in r.stdout.splitlines():
        if not line.strip():
            continue
        chrom = line.split("\t", 1)[0]
        if chrom == "*":
            continue
        return infer_chr_style(chrom)
    raise RuntimeError(f"Unable to determine chromosome style from bam: {bam_path}")


def sniff_header(path: Path) -> tuple[bool, list[str]]:
    with open_text(path) as fh:
        for line in fh:
            s = line.strip()
            if not s or s.startswith("#"):
                continue
            cols = s.split("\t")
            first = cols[0].lower()
            if first in {"chr", "chrom", "chromosome", "contig"}:
                return True, cols
            if first.startswith("rs"):
                return False, []
            return False, []
    return False, []


def load_loci(loci_path: Path, target_style: str, include_chroms: list[str], max_sites: int):
    include = {norm_chr(c) for c in include_chroms}
    has_header, _ = sniff_header(loci_path)
    rows = []
    with open_text(loci_path) as fh:
        for line in fh:
            s = line.strip()
            if not s or s.startswith("#"):
                continue
            cols = s.split("\t")
            if has_header:
                has_header = False
                continue
            if len(cols) < 2:
                continue
            chrom = cols[0]
            try:
                pos = int(cols[1])
            except ValueError:
                continue
            ccore = norm_chr(chrom)
            if ccore not in include or ccore not in VALID_AUTOSOMES:
                continue
            rows.append((convert_chr_style(chrom, target_style), pos))
            if len(rows) >= max_sites:
                break
    return rows


def find_col(header_map: dict, names: list[str]):
    for n in names:
        if n in header_map:
            return header_map[n]
    return None


def load_alleles(alleles_path: Path):
    has_header, header = sniff_header(alleles_path)
    allele_map = {}

    with open_text(alleles_path) as fh:
        if has_header:
            h = header
            header_map = {c.strip().lower(): i for i, c in enumerate(h)}
            chr_i = find_col(header_map, ["chromosome", "chrom", "chr", "contig"])
            pos_i = find_col(header_map, ["position", "pos"])
            ref_i = find_col(header_map, ["ref", "ref_allele", "a", "allele_a"])
            alt_i = find_col(header_map, ["alt", "alt_allele", "b", "allele_b"])
            if None in {chr_i, pos_i, ref_i, alt_i}:
                raise RuntimeError("Cannot parse alleles header; expected chr/pos/ref/alt columns")

            next(fh)
            for line in fh:
                s = line.strip()
                if not s or s.startswith("#"):
                    continue
                cols = s.split("\t")
                if len(cols) <= max(chr_i, pos_i, ref_i, alt_i):
                    continue
                chrom = cols[chr_i]
                try:
                    pos = int(cols[pos_i])
                except ValueError:
                    continue
                ref = cols[ref_i].upper()
                alt = cols[alt_i].upper()
                key = (norm_chr(chrom), pos)
                allele_map[key] = (ref, alt)
        else:
            for line in fh:
                s = line.strip()
                if not s or s.startswith("#"):
                    continue
                cols = s.split("\t")
                if len(cols) < 4:
                    continue
                chrom = cols[0]
                try:
                    pos = int(cols[1])
                except ValueError:
                    continue
                ref = cols[2].upper()
                alt = cols[3].upper()
                key = (norm_chr(chrom), pos)
                allele_map[key] = (ref, alt)

    if not allele_map:
        raise RuntimeError(f"No allele records parsed from: {alleles_path}")
    return allele_map


def load_gc(gc_path: Path):
    has_header, header = sniff_header(gc_path)
    gc_map = {}

    with open_text(gc_path) as fh:
        if has_header:
            hmap = {c.strip().lower(): i for i, c in enumerate(header)}
            chr_i = find_col(hmap, ["chromosome", "chrom", "chr", "contig"])
            pos_i = find_col(hmap, ["position", "pos"])
            gc_i = find_col(hmap, ["gc", "gc_content", "gcpct"])
            if None in {chr_i, pos_i, gc_i}:
                raise RuntimeError("Cannot parse gc header; expected chr/pos/gc columns")
            next(fh)
            for line in fh:
                cols = line.strip().split("\t")
                if len(cols) <= max(chr_i, pos_i, gc_i):
                    continue
                chrom = cols[chr_i]
                try:
                    pos = int(cols[pos_i])
                    gc = float(cols[gc_i])
                except ValueError:
                    continue
                gc_map[(norm_chr(chrom), pos)] = gc
        else:
            for line in fh:
                cols = line.strip().split("\t")
                if len(cols) < 3:
                    continue
                try:
                    pos = int(cols[1])
                    gc = float(cols[2])
                except ValueError:
                    continue
                gc_map[(norm_chr(cols[0]), pos)] = gc

    if not gc_map:
        raise RuntimeError(f"No gc records parsed from: {gc_path}")
    return gc_map


def parse_pileup_bases(ref_base: str, bases: str):
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


def run_mpileup_counts(tumor_bam: str, normal_bam: str, loci: list[tuple[str, int]], allele_map: dict):
    counts = {}
    with tempfile.NamedTemporaryFile("w", suffix=".loci", delete=True) as lf:
        for chrom, pos in loci:
            lf.write(f"{chrom}\t{pos}\n")
        lf.flush()

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
            lf.name,
            normal_bam,
            tumor_bam,
        ]
        proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        assert proc.stdout is not None
        for line in proc.stdout:
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 9:
                continue
            chrom = cols[0]
            try:
                pos = int(cols[1])
                n_dp = int(cols[3])
                t_dp = int(cols[6])
            except ValueError:
                continue
            key = (norm_chr(chrom), pos)
            ref_alt = allele_map.get(key)
            if not ref_alt:
                continue
            ref, alt = ref_alt
            n_counts = parse_pileup_bases(ref, cols[4])
            t_counts = parse_pileup_bases(ref, cols[7])
            counts[key] = {
                "normal_ref": n_counts.get(ref, 0),
                "normal_alt": n_counts.get(alt, 0),
                "tumor_ref": t_counts.get(ref, 0),
                "tumor_alt": t_counts.get(alt, 0),
                "normal_dp": n_dp,
                "tumor_dp": t_dp,
            }
        stderr = proc.stderr.read() if proc.stderr else ""
        rc = proc.wait()
        if rc != 0:
            raise RuntimeError(f"samtools mpileup failed (rc={rc}): {stderr.strip()}")
    return counts


def has_allelecounter(exe: str) -> bool:
    return shutil.which(exe) is not None


def run_allelecounter_single(exe: str, bam: str, loci: list[tuple[str, int]], out_path: Path):
    with tempfile.NamedTemporaryFile("w", suffix=".loci", delete=True) as lf:
        for chrom, pos in loci:
            lf.write(f"{chrom}\t{pos}\n")
        lf.flush()
        cmd = [exe, "-b", bam, "-l", lf.name, "-o", str(out_path), "-q", "20", "-m", "20"]
        r = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        if r.returncode != 0:
            raise RuntimeError(f"alleleCounter failed: {' '.join(cmd)}\n{r.stderr.strip()}")


def parse_allelecounter_output(path: Path):
    rows = {}
    with path.open("r") as fh:
        first = fh.readline().strip().split("\t")
        header_mode = len(first) > 0 and first[0].lower() in {"chromosome", "chrom", "chr"}
        if header_mode:
            hmap = {c.strip().lower(): i for i, c in enumerate(first)}
            chr_i = find_col(hmap, ["chromosome", "chrom", "chr"])
            pos_i = find_col(hmap, ["position", "pos"])
            ref_i = find_col(hmap, ["refcount", "ref_count", "acount", "a"])
            alt_i = find_col(hmap, ["altcount", "alt_count", "bcount", "b"])
            if None in {chr_i, pos_i, ref_i, alt_i}:
                raise RuntimeError(f"Unknown alleleCounter header in {path}")
        else:
            chr_i, pos_i, ref_i, alt_i = 0, 1, 2, 3
            cols = first
            if len(cols) >= 4:
                chrom = cols[chr_i]
                pos = int(cols[pos_i])
                ref_c = int(float(cols[ref_i]))
                alt_c = int(float(cols[alt_i]))
                rows[(norm_chr(chrom), pos)] = (ref_c, alt_c)

        for line in fh:
            cols = line.strip().split("\t")
            if len(cols) <= max(chr_i, pos_i, ref_i, alt_i):
                continue
            try:
                chrom = cols[chr_i]
                pos = int(cols[pos_i])
                ref_c = int(float(cols[ref_i]))
                alt_c = int(float(cols[alt_i]))
            except ValueError:
                continue
            rows[(norm_chr(chrom), pos)] = (ref_c, alt_c)
    return rows


def run_allelecounter_counts(exe: str, tumor_bam: str, normal_bam: str, loci: list[tuple[str, int]]):
    with tempfile.TemporaryDirectory() as td:
        n_out = Path(td) / "normal.tsv"
        t_out = Path(td) / "tumor.tsv"
        run_allelecounter_single(exe, normal_bam, loci, n_out)
        run_allelecounter_single(exe, tumor_bam, loci, t_out)
        n_counts = parse_allelecounter_output(n_out)
        t_counts = parse_allelecounter_output(t_out)

    out = {}
    keys = set(n_counts) | set(t_counts)
    for k in keys:
        nr, na = n_counts.get(k, (0, 0))
        tr, ta = t_counts.get(k, (0, 0))
        out[k] = {
            "normal_ref": nr,
            "normal_alt": na,
            "tumor_ref": tr,
            "tumor_alt": ta,
            "normal_dp": nr + na,
            "tumor_dp": tr + ta,
        }
    return out


def linear_gc_correct(rows, gc_map):
    xs = []
    ys = []
    for r in rows:
        gc = gc_map.get((norm_chr(r["chr"]), r["pos"]))
        if gc is None:
            continue
        xs.append(gc)
        ys.append(r["raw_logr"])

    if len(xs) < 100:
        return rows, "fallback_raw_logr(no_enough_gc_overlap)", True

    n = len(xs)
    sx = sum(xs)
    sy = sum(ys)
    sxx = sum(x * x for x in xs)
    sxy = sum(x * y for x, y in zip(xs, ys))
    den = n * sxx - sx * sx
    if den == 0:
        return rows, "fallback_raw_logr(gc_regression_singular)", True

    slope = (n * sxy - sx * sy) / den
    intercept = (sy - slope * sx) / n

    med = sorted(ys)[len(ys) // 2]
    corrected = []
    for r in rows:
        gc = gc_map.get((norm_chr(r["chr"]), r["pos"]))
        if gc is None:
            corr = r["raw_logr"]
        else:
            pred = intercept + slope * gc
            corr = (r["raw_logr"] - pred) + med
        x = dict(r)
        x["logr"] = corr
        corrected.append(x)
    return corrected, "log2_tn_depth_ratio+linear_gc_correction", False


def build_rows(loci, allele_map, counts, cfg):
    rows = []
    for chrom, pos in loci:
        key = (norm_chr(chrom), pos)
        alleles = allele_map.get(key)
        c = counts.get(key)
        if not alleles or not c:
            continue

        n_dp = c["normal_dp"]
        t_dp = c["tumor_dp"]
        if n_dp < cfg["min_normal_depth"] or t_dp < cfg["min_tumor_depth"]:
            continue

        n_ref_alt = c["normal_ref"] + c["normal_alt"]
        t_ref_alt = c["tumor_ref"] + c["tumor_alt"]
        if n_ref_alt == 0 or t_ref_alt == 0:
            continue

        n_baf = c["normal_alt"] / n_ref_alt
        if n_baf < cfg["normal_het_min_baf"] or n_baf > cfg["normal_het_max_baf"]:
            continue

        t_baf = c["tumor_alt"] / t_ref_alt
        raw_logr = math.log2((t_dp + 1.0) / (n_dp + 1.0))

        rows.append(
            {
                "chr": chrom,
                "pos": pos,
                "baf": t_baf,
                "raw_logr": raw_logr,
                "normal_dp": n_dp,
                "tumor_dp": t_dp,
            }
        )
    return rows


def write_outputs(outdir: Path, rows):
    with (outdir / "baf.tsv").open("w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["chr", "pos", "baf"])
        for r in rows:
            w.writerow([r["chr"], r["pos"], f"{r['baf']:.6f}"])

    with (outdir / "logr.tsv").open("w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["chr", "pos", "logr"])
        for r in rows:
            w.writerow([r["chr"], r["pos"], f"{r['logr']:.6f}"])


def write_run_info(outdir: Path, sample_id: str, cfg: dict, loci_n: int, rows_n: int, baf_method: str, logr_method: str, fallback: bool):
    with (outdir / "ascat_prepare.run_info.tsv").open("w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(
            [
                "sample_id",
                "status",
                "root_dir",
                "loci_path",
                "alleles_path",
                "gc_path",
                "n_sites_requested",
                "n_sites_retained",
                "baf_method",
                "logr_method",
                "fallback",
                "min_normal_depth",
                "min_tumor_depth",
                "normal_het_min_baf",
                "normal_het_max_baf",
                "notes",
            ]
        )
        w.writerow(
            [
                sample_id,
                "ok" if rows_n > 0 else "empty",
                str(cfg["root_dir"]),
                str(cfg["loci_path"]),
                str(cfg["alleles_path"]),
                str(cfg["gc_path"]),
                loci_n,
                rows_n,
                baf_method,
                logr_method,
                "yes" if fallback else "no",
                cfg["min_normal_depth"],
                cfg["min_tumor_depth"],
                cfg["normal_het_min_baf"],
                cfg["normal_het_max_baf"],
                "Sarek/ASCAT-aligned scaffold: fixed loci/alleles/GC resources; normal-guided het selection; tumor BAF + depth-based logR",
            ]
        )


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

    cfg = parse_config(Path(args.config))

    bam_style = detect_bam_chr_style(args.tumor_bam)
    if cfg["chr_style"] not in {"chr", "no_chr"}:
        raise RuntimeError("ascat_resources.chr_style must be chr or no_chr")
    if bam_style != cfg["chr_style"]:
        raise RuntimeError(
            f"BAM chromosome style ({bam_style}) inconsistent with ascat_resources.chr_style ({cfg['chr_style']})"
        )

    loci = load_loci(cfg["loci_path"], cfg["chr_style"], cfg["include_chroms"], cfg["max_sites"])
    if not loci:
        raise RuntimeError("No loci loaded from ascat_resources.loci_path")
    log(f"loaded loci: {len(loci)}")

    allele_map = load_alleles(cfg["alleles_path"])
    gc_map = load_gc(cfg["gc_path"])
    log(f"loaded alleles: {len(allele_map)}")
    log(f"loaded gc entries: {len(gc_map)}")

    fallback = False
    baf_method = "alleleCounter_fixed_loci_with_normal_het_filter"
    if has_allelecounter(cfg["allelecounter_exe"]):
        try:
            log(f"using alleleCounter: {cfg['allelecounter_exe']}")
            counts = run_allelecounter_counts(cfg["allelecounter_exe"], args.tumor_bam, args.normal_bam, loci)
        except Exception as e:
            log(f"alleleCounter failed; fallback to mpileup. reason: {e}")
            counts = run_mpileup_counts(args.tumor_bam, args.normal_bam, loci, allele_map)
            baf_method = "mpileup_fixed_loci_with_normal_het_filter(fallback)"
            fallback = True
    else:
        log("alleleCounter not found; fallback to mpileup")
        counts = run_mpileup_counts(args.tumor_bam, args.normal_bam, loci, allele_map)
        baf_method = "mpileup_fixed_loci_with_normal_het_filter(fallback)"
        fallback = True

    rows = build_rows(loci, allele_map, counts, cfg)
    log(f"retained loci after filters: {len(rows)}")
    rows, logr_method, gc_fallback = linear_gc_correct(rows, gc_map)
    fallback = fallback or gc_fallback

    write_outputs(outdir, rows)
    write_run_info(outdir, args.sample_id, cfg, len(loci), len(rows), baf_method, logr_method, fallback)


if __name__ == "__main__":
    try:
        main()
    except Exception as e:
        print(f"[ascat_prepare][ERROR] {e}", file=sys.stderr, flush=True)
        raise
