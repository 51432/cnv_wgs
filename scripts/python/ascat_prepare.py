#!/usr/bin/env python3
import argparse
import csv
import gzip
import math
import shutil
import subprocess
import sys
import tempfile
from collections import defaultdict
from pathlib import Path

import yaml


DEFAULT_CORE_CHROMS = [str(i) for i in range(1, 23)] + ["X"]
EXTRA_CHROMS = {"Y", "M", "MT"}
VALID_CHROMS = set(DEFAULT_CORE_CHROMS) | EXTRA_CHROMS


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


def normalize_include_chroms(chroms: list[str], effective_style: str) -> list[str]:
    out = []
    for c in chroms:
        core = norm_chr(str(c))
        if core in VALID_CHROMS:
            out.append(convert_chr_style(core, effective_style))
    # preserve order, dedup
    seen = set()
    uniq = []
    for c in out:
        if c not in seen:
            seen.add(c)
            uniq.append(c)
    return uniq


def parse_config(config_path: Path) -> dict:
    cfg = load_yaml(config_path)

    ascat_resources = cfg.get("ascat_resources") or {}
    if not ascat_resources:
        raise RuntimeError("Missing required config section: ascat_resources")

    root_dir_v = ascat_resources.get("root_dir")
    if not root_dir_v:
        raise RuntimeError("Missing required config: ascat_resources.root_dir")
    root_dir = Path(str(root_dir_v))

    loci_dir_v = ascat_resources.get("loci_dir") or ascat_resources.get("loci_path")
    alleles_dir_v = ascat_resources.get("alleles_dir") or ascat_resources.get("alleles_path")
    gc_path_v = ascat_resources.get("gc_path")

    if not loci_dir_v or not alleles_dir_v or not gc_path_v:
        raise RuntimeError(
            "Missing resource path(s). Please provide ascat_resources.loci_dir, ascat_resources.alleles_dir and ascat_resources.gc_path"
        )

    loci_dir = Path(str(loci_dir_v))
    alleles_dir = Path(str(alleles_dir_v))
    gc_path = Path(str(gc_path_v))

    if not loci_dir.exists() or not loci_dir.is_dir():
        if loci_dir.exists() and loci_dir.is_file():
            raise RuntimeError(
                f"ascat_resources.loci_dir points to a file ({loci_dir}). Directory mode is required unless explicit one-file mode is enabled."
            )
        raise RuntimeError(f"loci_dir not found or not a directory: {loci_dir}")
    if not alleles_dir.exists() or not alleles_dir.is_dir():
        if alleles_dir.exists() and alleles_dir.is_file():
            raise RuntimeError(
                f"ascat_resources.alleles_dir points to a file ({alleles_dir}). Directory mode is required unless explicit one-file mode is enabled."
            )
        raise RuntimeError(f"alleles_dir not found or not a directory: {alleles_dir}")
    if not gc_path.exists() or not gc_path.is_file():
        raise RuntimeError(f"gc_path not found or not a file: {gc_path}")

    prep_cfg = cfg.get("ascat_prepare") or {}
    max_sites = prep_cfg.get("max_sites")
    max_sites = int(max_sites) if max_sites is not None else None

    include_raw = prep_cfg.get("include_chroms")
    if include_raw:
        include_cores = [norm_chr(str(c)) for c in include_raw]
    else:
        include_cores = list(DEFAULT_CORE_CHROMS)

    return {
        "root_dir": root_dir,
        "loci_dir": loci_dir,
        "alleles_dir": alleles_dir,
        "gc_path": gc_path,
        "gc_window": str(ascat_resources.get("gc_window", "1kb")),
        "allelecounter_exe": str(ascat_resources.get("allelecounter_exe", "alleleCounter") or "alleleCounter"),
        "min_normal_depth": int(prep_cfg.get("min_normal_depth", 12)),
        "min_tumor_depth": int(prep_cfg.get("min_tumor_depth", 8)),
        "normal_het_min_baf": float(prep_cfg.get("normal_het_min_baf", 0.30)),
        "normal_het_max_baf": float(prep_cfg.get("normal_het_max_baf", 0.70)),
        "max_sites": max_sites,
        "include_chrom_cores": include_cores,
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


def find_col(header_map: dict, names: list[str]):
    for n in names:
        if n in header_map:
            return header_map[n]
    return None


def is_chrom_file(name: str, chrom_core: str) -> bool:
    low = name.lower()
    c = chrom_core.lower()
    return (f"chr{c}" in low) or (f"_{c}." in low) or low.endswith(f"_{c}")


def list_resource_chrom_files(resource_dir: Path, kind: str, include_chrom_cores: list[str]) -> dict[str, Path]:
    files = [p for p in resource_dir.iterdir() if p.is_file() and not p.name.startswith(".")]
    out = {}
    unmatched = []
    for p in sorted(files):
        matched = None
        for core in include_chrom_cores:
            if is_chrom_file(p.name, core):
                matched = core
                break
        if matched:
            out[matched] = p
        else:
            unmatched.append(p.name)
    if unmatched:
        log(f"{kind}: files with unmatched chromosome token: {', '.join(unmatched[:20])}")
    if len(out) < 10:
        raise RuntimeError(
            f"{kind} directory appears incomplete: matched chromosome files={len(out)} (<10). dir={resource_dir}"
        )
    return out


def sniff_file_chr_style(path: Path) -> str | None:
    with open_text(path) as fh:
        for line in fh:
            s = line.strip()
            if not s or s.startswith("#"):
                continue
            cols = s.split("\t")
            if len(cols) < 2:
                continue
            if cols[0].lower() in {"chrom", "chromosome", "chr", "contig", "id"}:
                continue
            return infer_chr_style(cols[0])
    return None


def detect_resource_chr_style(loci_files: dict[str, Path], alleles_files: dict[str, Path], gc_path: Path) -> str:
    styles = []
    for d in (loci_files, alleles_files):
        for p in d.values():
            st = sniff_file_chr_style(p)
            if st:
                styles.append(st)
                break
    # GC_G1000_hg38 第一列通常是行ID，不能直接用于 chr 风格判定；此处仅在前两者都缺失时尝试。
    if not styles:
        with open_text(gc_path) as fh:
            header = fh.readline().strip().split("\t")
            hmap = {c.strip().lower(): i for i, c in enumerate(header)}
            chr_i = find_col(hmap, ["chromosome", "chrom", "chr", "contig", "chrname"])
            if chr_i is not None:
                for line in fh:
                    s = line.strip()
                    if not s:
                        continue
                    cols = s.split("\t")
                    if len(cols) <= chr_i:
                        continue
                    styles.append(infer_chr_style(cols[chr_i]))
                    break
    if not styles:
        return "no_chr"
    return "chr" if styles.count("chr") >= styles.count("no_chr") else "no_chr"


def load_loci_one_file(path: Path, target_style: str, include_chrom_cores: set[str]) -> dict[str, list[int]]:
    per_chr = defaultdict(list)
    with open_text(path) as fh:
        for line in fh:
            s = line.strip()
            if not s or s.startswith("#"):
                continue
            cols = s.split("\t")
            if len(cols) < 2:
                continue
            if cols[0].lower() in {"chrom", "chromosome", "chr", "contig"}:
                continue
            try:
                pos = int(cols[1])
            except ValueError:
                continue
            core = norm_chr(cols[0])
            if core not in include_chrom_cores:
                continue
            per_chr[convert_chr_style(core, target_style)].append(pos)
    return per_chr


def balanced_sample_loci(per_chr_positions: dict[str, list[int]], max_sites: int | None):
    chroms = sorted(per_chr_positions)
    total_requested = sum(len(v) for v in per_chr_positions.values())
    if max_sites is None or total_requested <= max_sites:
        kept = {c: sorted(v) for c, v in per_chr_positions.items()}
        return kept, total_requested

    remaining = max_sites
    counts = {c: 0 for c in chroms}
    available = {c: len(per_chr_positions[c]) for c in chroms}
    alive = set(chroms)
    while remaining > 0 and alive:
        fair = max(1, remaining // len(alive))
        done = set()
        for c in sorted(alive):
            take = min(fair, available[c] - counts[c])
            counts[c] += take
            remaining -= take
            if counts[c] >= available[c]:
                done.add(c)
            if remaining == 0:
                break
        alive -= done
        if fair == 0:
            break

    kept = {}
    for c in chroms:
        arr = sorted(per_chr_positions[c])
        n = counts[c]
        if n <= 0:
            continue
        if n >= len(arr):
            kept[c] = arr
            continue
        step = len(arr) / n
        idxs = [min(len(arr) - 1, int(i * step)) for i in range(n)]
        kept[c] = [arr[i] for i in idxs]
    return kept, total_requested


def load_loci_dir(loci_dir: Path, target_style: str, include_chroms: list[str], max_sites: int | None):
    include_cores = [norm_chr(c) for c in include_chroms]
    loci_files = list_resource_chrom_files(loci_dir, "loci", include_cores)

    per_chr_positions = defaultdict(list)
    for core in include_cores:
        p = loci_files.get(core)
        if not p:
            continue
        loaded = load_loci_one_file(p, target_style, {core})
        chrom = convert_chr_style(core, target_style)
        per_chr_positions[chrom].extend(loaded.get(chrom, []))
        log(f"loci loaded for {chrom}: {len(loaded.get(chrom, []))}")

    kept_map, total_requested = balanced_sample_loci(per_chr_positions, max_sites)
    for c in sorted(kept_map):
        log(f"loci retained for {c}: {len(kept_map[c])}")

    loci = []
    for chrom in sorted(kept_map):
        for pos in kept_map[chrom]:
            loci.append((chrom, pos))
    per_chr_retained = {c: len(v) for c, v in kept_map.items()}
    return loci, loci_files, total_requested, per_chr_retained


def parse_allele_line(cols: list[str], idx: tuple[int, int, int, int]):
    chr_i, pos_i, ref_i, alt_i = idx
    if len(cols) <= max(idx):
        return None
    try:
        chrom = cols[chr_i]
        pos = int(cols[pos_i])
    except ValueError:
        return None
    ref = cols[ref_i].upper()
    alt = cols[alt_i].upper()
    if not ref or not alt:
        return None
    return norm_chr(chrom), pos, ref, alt


def load_alleles_one_file(path: Path) -> dict:
    allele_map = {}
    with open_text(path) as fh:
        first = fh.readline().strip().split("\t")
        if first and first[0].lower() in {"chromosome", "chrom", "chr", "contig"}:
            hmap = {c.strip().lower(): i for i, c in enumerate(first)}
            idx = (
                find_col(hmap, ["chromosome", "chrom", "chr", "contig"]),
                find_col(hmap, ["position", "pos"]),
                find_col(hmap, ["ref", "ref_allele", "a", "allele_a"]),
                find_col(hmap, ["alt", "alt_allele", "b", "allele_b"]),
            )
            if None in idx:
                raise RuntimeError(f"Cannot parse alleles header in {path}; expected chr/pos/ref/alt columns")
        else:
            idx = (0, 1, 2, 3)
            row = parse_allele_line(first, idx)
            if row:
                allele_map[(row[0], row[1])] = (row[2], row[3])

        for line in fh:
            s = line.strip()
            if not s or s.startswith("#"):
                continue
            row = parse_allele_line(s.split("\t"), idx)
            if not row:
                continue
            allele_map[(row[0], row[1])] = (row[2], row[3])
    return allele_map


def load_alleles_dir(alleles_dir: Path, include_chroms: list[str]):
    include_cores = [norm_chr(c) for c in include_chroms]
    allele_files = list_resource_chrom_files(alleles_dir, "alleles", include_cores)
    allele_map = {}
    for core in include_cores:
        p = allele_files.get(core)
        if not p:
            continue
        m = load_alleles_one_file(p)
        allele_map.update(m)
        log(f"alleles loaded for {core}: {len(m)}")
    if not allele_map:
        raise RuntimeError(f"No allele records parsed from directory: {alleles_dir}")
    return allele_map, allele_files


def load_gc(gc_path: Path, gc_window: str):
    gc_map = {}
    with open_text(gc_path) as fh:
        header = fh.readline().strip().split("\t")
        if not header:
            raise RuntimeError(f"GC file has empty header: {gc_path}")
        hmap = {c.strip().lower(): i for i, c in enumerate(header)}

        # Sarek/ASCAT 的 GC_G1000_hg38.txt 第一列常标记为 Chr，但实际是行ID（例如 1_809641）；
        # 真正染色体与位置通常在后续列，需要按表头定位。
        chr_i = find_col(hmap, ["chromosome", "chrom", "contig", "chrname"])
        if chr_i is None:
            if "chr" in hmap and hmap["chr"] > 0:
                chr_i = hmap["chr"]
            else:
                raise RuntimeError("Cannot locate true chromosome column in GC table")
        pos_i = find_col(hmap, ["position", "pos"])
        if pos_i is None:
            raise RuntimeError("Cannot locate Position column in GC table")

        win_key = gc_window.lower()
        if win_key not in hmap:
            windows = [c for c in header if c.strip().lower() not in {"chr", "chrom", "chromosome", "contig", "position", "pos", "id"}]
            raise RuntimeError(
                f"Requested gc_window '{gc_window}' not found in GC table. Available windows: {', '.join(windows)}"
            )
        gc_i = hmap[win_key]

        for line in fh:
            s = line.strip()
            if not s:
                continue
            cols = s.split("\t")
            if len(cols) <= max(chr_i, pos_i, gc_i):
                continue
            try:
                chrom = norm_chr(cols[chr_i])
                pos = int(cols[pos_i])
                gc = float(cols[gc_i])
            except ValueError:
                continue
            gc_map[(chrom, pos)] = gc

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


def group_loci_by_chr(loci: list[tuple[str, int]]) -> dict[str, list[int]]:
    out = defaultdict(list)
    for chrom, pos in loci:
        out[chrom].append(pos)
    return out


def run_mpileup_counts(tumor_bam: str, normal_bam: str, loci: list[tuple[str, int]], allele_map: dict):
    counts = {}
    loci_by_chr = group_loci_by_chr(loci)

    for chrom in sorted(loci_by_chr):
        with tempfile.NamedTemporaryFile("w", suffix=f".{chrom}.loci", delete=True) as lf:
            for pos in sorted(loci_by_chr[chrom]):
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
                try:
                    pos = int(cols[1])
                    n_dp = int(cols[3])
                    t_dp = int(cols[6])
                except ValueError:
                    continue
                key = (norm_chr(cols[0]), pos)
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
                raise RuntimeError(f"samtools mpileup failed at {chrom} (rc={rc}): {' '.join(cmd)}\n{stderr.strip()}")
    return counts


def has_allelecounter(exe: str) -> bool:
    return shutil.which(exe) is not None


def run_allelecounter_single(exe: str, bam: str, chrom: str, loci_pos: list[int], out_path: Path):
    # 按染色体分批运行，贴近 Sarek/ASCAT 按 chr 的资源组织，也便于定位失败染色体。
    with tempfile.NamedTemporaryFile("w", suffix=f".{chrom}.loci", delete=True) as lf:
        for pos in sorted(loci_pos):
            lf.write(f"{chrom}\t{pos}\n")
        lf.flush()
        cmd = [exe, "-b", bam, "-l", lf.name, "-o", str(out_path), "-q", "20", "-m", "20"]
        r = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        if r.returncode != 0:
            raise RuntimeError(f"alleleCounter failed at {chrom}: {' '.join(cmd)}\n{r.stderr.strip()}")


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
    out = {}
    loci_by_chr = group_loci_by_chr(loci)

    with tempfile.TemporaryDirectory() as td:
        tdp = Path(td)
        for chrom in sorted(loci_by_chr):
            n_out = tdp / f"normal.{chrom}.tsv"
            t_out = tdp / f"tumor.{chrom}.tsv"
            run_allelecounter_single(exe, normal_bam, chrom, loci_by_chr[chrom], n_out)
            run_allelecounter_single(exe, tumor_bam, chrom, loci_by_chr[chrom], t_out)
            n_counts = parse_allelecounter_output(n_out)
            t_counts = parse_allelecounter_output(t_out)
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
    # 轻量近似：logR 对 GC 做一元线性回归并残差中心化，
    # 用于工程化替代，不等价于 ASCAT 官方 ascat.correctLogR()。
    xs = []
    ys = []
    for r in rows:
        gc = gc_map.get((norm_chr(r["chr"]), r["pos"]))
        if gc is None:
            continue
        xs.append(gc)
        ys.append(r["raw_logr"])

    if len(xs) < 100:
        return rows, f"fallback_raw_logr(no_enough_gc_overlap:{len(xs)}/{len(rows)})", True

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


def summarize_per_chr_counts(loci: list[tuple[str, int]]) -> str:
    d = defaultdict(int)
    for c, _ in loci:
        d[c] += 1
    return ";".join(f"{c}:{d[c]}" for c in sorted(d))


def write_run_info(
    outdir: Path,
    sample_id: str,
    cfg: dict,
    loci_n: int,
    rows_n: int,
    baf_method: str,
    logr_method: str,
    fallback: bool,
    bam_chr_style: str,
    resource_chr_style: str,
    effective_chr_style: str,
    n_sites_requested_total: int,
    per_chr_counts_summary: str,
):
    with (outdir / "ascat_prepare.run_info.tsv").open("w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(
            [
                "sample_id",
                "status",
                "root_dir",
                "loci_dir",
                "alleles_dir",
                "gc_path",
                "gc_window",
                "bam_chr_style",
                "resource_chr_style",
                "effective_chr_style",
                "n_sites_requested_total",
                "n_sites_retained_total",
                "per_chr_counts_summary",
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
                str(cfg["loci_dir"]),
                str(cfg["alleles_dir"]),
                str(cfg["gc_path"]),
                str(cfg["gc_window"]),
                bam_chr_style,
                resource_chr_style,
                effective_chr_style,
                n_sites_requested_total,
                loci_n,
                per_chr_counts_summary,
                baf_method,
                logr_method,
                "yes" if fallback else "no",
                cfg["min_normal_depth"],
                cfg["min_tumor_depth"],
                cfg["normal_het_min_baf"],
                cfg["normal_het_max_baf"],
                "Sarek/ASCAT-aligned scaffold; GC correction is lightweight approximation (not ascat.correctLogR)",
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

    tumor_style = detect_bam_chr_style(args.tumor_bam)
    normal_style = detect_bam_chr_style(args.normal_bam)
    if tumor_style != normal_style:
        raise RuntimeError(
            f"Tumor/normal BAM chromosome styles differ: tumor={tumor_style}, normal={normal_style}"
        )
    bam_style = tumor_style

    effective_style = bam_style
    include_chroms = normalize_include_chroms(cfg["include_chrom_cores"], effective_style)
    if not include_chroms:
        raise RuntimeError("No valid chromosomes in include_chroms after normalization")

    log(f"resolved loci_dir={cfg['loci_dir']}")
    log(f"resolved alleles_dir={cfg['alleles_dir']}")
    log(f"resolved gc_path={cfg['gc_path']}")

    loci, loci_files, n_sites_requested_total, per_chr_retained = load_loci_dir(
        cfg["loci_dir"], effective_style, include_chroms, cfg["max_sites"]
    )
    if not loci:
        raise RuntimeError("No loci loaded from loci_dir")

    allele_map, allele_files = load_alleles_dir(cfg["alleles_dir"], include_chroms)
    gc_map = load_gc(cfg["gc_path"], cfg["gc_window"])

    resource_style = detect_resource_chr_style(loci_files, allele_files, cfg["gc_path"])
    log(f"detected bam style={bam_style}, resource style={resource_style}, effective style={effective_style}")
    log(f"loaded loci total(retained): {len(loci)}")
    log(f"loaded alleles total: {len(allele_map)}")
    log(f"loaded gc entries: {len(gc_map)} using window={cfg['gc_window']}")

    locus_keys = {(norm_chr(c), p) for c, p in loci}
    overlap = len(locus_keys & set(allele_map.keys()))
    if overlap < max(1000, int(0.1 * len(locus_keys))):
        log(f"WARNING: low loci/alleles overlap: {overlap}/{len(locus_keys)}")

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

    gc_overlap = sum(1 for r in rows if (norm_chr(r["chr"]), r["pos"]) in gc_map)
    if rows and gc_overlap < max(100, int(0.05 * len(rows))):
        log(f"WARNING: low GC overlap before correction: {gc_overlap}/{len(rows)}")

    rows, logr_method, gc_fallback = linear_gc_correct(rows, gc_map)
    fallback = fallback or gc_fallback

    write_outputs(outdir, rows)
    write_run_info(
        outdir,
        args.sample_id,
        cfg,
        len(loci),
        len(rows),
        baf_method,
        logr_method,
        fallback,
        bam_style,
        resource_style,
        effective_style,
        n_sites_requested_total,
        summarize_per_chr_counts(loci),
    )


if __name__ == "__main__":
    try:
        main()
    except Exception as e:
        print(f"[ascat_prepare][ERROR] {e}", file=sys.stderr, flush=True)
        raise
