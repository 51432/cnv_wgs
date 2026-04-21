#!/usr/bin/env python3
import argparse
import csv
from pathlib import Path


def read_tsv(path):
    p = Path(path)
    if not p.exists() or p.stat().st_size == 0:
        return []
    with p.open() as f:
        return list(csv.DictReader(f, delimiter='\t'))


def read_segments(ascat_dir):
    segs = []
    for r in read_tsv(Path(ascat_dir) / 'segments.tsv'):
        try:
            segs.append((r.get('chr', ''), int(r.get('start', '0')), int(r.get('end', '0'))))
        except ValueError:
            continue
    return segs


def read_hpv_positions(hpv_file):
    if not hpv_file:
        return []
    p = Path(hpv_file)
    if not p.exists():
        return []

    out = []
    with p.open() as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            for token in line.replace('|', ' ').split():
                if ':' not in token:
                    continue
                chrom, pos = token.split(':', 1)
                chrom = chrom.replace('chr', '')
                try:
                    ipos = abs(int(pos))
                    out.append((f'chr{chrom}', ipos))
                except ValueError:
                    continue
    return out


def resolve_chr_pos(row):
    chr_keys = ['chr', 'chrom', 'chrom1', 'chr1']
    pos_keys = ['pos', 'start', 'position', 'pos1', 'start1']

    chrom = None
    pos = None
    for k in chr_keys:
        if row.get(k):
            chrom = row[k]
            break
    for k in pos_keys:
        if row.get(k):
            try:
                pos = int(float(row[k]))
                break
            except ValueError:
                pass

    if chrom and not chrom.startswith('chr'):
        chrom = f'chr{chrom}'
    return chrom, pos


def near_cn_breakpoint(chrom, pos, segs, window=100000):
    if not chrom or pos is None:
        return 'NA'
    for c, s, e in segs:
        if c != chrom:
            continue
        if abs(pos - s) <= window or abs(pos - e) <= window:
            return 'YES'
    return 'NO'


def near_hpv(chrom, pos, hpv_sites, window=100000):
    if not chrom or pos is None:
        return 'NA'
    for c, p in hpv_sites:
        if c == chrom and abs(pos - p) <= window:
            return 'YES'
    return 'NO'


def cytoband_proxy(chrom, pos):
    if not chrom or pos is None:
        return 'NA'
    arm = 'p' if pos < 50000000 else 'q'
    return f'{chrom}{arm}'


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--sample-id', required=True)
    ap.add_argument('--merged-tsv', required=True)
    ap.add_argument('--ascat-dir', required=True)
    ap.add_argument('--hpv-breakpoints', required=False, default='')
    ap.add_argument('--outdir', required=True)
    args = ap.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    merged_rows = read_tsv(args.merged_tsv)
    segs = read_segments(args.ascat_dir)
    hpv_sites = read_hpv_positions(args.hpv_breakpoints)

    out_path = outdir / 'annotated.sv.tsv'
    with out_path.open('w', newline='') as f:
        writer = csv.writer(f, delimiter='\t')
        writer.writerow(['sample_id', 'sv_id', 'gene', 'cytoband', 'sv_type', 'near_cn_breakpoint', 'near_hpv'])

        if not merged_rows:
            return

        for i, row in enumerate(merged_rows, start=1):
            chrom, pos = resolve_chr_pos(row)
            sv_id = row.get('sv_id') or row.get('id') or f'{args.sample_id}_SV{i}'
            sv_type = row.get('sv_type') or row.get('type') or 'NA'
            gene = row.get('gene') or 'NA'
            writer.writerow([
                args.sample_id,
                sv_id,
                gene,
                cytoband_proxy(chrom, pos),
                sv_type,
                near_cn_breakpoint(chrom, pos, segs),
                near_hpv(chrom, pos, hpv_sites),
            ])


if __name__ == '__main__':
    main()
