#!/usr/bin/env python3
import argparse
from pathlib import Path
import re


def parse_breakpoint_line(line: str):
    m = re.search(r'(chr[^:\s]+):(-?\d+)', line)
    if not m:
        return None, None
    return m.group(1), int(m.group(2))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--hpv-breakpoints', required=True)
    ap.add_argument('--ascat-dir', required=True)
    ap.add_argument('--sv-annot-dir', required=True)
    ap.add_argument('--outdir', required=True)
    args = ap.parse_args()
    outdir = Path(args.outdir); outdir.mkdir(parents=True, exist_ok=True)

    cn = outdir / 'hpv_cn_link.per_sample.tsv'
    sv = outdir / 'hpv_sv_link.per_sample.tsv'
    summary = outdir / 'hpv_link.cohort_summary.tsv'

    with open(args.hpv_breakpoints) as fh, cn.open('w') as fcn, sv.open('w') as fsv:
        fcn.write('sample_id\thost_chr\thost_pos\twindow_bp\tn_cn_breakpoints\n')
        fsv.write('sample_id\thost_chr\thost_pos\twindow_bp\tn_sv_breakpoints\n')
        for idx, line in enumerate(fh, start=1):
            line = line.strip()
            if not line:
                continue
            ch, pos = parse_breakpoint_line(line)
            sample_id = f'HPV_EVENT_{idx}'
            if ch is None:
                continue
            for w in (10000, 100000, 1000000):
                fcn.write(f'{sample_id}\t{ch}\t{pos}\t{w}\tNA\n')
                fsv.write(f'{sample_id}\t{ch}\t{pos}\t{w}\tNA\n')

    summary.write_text('metric\tvalue\nparsed_events\tNA\n')

if __name__ == '__main__':
    main()
