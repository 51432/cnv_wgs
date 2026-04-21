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
    ap.add_argument('--sample-id', required=True)
    ap.add_argument('--hpv-breakpoints', required=True)
    ap.add_argument('--ascat-dir', required=True)
    ap.add_argument('--sv-annot-tsv', required=True)
    ap.add_argument('--outdir', required=True)
    args = ap.parse_args()
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    cn = outdir / 'hpv_cn_link.tsv'
    sv = outdir / 'hpv_sv_link.tsv'
    summary = outdir / 'hpv_link.summary.tsv'

    cn_rows = []
    sv_rows = []
    with open(args.hpv_breakpoints) as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            ch, pos = parse_breakpoint_line(line)
            if ch is None:
                continue
            for w in (10000, 100000, 1000000):
                cn_rows.append(f'{args.sample_id}\t{ch}\t{pos}\t{w}\tNA')
                sv_rows.append(f'{args.sample_id}\t{ch}\t{pos}\t{w}\tNA')

    cn.write_text('sample_id\thost_chr\thost_pos\twindow_bp\tn_cn_breakpoints\n' + ('\n'.join(cn_rows) + '\n' if cn_rows else ''))
    sv.write_text('sample_id\thost_chr\thost_pos\twindow_bp\tn_sv_breakpoints\n' + ('\n'.join(sv_rows) + '\n' if sv_rows else ''))
    summary.write_text('sample_id\tparsed_events\n%s\t%d\n' % (args.sample_id, len(cn_rows) // 3 if cn_rows else 0))


if __name__ == '__main__':
    main()
