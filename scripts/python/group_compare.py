#!/usr/bin/env python3
import argparse
from pathlib import Path
import csv


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--group-tsv', required=True)
    ap.add_argument('--cohort-dir', required=True)
    ap.add_argument('--outdir', required=True)
    args = ap.parse_args()
    outdir = Path(args.outdir); outdir.mkdir(parents=True, exist_ok=True)
    with open(args.group_tsv) as fh:
        rd = csv.reader(fh, delimiter='\t')
        header = next(rd)
    for col in header[1:]:
        (outdir / f'{col}.summary.tsv').write_text(f'group\tn\tmetric\n')
        (outdir / f'{col}.plot.txt').write_text('placeholder plot\n')

if __name__ == '__main__':
    main()
