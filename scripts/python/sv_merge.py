#!/usr/bin/env python3
import argparse
from pathlib import Path


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--sample-id', required=True)
    ap.add_argument('--postfilter-dir', required=True)
    ap.add_argument('--outdir', required=True)
    args = ap.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    merged = outdir / 'merged.sv.tsv'

    # skeleton: union + source label placeholder
    merged.write_text('sample_id\tsv_id\tsource\n%s\tNA\tBOTH\n' % args.sample_id)


if __name__ == '__main__':
    main()
