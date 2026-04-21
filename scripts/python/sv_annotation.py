#!/usr/bin/env python3
import argparse
from pathlib import Path


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
    (outdir / 'annotated.sv.tsv').write_text(
        'sample_id\tsv_id\tgene\tcytoband\tsv_type\tnear_cn_breakpoint\tnear_hpv\n'
        f'{args.sample_id}\tNA\tNA\tNA\tNA\tNA\tNA\n'
    )


if __name__ == '__main__':
    main()
