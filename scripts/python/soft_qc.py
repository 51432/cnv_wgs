#!/usr/bin/env python3
import argparse
from pathlib import Path

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--sample-id', required=True)
    ap.add_argument('--tumor-bam', required=True)
    ap.add_argument('--normal-bam', required=True)
    ap.add_argument('--outdir', required=True)
    args = ap.parse_args()
    outdir = Path(args.outdir); outdir.mkdir(parents=True, exist_ok=True)
    (outdir / 'soft_qc.summary.tsv').write_text(
        'sample_id\ttotal_reads_tumor\tmapped_rate_tumor\tduplicate_rate_tumor\tmean_cov_tumor\tmean_cov_normal\ttn_cov_ratio\tqc_flag\n'
        f'{args.sample_id}\tNA\tNA\tNA\tNA\tNA\tNA\tREVIEW\n'
    )

if __name__ == '__main__':
    main()
