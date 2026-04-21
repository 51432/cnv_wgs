#!/usr/bin/env python3
import argparse
from pathlib import Path

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--sample-id', required=True)
    ap.add_argument('--manta-vcf', required=True)
    ap.add_argument('--gridss-vcf', required=True)
    ap.add_argument('--outdir', required=True)
    args = ap.parse_args()
    outdir = Path(args.outdir); outdir.mkdir(parents=True, exist_ok=True)
    (outdir/'manta.postfilter.tsv').write_text('sample_id\tsv_id\tsource\tpe\tsr\tqual\tnormal_support\n')
    (outdir/'gridss.postfilter.tsv').write_text('sample_id\tsv_id\tsource\tpe\tsr\tqual\tnormal_support\n')

if __name__ == '__main__':
    main()
