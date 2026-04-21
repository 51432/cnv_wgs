#!/usr/bin/env python3
import argparse
from pathlib import Path

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--sample-id', required=True)
    ap.add_argument('--tumor-bam', required=True)
    ap.add_argument('--normal-bam', required=True)
    ap.add_argument('--config', required=True)
    ap.add_argument('--outdir', required=True)
    args = ap.parse_args()
    outdir = Path(args.outdir); outdir.mkdir(parents=True, exist_ok=True)
    (outdir/'somaticSV.vcf.gz').write_text('##placeholder manta vcf\n')
    (outdir/'manta.summary.tsv').write_text('sample_id\tn_events\n%s\t0\n' % args.sample_id)

if __name__ == '__main__':
    main()
