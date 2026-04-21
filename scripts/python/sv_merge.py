#!/usr/bin/env python3
import argparse
from pathlib import Path
import csv

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--pairbam', required=True)
    ap.add_argument('--postfilter-dir', required=True)
    ap.add_argument('--outdir', required=True)
    args = ap.parse_args()
    outdir = Path(args.outdir); outdir.mkdir(parents=True, exist_ok=True)
    merged = outdir / 'merged.sv.tsv'
    with open(args.pairbam) as fh, merged.open('w') as fw:
      rd = csv.DictReader(fh, delimiter='\t')
      fw.write('sample_id\tsv_id\tsource\n')
      for row in rd:
        fw.write(f"{row['sample_id']}\tNA\tBOTH\n")

if __name__ == '__main__':
    main()
