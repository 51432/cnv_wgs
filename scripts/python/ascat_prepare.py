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
    (outdir / 'baf.tsv').write_text('chr\tpos\tbaf\n')
    (outdir / 'logr.tsv').write_text('chr\tpos\tlogr\n')
    (outdir / 'ascat_prepare.run_info.tsv').write_text('sample_id\tstatus\n%s\tplaceholder\n' % args.sample_id)

if __name__ == '__main__':
    main()
