#!/usr/bin/env python3
import argparse
from pathlib import Path


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--ascat-dir', required=True)
    ap.add_argument('--sv-annot-dir', required=True)
    ap.add_argument('--outdir', required=True)
    args = ap.parse_args()
    outdir = Path(args.outdir); outdir.mkdir(parents=True, exist_ok=True)
    (outdir/'cohort.purity_ploidy.tsv').write_text('sample_id\tpurity\tploidy\n')
    (outdir/'cohort.arm_level_cn.matrix.tsv').write_text('arm\n')
    (outdir/'cohort.gene_level_cn.matrix.tsv').write_text('gene\n')
    (outdir/'cohort.loh.matrix.tsv').write_text('region\n')
    (outdir/'cohort.sv_burden.tsv').write_text('sample_id\tsv_count\n')
    (outdir/'cohort.sv_event_table.tsv').write_text('sample_id\tsv_id\tsource\n')

if __name__ == '__main__':
    main()
