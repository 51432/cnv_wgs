#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(k){
  idx <- match(k, args)
  if (is.na(idx)) return(NA)
  args[idx+1]
}
outdir <- get_arg("--outdir")
sample_id <- get_arg("--sample_id")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
write.table(data.frame(sample_id=sample_id,purity="NA",ploidy="NA"),
            file=file.path(outdir,"purity_ploidy.tsv"),sep="\t",quote=FALSE,row.names=FALSE)
write.table(data.frame(chr=character(),start=integer(),end=integer(),total_cn=integer()),
            file=file.path(outdir,"segments.tsv"),sep="\t",quote=FALSE,row.names=FALSE)
write.table(data.frame(chr=character(),start=integer(),end=integer(),major_cn=integer(),minor_cn=integer()),
            file=file.path(outdir,"allele_specific_cn.tsv"),sep="\t",quote=FALSE,row.names=FALSE)
write.table(data.frame(chr=character(),start=integer(),end=integer(),loh=character()),
            file=file.path(outdir,"loh.tsv"),sep="\t",quote=FALSE,row.names=FALSE)
write.table(data.frame(arm=character(),sample_id=character(),cn_state=character()),
            file=file.path(outdir,"arm_level_cn.tsv"),sep="\t",quote=FALSE,row.names=FALSE)
write.table(data.frame(gene=character(),sample_id=character(),cn_state=character()),
            file=file.path(outdir,"gene_level_cn.tsv"),sep="\t",quote=FALSE,row.names=FALSE)
writeLines("placeholder plot", con=file.path(outdir, "genome_wide_plot.txt"))
write.table(data.frame(sample_id=sample_id,converged="NO",message="skeleton"),
            file=file.path(outdir,"run_info.tsv"),sep="\t",quote=FALSE,row.names=FALSE)
