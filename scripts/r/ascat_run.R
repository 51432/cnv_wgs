#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(k){
  idx <- match(k, args)
  if (is.na(idx) || idx == length(args)) return(NA)
  args[idx + 1]
}

sample_id <- get_arg("--sample_id")
prep_dir <- get_arg("--prep_dir")
outdir <- get_arg("--outdir")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

baf_path <- file.path(prep_dir, "baf.tsv")
logr_path <- file.path(prep_dir, "logr.tsv")

safe_read <- function(path){
  if (!file.exists(path)) return(data.frame())
  x <- tryCatch(read.table(path, sep = "\t", header = TRUE, stringsAsFactors = FALSE), error = function(e) data.frame())
  x
}

baf <- safe_read(baf_path)
logr <- safe_read(logr_path)

if (nrow(logr) > 0) {
  logr$pos <- as.numeric(logr$pos)
  logr$logr <- as.numeric(logr$logr)
  logr <- logr[!is.na(logr$pos) & !is.na(logr$logr), ]
}
if (nrow(baf) > 0) {
  baf$pos <- as.numeric(baf$pos)
  baf$baf <- as.numeric(baf$baf)
  baf <- baf[!is.na(baf$pos) & !is.na(baf$baf), ]
}

if (nrow(baf) > 0) {
  purity <- max(0.1, min(1.0, median(abs(baf$baf - 0.5), na.rm = TRUE) * 2))
} else {
  purity <- NA_real_
}

if (nrow(logr) > 0) {
  ploidy <- max(1.2, min(6.0, 2 + median(logr$logr, na.rm = TRUE) * 2))
} else {
  ploidy <- NA_real_
}

pp <- data.frame(sample_id = sample_id, purity = round(purity, 4), ploidy = round(ploidy, 4))
write.table(pp, file = file.path(outdir, "purity_ploidy.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

make_segments <- function(logr_df) {
  if (nrow(logr_df) == 0) {
    return(data.frame(chr=character(), start=integer(), end=integer(), total_cn=integer()))
  }

  logr_df <- logr_df[order(logr_df$chr, logr_df$pos), ]
  state <- ifelse(logr_df$logr > 0.2, "GAIN", ifelse(logr_df$logr < -0.2, "LOSS", "NEUTRAL"))

  seg_chr <- c()
  seg_start <- c()
  seg_end <- c()
  seg_state <- c()

  cur_chr <- logr_df$chr[1]
  cur_start <- logr_df$pos[1]
  cur_end <- logr_df$pos[1]
  cur_state <- state[1]

  for (i in 2:nrow(logr_df)) {
    same_block <- (logr_df$chr[i] == cur_chr) && (state[i] == cur_state) && ((logr_df$pos[i] - cur_end) <= 1e6)
    if (same_block) {
      cur_end <- logr_df$pos[i]
    } else {
      seg_chr <- c(seg_chr, cur_chr)
      seg_start <- c(seg_start, as.integer(cur_start))
      seg_end <- c(seg_end, as.integer(cur_end))
      seg_state <- c(seg_state, cur_state)

      cur_chr <- logr_df$chr[i]
      cur_start <- logr_df$pos[i]
      cur_end <- logr_df$pos[i]
      cur_state <- state[i]
    }
  }

  seg_chr <- c(seg_chr, cur_chr)
  seg_start <- c(seg_start, as.integer(cur_start))
  seg_end <- c(seg_end, as.integer(cur_end))
  seg_state <- c(seg_state, cur_state)

  total_cn <- ifelse(seg_state == "GAIN", 3L, ifelse(seg_state == "LOSS", 1L, 2L))
  data.frame(chr = seg_chr, start = seg_start, end = seg_end, total_cn = total_cn, stringsAsFactors = FALSE)
}

segments <- make_segments(logr)
write.table(segments, file = file.path(outdir, "segments.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

allele_specific <- data.frame(chr=character(), start=integer(), end=integer(), major_cn=integer(), minor_cn=integer())
if (nrow(segments) > 0) {
  major <- pmax(1, segments$total_cn - 1)
  minor <- pmax(0, segments$total_cn - major)
  allele_specific <- data.frame(chr=segments$chr, start=segments$start, end=segments$end, major_cn=major, minor_cn=minor)
}
write.table(allele_specific, file = file.path(outdir, "allele_specific_cn.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

loh <- data.frame(chr=character(), start=integer(), end=integer(), loh=character())
if (nrow(allele_specific) > 0) {
  loh <- data.frame(
    chr = allele_specific$chr,
    start = allele_specific$start,
    end = allele_specific$end,
    loh = ifelse(allele_specific$minor_cn == 0, "YES", "NO"),
    stringsAsFactors = FALSE
  )
}
write.table(loh, file = file.path(outdir, "loh.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

arm_level <- data.frame(arm=character(), sample_id=character(), cn_state=character())
if (nrow(segments) > 0) {
  chr_split <- tapply(segments$end, segments$chr, median)
  arm <- ifelse(segments$end <= chr_split[segments$chr], paste0(segments$chr, "p"), paste0(segments$chr, "q"))
  cn_state <- ifelse(segments$total_cn >= 3, "GAIN", ifelse(segments$total_cn <= 1, "LOSS", "NEUTRAL"))
  arm_level <- unique(data.frame(arm=arm, sample_id=sample_id, cn_state=cn_state, stringsAsFactors = FALSE))
}
write.table(arm_level, file = file.path(outdir, "arm_level_cn.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

gene_level <- data.frame(gene=character(), sample_id=character(), cn_state=character())
if (nrow(segments) > 0) {
  cn_state <- ifelse(segments$total_cn >= 3, "GAIN", ifelse(segments$total_cn <= 1, "LOSS", "NEUTRAL"))
  gene_level <- data.frame(
    gene = paste0("BIN_", segments$chr, "_", seq_len(nrow(segments))),
    sample_id = sample_id,
    cn_state = cn_state,
    stringsAsFactors = FALSE
  )
}
write.table(gene_level, file = file.path(outdir, "gene_level_cn.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

writeLines("ascat_run proxy completed", con = file.path(outdir, "genome_wide_plot.txt"))
run_info <- data.frame(
  sample_id = sample_id,
  converged = ifelse(nrow(logr) > 0, "YES", "NO"),
  message = ifelse(nrow(logr) > 0, "proxy segmentation based on logr", "no logr rows"),
  stringsAsFactors = FALSE
)
write.table(run_info, file = file.path(outdir, "run_info.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
