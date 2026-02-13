#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(tidyverse)
  library(yaml)
})

args <- commandArgs(trailingOnly = TRUE)
config <- read_yaml(args[1])

# --------- helpers ----------
read_samples <- function(samples_path){
  lines <- readLines(samples_path)
  lines <- lines[grepl("\\S", lines)]
  split_lines <- strsplit(lines, "[[:space:]]+")
  header <- split_lines[[1]]
  rows <- split_lines[-1]

  df <- as.data.frame(do.call(rbind, rows), stringsAsFactors = FALSE)
  colnames(df) <- header

  if (!("sample" %in% colnames(df))) stop("samples.tsv must have a column called 'sample'")
  df
}

count_diamond_passing <- function(diamond_path, min_cov, min_pident){
  if (!file.exists(diamond_path)) return(NA_integer_)
  if (file.size(diamond_path) == 0) return(0L)

  dmd <- readr::read_tsv(
    diamond_path,
    col_names = c("qseqid","sseqid","pident","length","slen","qlen","gaps","evalue"),
    col_types = readr::cols(
      qseqid = readr::col_character(),
      sseqid = readr::col_character(),
      pident = readr::col_double(),
      length = readr::col_integer(),
      slen   = readr::col_integer(),
      qlen   = readr::col_integer(),
      gaps   = readr::col_integer(),
      evalue = readr::col_double()
    )
  )

  if (nrow(dmd) == 0) return(0L)

  dmd <- dmd %>% mutate(cov = 100 * (length - gaps) / slen)
  sum(dmd$cov >= min_cov & dmd$pident >= min_pident, na.rm = TRUE)
}

count_rows_minus_header <- function(path){
  if (!file.exists(path)) return(NA_integer_)
  n <- tryCatch(nrow(readr::read_tsv(path, show_col_types = FALSE)), error = function(e) NA_integer_)
  if (is.na(n)) return(NA_integer_)
  n
}

# --------- inputs from config ----------
if (is.null(config$samples_tsv) || config$samples_tsv == "") stop("samples_tsv is missing in run.yaml")

tmp_dir  <- config$tmp_dir
out_base <- config$cf_classification_outdir

min_cov    <- config$min_cov
min_pident <- config$min_pident
min_score  <- config$cf_classify_min_score_fraction

# diamond dir (same as your other scripts)
diamond_dir <- file.path(tmp_dir, "diamond_aligned_amino_cf_db")

samples_df  <- read_samples(config$samples_tsv)
genomes     <- unique(samples_df$sample)

# --------- per-genome summary ----------
summary_df <- purrr::map_dfr(genomes, function(g){

  diamond_path <- file.path(diamond_dir, paste0("diamond_align_", g, ".tsv"))
  n_diamond_passing <- count_diamond_passing(diamond_path, min_cov, min_pident)

  cf_path <- file.path(out_base, g, paste0("cf_classification_min_classify_score_fraction_", min_score, ".tsv"))
  cf_rows <- count_rows_minus_header(cf_path)
  cf_operons_called <- ifelse(is.na(cf_rows), NA_integer_, cf_rows)

  diamond_any_passing <- ifelse(is.na(n_diamond_passing), NA, n_diamond_passing > 0)
  cf_called <- ifelse(is.na(cf_operons_called), NA, cf_operons_called > 0)

  reason_cf_negative <- case_when(
    is.na(diamond_any_passing) | is.na(cf_called) ~ NA_character_,
    cf_called ~ "cf_called",
    !diamond_any_passing ~ "no_diamond_hits",
    diamond_any_passing & !cf_called ~ "diamond_hits_but_no_operon_call",
    TRUE ~ NA_character_
  )

  tibble(
    genome_id = g,
    min_cov = min_cov,
    min_pident = min_pident,
    cf_classify_min_score_fraction = min_score,
    diamond_hits_passing = n_diamond_passing,
    diamond_any_passing = diamond_any_passing,
    cf_operons_called = cf_operons_called,
    cf_called = cf_called,
    reason_cf_negative = reason_cf_negative,
    diamond_file = diamond_path,
    cf_classification_file = cf_path
  )
})

# write per-genome table
dir.create(out_base, recursive = TRUE, showWarnings = FALSE)
out_path <- file.path(out_base, "ecgc_run_summary.tsv")
readr::write_tsv(summary_df, out_path)

# also print a short run summary to stdout (goes into snakemake logs)
n_total <- nrow(summary_df)

n_dpos <- sum(summary_df$diamond_any_passing %in% TRUE, na.rm = TRUE)
n_dneg <- sum(summary_df$diamond_any_passing %in% FALSE, na.rm = TRUE)

n_cfpos <- sum(summary_df$cf_called %in% TRUE, na.rm = TRUE)
n_cfneg <- sum(summary_df$cf_called %in% FALSE, na.rm = TRUE)

n_dpos_cfneg <- sum(summary_df$reason_cf_negative == "diamond_hits_but_no_operon_call", na.rm = TRUE)

message("Run summary:")
message("  genomes_total = ", n_total)
message("  diamond_positive = ", n_dpos, " | diamond_negative = ", n_dneg)
message("  cf_called = ", n_cfpos, " | cf_negative = ", n_cfneg)
message("  diamond_positive_but_cf_negative = ", n_dpos_cfneg)
message("Wrote: ", out_path)