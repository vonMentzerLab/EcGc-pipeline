#!/usr/bin/env Rscript

library(tidyverse)
library(plyranges)
library(BSgenome)
library(Biostrings)
library(Rsamtools)
library(rtracklayer)
library(foreach)
library(doParallel)
library(yaml)

load_fa_seqs <- function(fa_file){
  Biostrings::getSeq(Rsamtools::FaFile(fa_file))
}

make_bakta_contignames <- function(dnastringset){
  names(dnastringset) <- paste0("contig_", seq_along(dnastringset))
  dnastringset
}

args <- commandArgs(trailingOnly = TRUE)
config <- read_yaml(args[1])

# Backwards/forwards compatible config keys
if (is.null(config$genome_link_dir) && !is.null(config$genome_seqs_dir)) {
  config$genome_link_dir <- config$genome_seqs_dir
}
if (is.null(config$gff3_suffix)) config$gff3_suffix <- ".gff3"

# output dirs from config
nucl_dir        <- file.path(config$cds_features_outdir, "cds_nucl_seqs")
prot_dir        <- file.path(config$cds_features_outdir, "cds_amino_seqs")
cds_coords_dir  <- file.path(config$cds_features_outdir, "cds_coords")

dir.create(nucl_dir,       recursive = TRUE, showWarnings = FALSE)
dir.create(prot_dir,       recursive = TRUE, showWarnings = FALSE)
dir.create(cds_coords_dir, recursive = TRUE, showWarnings = FALSE)

########################################
# read samples.tsv robustly (tabs OR spaces)
########################################
samples_path <- "config/samples.tsv"
lines <- readLines(samples_path)
lines <- lines[grepl("\\S", lines)]
split_lines <- strsplit(lines, "[[:space:]]+")
header <- split_lines[[1]]
rows <- split_lines[-1]

samples_df <- as.data.frame(do.call(rbind, rows), stringsAsFactors = FALSE)
colnames(samples_df) <- header

if (!("sample" %in% colnames(samples_df))) {
  stop("samples.tsv must have a column called 'sample'")
}
if (!("fasta" %in% colnames(samples_df))) {
  stop("samples.tsv must have a column called 'fasta'")
}

genome_ids <- samples_df$sample
message("Samples detected: ", paste(genome_ids, collapse = ", "))

########################################
# parallel extraction loop
########################################
cores  <- config$cores
cl <- makeCluster(cores)
registerDoParallel(cl)

bacterial_gen_code <- getGeneticCode("11")

silence <- foreach(
  genome_id = genome_ids,
  .inorder = FALSE,
  .packages = c("tidyverse", "plyranges", "BSgenome", "Biostrings", "Rsamtools", "rtracklayer")
) %dopar% {

  # paths for THIS sample produced by Snakefile rules
  gff_path <- file.path(config$bakta_outdir, genome_id, paste0(genome_id, config$gff3_suffix))

  # FASTA can be either .fa.gz or .fa (we prefer .fa.gz if both exist)
  fa_gz <- file.path(config$genome_link_dir, paste0(genome_id, ".fa.gz"))
  fa_fa <- file.path(config$genome_link_dir, paste0(genome_id, ".fa"))

  fa_path <- if (file.exists(fa_gz)) fa_gz else fa_fa
  if (!file.exists(fa_path)) {
    stop("FASTA not found for sample ", genome_id,
       ". Expected one of: ", fa_gz, " or ", fa_fa)
  }
  
  # load CDS features from Bakta GFF
  cds <- rtracklayer::import.gff3(gff_path) %>%
    dplyr::filter(type == "CDS") %>%
    dplyr::select(phase, ID) %>%
    dplyr::mutate(
      genome_id = genome_id,
      my_cds_id = paste0(genome_id, ".", ID)
    )

  # adjust width if phase != 0
  if (any(cds$phase != 0, na.rm = TRUE)) {
    cds <- cds %>%
      anchor_3p() %>%
      mutate(width = width - phase)
  }

  # load contig sequences from the symlinked FASTA Bakta actually annotated
  contig_seqs <- load_fa_seqs(fa_path) %>%
    make_bakta_contignames()

  # extract CDS nucleotide seqs and translate
  nucl_seqs <- BSgenome::getSeq(contig_seqs, cds)
  names(nucl_seqs) <- cds$my_cds_id

  aa_seqs <- Biostrings::translate(
    nucl_seqs,
    genetic.code = bacterial_gen_code,
    if.fuzzy.codon = "solve"
  )

  # write outputs for THIS genome_id
  Biostrings::writeXStringSet(
    nucl_seqs,
    file.path(nucl_dir, paste0("cds_nucl_seqs_", genome_id, ".fa.gz")),
    compress = TRUE
  )

  Biostrings::writeXStringSet(
    aa_seqs,
    file.path(prot_dir, paste0("cds_amino_seqs_", genome_id, ".fa.gz")),
    compress = TRUE
  )

  cds %>%
    as_tibble() %>%
    readr::write_tsv(
      file.path(cds_coords_dir, paste0("cds_features_", genome_id, ".tsv.gz"))
    )

  NULL
}

message("CDS extraction finished.")