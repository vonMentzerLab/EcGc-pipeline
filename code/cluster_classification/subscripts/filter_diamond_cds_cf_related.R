library(tidyverse)
library(Biostrings)
library(doParallel)
library(foreach)
library(yaml)
library(Rsamtools)

load_fa_seqs <- function(fa_file){
  Biostrings::getSeq(Rsamtools::FaFile(fa_file))
}

myalign <- function(seqs, seq2){
  pairwiseAlignment(seqs, seq2, substitutionMatrix = 'BLOSUM62', scoreOnly = TRUE, gapOpening = 10, gapExtension = 4)
}

myalign2 <- function(seqs, seq2){
  submat = 'BLOSUM62'
  gapopen = 10
  gapext = 4

  alignment <- pairwiseAlignment(seqs, seq2, substitutionMatrix = submat, scoreOnly = FALSE, gapOpening = gapopen, gapExtension = gapext)
  tibble(
    seq1_name = names(seqs),
    seq2_name= names(seq2),
    score = score(alignment),
    max_score = pairwiseAlignment(seq2, seq2, substitutionMatrix = submat, scoreOnly = TRUE, gapOpening = gapopen, gapExtension = gapext),
    pid = pid(alignment))
}

bacterial_gen_code <- getGeneticCode('11')

args <- commandArgs(trailingOnly = TRUE)
config <- read_yaml(args[1])

# map new config keys -> internal variables this script expects
# tmp_dir replaces old gene_search_dir
# cds_features_dir does not exist in YAML, but we know where cds outputs go
config$gene_search_dir <- config$tmp_dir
config$cds_features_dir <- config$cds_features_outdir
config$cores <- config$classification_threads

# new reference key names
config$clustered_genes_reference_db_nucl_seq <- config$cf_reference_nt_db
config$clustered_genes_reference_db_aa_seq <- config$cf_reference_aa_db

# build working subdirs under tmp_dir instead of /
config$diamond_prot_dir <- file.path(config$gene_search_dir, "diamond_aligned_amino_cf_db")
config$cf_filtered_cds_coords_dir  <- file.path(config$gene_search_dir, "cf_cds_coords")
config$global_alignment_scores_dir  <- file.path(config$gene_search_dir, "cds_global_alignment_scores")

config$cds_coords_dir <- file.path(config$cds_features_dir, "cds_coords")
config$nucl_seqs_dir <- file.path(config$cds_features_dir, "cds_nucl_seqs")

# make sure dirs exist
dir.create(config$cf_filtered_cds_coords_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(config$global_alignment_scores_dir, recursive = TRUE, showWarnings = FALSE)

files <- list.files(config$diamond_prot_dir)

# load AA reference sequences directly from provided translated DB
CF_reference_aa_seqs <- Biostrings::readAAStringSet(config$clustered_genes_reference_db_aa_seq) 

cores <- config$cores
cl <- makeCluster(cores)
registerDoParallel(cl)

silence <- foreach(
  file = files,
  .packages = c("tidyverse", "Biostrings", "Rsamtools")
) %dopar% {

  dmd_path <- file.path(config$diamond_prot_dir, file)

  # skip if file is empty or 0 bytes
  if (file.size(dmd_path) == 0) {
    return(NULL)
  }

  # read diamond table as exactly 8 tab-delimited columns
  dmd <- readr::read_tsv(
    dmd_path,
    col_names = c(
      "qseqid","sseqid","pident","length","slen","qlen","gaps","evalue"
    ),
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

  # if after reading there are no rows, bail out gracefully
  if (nrow(dmd) == 0) {
    return(NULL)
  }

  dmd <- dmd %>%
    mutate(cov = 100 * (length - gaps) / slen)

  # Filter cf diamond hits with minimum pident and coverage
  dmd_qseqid_filtered <- dmd %>% 
    filter(
      cov >= config$min_cov,
      pident >= config$min_pident
    ) %>% 
    dplyr::select(qseqid) %>% 
    distinct()

  # foreach "continue" if there are no CDSs remaining after filtering
  if (nrow(dmd_qseqid_filtered) == 0) {
    return(NULL)
  }

  # Get genome_id so that we can load the cds features file
  genome_id <- dmd_qseqid_filtered %>% 
    tidyr::separate(qseqid, into = c('genome_id', 'cds_id'), sep = '\\.', remove = FALSE) %>% 
    dplyr::pull(genome_id) %>% 
    unique()

  cds_features <- readr::read_tsv(
    file.path(config$cds_coords_dir, paste0("cds_features_", genome_id, ".tsv.gz")),
    col_types = "ciiiciccc"
  )

  # Get first start and last end of diamond hit genes
  cf_cds_start_end <- cds_features %>% 
    dplyr::filter(my_cds_id %in% dmd_qseqid_filtered$qseqid) %>% 
    dplyr::group_by(seqnames) %>% 
    dplyr::summarise(
      first_start = min(start),
      last_end = max(end),
      .groups = "drop"
    )

  # keep any CDS near those hits within max_dist
  cds_filtered <- cds_features %>%
    dplyr::inner_join(cf_cds_start_end, by = "seqnames") %>%
    dplyr::left_join(
      dmd_qseqid_filtered %>%
        dplyr::rename(my_cds_id = qseqid) %>%
        dplyr::mutate(diamond_cf_candidate = TRUE),
      by = "my_cds_id"
    ) %>%
    tidyr::replace_na(list(diamond_cf_candidate = FALSE)) %>%
    dplyr::group_by(seqnames) %>%
    dplyr::filter(first_start - end <= config$max_dist & start - last_end <= config$max_dist) %>%
    dplyr::ungroup() %>%
    dplyr::select(-first_start, -last_end)

  # write coords around CF for this genome
  readr::write_tsv(
    cds_filtered,
    file.path(
      config$cf_filtered_cds_coords_dir,
      paste0("cds_around_cf_", genome_id, ".tsv.gz")
    )
  )

  # get their nucleotide seqs (per-genome file from earlier pipeline step)
  nucl_seqs <- load_fa_seqs(
    file.path(
      config$nucl_seqs_dir,
      paste0("cds_nucl_seqs_", genome_id, ".fa.gz")
    )
  )
  nucl_seqs <- nucl_seqs[cds_filtered$my_cds_id]

  aa_seqs <- Biostrings::translate(
    nucl_seqs,
    genetic.code = bacterial_gen_code,
    if.fuzzy.codon = "solve"
  )

if (length(aa_seqs) == 0) {
  tibble(
    my_cds_id      = character(),
    ref_cf_full    = character(),
    ref_cf_base    = character(),
    score          = numeric(),
    max_score      = numeric(),
    pid            = numeric()
  ) %>%
    readr::write_tsv(
      file.path(
        config$global_alignment_scores_dir,
        paste0("global_alignment_scores_cds_around_cf_", genome_id, ".tsv.gz")
      )
    )

  return(NULL)
}

  # global alignments: each CDS amino acid vs every reference AA
  alignment_scores <- bind_rows(
  purrr::map(
    seq_along(aa_seqs),
    function(x) myalign2(CF_reference_aa_seqs, aa_seqs[x])
  )
) %>%
  transmute(
    my_cds_id      = seq2_name,              # CDS from genome
    ref_cf_full    = seq1_name,              # full DIAMOND subject id (may have _1, _2,...)
    ref_cf_base    = sub("_[0-9]+$", "", seq1_name),  # strip trailing _# for lookup
    score          = score,
    max_score      = max_score,
    pid            = pid
  )

  readr::write_tsv(
    alignment_scores,
    file.path(
      config$global_alignment_scores_dir,
      paste0("global_alignment_scores_cds_around_cf_", genome_id, ".tsv.gz")
    )
  )

  NULL
}

cat("Done!\n")
