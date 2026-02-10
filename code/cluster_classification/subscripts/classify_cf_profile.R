library(tidyverse)
library(cowplot)
library(doParallel)
library(foreach)
library(yaml)

read_cf_rules <- function(cf_rules_path){
  read_tsv(cf_rules_path, col_types = "cc") %>% 
    group_by(cf) %>% 
    mutate(cf_gene_order = row_number()) %>% 
    ungroup() %>% 
    mutate(gene = ifelse(str_detect(gene, ":"), gene, paste0(gene, ":", gene))) %>% 
    separate(gene, sep = ":", into = c("equiv_gene", "gene")) %>% 
    mutate(gene = str_split(gene, ",")) %>% 
    unnest(gene)
}

get_cf_required_gene_counts <- function(cf_rules){
  cf_rules %>%
    filter(required) %>% 
    select(cf, equiv_gene) %>% 
    distinct() %>%
    count(cf, name = "req_gene_count_in_cf_rules")
}

find_prob_cf_coords <- function(cf_coords){
  operon_zone <- cf_coords %>%
    filter(diamond_cf_candidate) %>% 
    mutate(first_cf_in_operon = (start - lag(end)) > config$max_dist) %>% 
    replace_na(list(first_cf_in_operon = TRUE)) %>% 
    mutate(operon_id = paste0(seqnames, ".", cumsum(first_cf_in_operon))) %>% 
    group_by(operon_id, seqnames) %>% 
    summarise(
      operon_zone_start = min(start),
      operon_zone_end   = max(end),
      .groups = "drop"
    ) %>% 
    rowwise() %>% 
    mutate(
      operon_zone_start = max(operon_zone_start - config$max_dist, 1),
      operon_zone_end   = operon_zone_end + config$max_dist
    ) %>% 
    ungroup()

  # include nearby genes in that zone
  operon_cf_coords <- bind_rows(
    purrr::map(
      unique(operon_zone$operon_id),
      function(x)
        inner_join(cf_coords, filter(operon_zone, operon_id == x), by = "seqnames") %>%
          filter(start >= operon_zone_start, end <= operon_zone_end)
    )
  )

  prob_cf_coords <- operon_cf_coords %>% 
    group_by(seqnames, operon_id, diamond_cf_candidate) %>% 
    mutate(
      first_cf_gene_start = ifelse(diamond_cf_candidate & start == min(start), start, 0),
      last_cf_gene_end    = ifelse(diamond_cf_candidate & end   == max(end),   end,   0)
    )  %>% 
    group_by(seqnames, operon_id) %>% 
    mutate(
      first_cf_gene_start = sum(first_cf_gene_start),
      last_cf_gene_end    = sum(last_cf_gene_end)
    ) %>% 
    ungroup() %>% 
    mutate(prob_cf_operon_gene = start >= first_cf_gene_start & end <= last_cf_gene_end) %>% 
    filter(prob_cf_operon_gene)

  operon_strands <- prob_cf_coords %>% 
    count(operon_id, strand) %>% 
    arrange(operon_id, strand) %>% 
    group_by(operon_id) %>% 
    slice_max(order_by = n, with_ties = FALSE) %>% 
    ungroup() %>% 
    select(operon_id, operon_strand = strand)

  prob_cf_coords <- prob_cf_coords %>% 
    left_join(operon_strands, by = "operon_id") %>% 
    group_by(operon_id) %>% 
    mutate(gene_order = ifelse(operon_strand == "-", dplyr::n() - row_number() + 1, row_number()))  %>% 
    ungroup()

  prob_cf_coords
}

make_cf_score_df <- function(prob_cf_coords){

  prob_cf_coords %>%
    # bring in alignment info per CDS
    left_join(cf_align_scores, by = "my_cds_id") %>%
    # per-hit metrics
    mutate(score_fraction = score / max_score) %>%
    # map reference sequence name -> gene short_name
    # cf_names: short_name, name
    left_join(
      cf_names %>% 
        dplyr::rename(
          ref_cf_base = name,
          gene        = short_name
        ),
      by = "ref_cf_base"
    ) %>%
    # map gene -> CF & required/etc
    left_join(cf_rules, by = "gene")
}

classify_cfs <- function(cf_score_df, cf_classify_min_score_fraction){

  genome_id <- unique(cf_score_df$genome_id)

  df_max_per_cf <- cf_score_df %>% 
    group_by(operon_id, cf, gene_order) %>%
    slice_max(order_by = score_fraction, with_ties = FALSE) %>% 
    ungroup()

  operon_gene_count <- cf_score_df %>%
    group_by(operon_id) %>%
    slice_max(gene_order, with_ties = FALSE) %>%
    ungroup() %>%
    select(operon_id, operon_gene_count = gene_order)

  df_all_genes_present <- df_max_per_cf  %>%
    filter(
      score_fraction >= cf_classify_min_score_fraction,
      required
    ) %>%
    select(operon_id, cf, equiv_gene) %>% 
    distinct() %>% 
    group_by(operon_id, cf) %>%
    arrange(operon_id, cf, equiv_gene) %>% 
    summarise(cf_present_genes = paste0(equiv_gene, collapse = " + "), .groups = "drop") %>% 
    left_join(
      cf_rules %>%
        filter(required) %>% 
        select(cf, equiv_gene) %>% 
        distinct() %>% 
        group_by(cf) %>%
        arrange(cf, equiv_gene) %>% 
        summarise(cf_required_genes = paste0(equiv_gene, collapse = " + ")),
      by = "cf"
    ) %>% 
    mutate(all_genes_present = cf_present_genes == cf_required_genes) %>% 
    select(operon_id, cf, all_genes_present)

  df_cf_gene_count <- df_max_per_cf %>% 
    filter(score_fraction >= cf_classify_min_score_fraction) %>%
    count(operon_id, cf, name = "gene_count") %>% 
    bind_rows(
      df_max_per_cf %>%
        select(operon_id, cf) %>%
        mutate(gene_count = 0)
    ) %>% 
    group_by(operon_id, cf) %>% 
    summarise(gene_count = sum(gene_count), .groups = "drop") %>% 
    left_join(cf_req_gene_counts, by = "cf")

  df_cf_in_order <- df_max_per_cf %>%
    group_by(operon_id, cf) %>% 
    filter(required) %>% 
    mutate(equiv_gene = ifelse(score_fraction >= cf_classify_min_score_fraction, equiv_gene, "too_low_score")) %>% 
    summarise(
      obs_gene_order = paste0(equiv_gene, collapse = "_"),
      .groups = "drop"
    ) %>%
    mutate(
      obs_gene_order_trimmed = str_replace(obs_gene_order, "^(too_low_score_)+", ""),
      obs_gene_order_trimmed = str_replace(obs_gene_order_trimmed, "(_too_low_score)+$", "")
    ) %>% 
    left_join(
      cf_rules %>%
        filter(required) %>%
        select(cf, equiv_gene, cf_gene_order) %>%
        distinct() %>%
        group_by(cf) %>%
        arrange(cf, cf_gene_order) %>%
        summarise(req_gene_order = paste0(equiv_gene, collapse = "_")),
      by = "cf"
    ) %>% 
    mutate(genes_in_order = str_detect(req_gene_order, obs_gene_order_trimmed)) %>%
    select(operon_id, cf, genes_in_order)

  df_score_fraction_summary <- df_max_per_cf %>% 
    select(operon_id, cf, gene, score_fraction) %>% 
    group_by(operon_id, cf) %>% 
    summarise(mean_score_fraction = mean(score_fraction), .groups = "drop")

  df_summary <- df_cf_gene_count %>% 
    left_join(df_cf_in_order, by = c("operon_id", "cf")) %>% 
    left_join(df_score_fraction_summary, by = c("operon_id", "cf")) %>% 
    left_join(df_all_genes_present, by = c("operon_id", "cf")) %>% 
    replace_na(list(genes_in_order = FALSE, all_genes_present = FALSE))

  complete_cfs <- df_summary %>% 
    filter(all_genes_present, genes_in_order) %>% 
    group_by(operon_id) %>% 
    slice_max(order_by = mean_score_fraction, with_ties = FALSE) %>% 
    ungroup()

  complete_cf_best_hits <- df_max_per_cf %>% 
    inner_join(
      complete_cfs %>% select(operon_id, cf),
      by = c("operon_id", "cf")
    ) %>% 
    select(
      seqnames, start, end, strand,
      my_cds_id, operon_id, operon_strand,
      genome_id, cf, gene, score_fraction, gene_order
    ) %>% 
    mutate(classification_type = "complete")

  incomplete_cfs <- df_summary %>% 
    filter(!(operon_id %in% complete_cfs$operon_id)) %>% 
    group_by(operon_id) %>% 
    slice_max(order_by = mean_score_fraction, with_ties = FALSE) %>% 
    ungroup()

  initial_cf_classification <- complete_cfs %>% 
    bind_rows(incomplete_cfs) %>% 
    mutate(
      all_genes_and_ordered = all_genes_present & genes_in_order,
      genome_id = genome_id
    )

  no_match_cfs <- initial_cf_classification %>%
    anti_join(
      df_max_per_cf %>%
        filter(score_fraction >= cf_classify_min_score_fraction) %>% 
        select(operon_id) %>% 
        distinct(),
      by = "operon_id"
    )

  if (nrow(initial_cf_classification) > 0) {

    hybrid_cf_tmp <- initial_cf_classification %>% 
      rename(best_operon_cf_match = cf) %>% 
      left_join(df_max_per_cf, by = c("operon_id", "genome_id")) %>% 
      rowwise() %>% 
      mutate(cf_matching_gene = (
        score_fraction >= !!cf_classify_min_score_fraction &&
        best_operon_cf_match == cf
      )) %>% 
      ungroup()

    hybrid_cf_best_cf_matching_genes <- hybrid_cf_tmp %>%
      filter(cf_matching_gene)

    hybrid_cf_best_cf_no_pass_matching <- hybrid_cf_tmp %>% 
      anti_join(
        hybrid_cf_best_cf_matching_genes %>% 
          select(operon_id, gene_order),
        by = c("operon_id", "gene_order")
      ) %>% 
      group_by(operon_id, gene_order) %>% 
      slice_max(order_by = score_fraction, with_ties = FALSE)

    hybrid_cf_genes <- hybrid_cf_best_cf_matching_genes %>% 
      bind_rows(hybrid_cf_best_cf_no_pass_matching) %>% 
      arrange(operon_id, gene_order) %>% 
      group_by(operon_id) %>% 
      mutate(hybrid_min_score_qualifying = all(score_fraction >= cf_classify_min_score_fraction)) %>% 
      ungroup()

    hybrid_cfs <- bind_rows(
      hybrid_cf_genes %>%
        filter(hybrid_min_score_qualifying) %>%
        group_by(operon_id, best_operon_cf_match) %>% 
        summarise(
          cf_classification = paste0(unique(cf), collapse = " / "),
          gene_profile = paste0(gene, collapse = " + "),
          .groups = "drop"
        ),
      hybrid_cf_genes %>% 
        filter(!hybrid_min_score_qualifying) %>% 
        group_by(operon_id, best_operon_cf_match) %>% 
        filter(score_fraction >= cf_classify_min_score_fraction) %>% 
        summarise(
          cf_classification = paste0(c(unique(cf), "Unknown"), collapse = " / "),
          .groups = "drop"
        ) %>% 
        bind_rows(
          no_match_cfs %>%
            select(operon_id, best_operon_cf_match = cf) %>% 
            mutate(cf_classification = "Unknown")
        ) %>% 
        left_join(
          hybrid_cf_genes %>% 
            filter(!hybrid_min_score_qualifying) %>% 
            mutate(gene = ifelse(score_fraction >= cf_classify_min_score_fraction, gene, "Unknown")) %>% 
            group_by(operon_id) %>% 
            summarise(gene_profile = paste0(gene, collapse = " + "), .groups = "drop"),
          by = "operon_id"
        )
    ) %>% mutate(genome_id = genome_id)

    cf_classification <- hybrid_cfs %>% 
      left_join(
        hybrid_cf_genes %>% 
          group_by(operon_id) %>% 
          summarise(mean_score_fraction = mean(score_fraction), .groups = "drop"),
        by = c("operon_id")
      ) %>% 
      left_join(
        initial_cf_classification %>%
          rename(
            best_operon_cf_match = cf,
            best_operon_cf_match_gene_count = gene_count,
            best_operon_cf_match_all_genes_present = all_genes_present,
            best_operon_cf_match_genes_in_order = genes_in_order,
            best_operon_cf_match_mean_score_fraction = mean_score_fraction,
            best_operon_cf_match_all_genes_and_ordered = all_genes_and_ordered
          ),
        by = c("operon_id", "genome_id", "best_operon_cf_match")
      ) %>% 
      left_join(operon_gene_count, by = "operon_id") %>% 
      relocate(
        genome_id,
        operon_id,
        cf_classification,
        gene_profile,
        mean_score_fraction,
        best_operon_cf_match,
        best_operon_cf_match_mean_score_fraction,
        operon_gene_count,
        best_operon_cf_match_gene_count,
        req_gene_count_in_cf_rules,
        best_operon_cf_match_all_genes_present,
        best_operon_cf_match_genes_in_order,
        best_operon_cf_match_all_genes_and_ordered
      )

    classifications_by_gene <- hybrid_cf_genes %>% 
      select(
        seqnames, start, end, strand,
        genome_id, operon_id, operon_strand, my_cds_id,
        best_operon_cf_match, cf, gene, score_fraction,
        score, max_score, pid, gene_order
      )

    operon_coords <- classifications_by_gene %>% 
      select(-strand) %>% 
      group_by(genome_id, operon_id, seqnames, strand = operon_strand) %>% 
      summarise(
        start = min(start),
        end   = max(end),
        .groups = "drop"
      )

    cf_classification <- cf_classification %>% 
      left_join(operon_coords, by = c("genome_id", "operon_id"))

  } else {
    cf_classification       <- tibble()
    classifications_by_gene <- tibble()
  }

  best_matching_gene_any_cf <- cf_score_df %>%
    group_by(operon_id, gene_order) %>%
    slice_max(score_fraction, with_ties = FALSE) %>% 
    ungroup() %>% 
    select(
      seqnames, start, end, strand,
      genome_id, operon_id, operon_strand, my_cds_id,
      cf, gene, score_fraction,
      score, max_score, pid, gene_order
    )

  list(
    cf_classification = cf_classification,
    classifications_by_gene = classifications_by_gene,
    best_matching_gene_any_cf = best_matching_gene_any_cf
  )
}

##############################
# Config
args <- commandArgs(trailingOnly = TRUE)
config <- read_yaml(args[1])

# remap new config keys -> older expectations
config$gene_search_dir    <- config$tmp_dir
config$out_base           <- config$cf_classification_outdir
config$cores              <- config$classification_threads

config$cf_filtered_cds_coords_dir  <- file.path(config$gene_search_dir, "cf_cds_coords")
config$global_alignment_scores_dir <- file.path(config$gene_search_dir, "cds_global_alignment_scores")

cf_names <- read_tsv(config$cf_lookup_table)
cf_rules <- read_cf_rules(config$cf_rules_table)
cf_req_gene_counts <- get_cf_required_gene_counts(cf_rules)

files <- list.files(config$cf_filtered_cds_coords_dir)

comb <- function(x, y){
  list(
    cf_classification         = bind_rows(x[[1]], y[[1]]),
    classifications_by_gene   = bind_rows(x[[2]], y[[2]]),
    best_matching_gene_any_cf = bind_rows(x[[3]], y[[3]])
  )
}

cores  <- config$cores
cl <- makeCluster(cores)
registerDoParallel(cl)

results <- foreach(
  i = seq_along(files),
  .combine = "comb",
  .init = list(tibble(), tibble(), tibble()),
  .inorder = FALSE,
  .packages = c("tidyverse")
) %dopar% {

  file <- files[i]

  cf_coords <- read_tsv(file.path(config$cf_filtered_cds_coords_dir, file))
  genome_id <- unique(cf_coords$genome_id)

  cf_align_scores <- read_tsv(
    file.path(
      config$global_alignment_scores_dir,
      paste0("global_alignment_scores_cds_around_cf_", genome_id, ".tsv.gz")
    )
  )

  prob_cf_coords <- find_prob_cf_coords(cf_coords)

  df <- make_cf_score_df(prob_cf_coords)

  classify_cfs(df, config$cf_classify_min_score_fraction)
}

# now write outputs PER GENOME inside cf_classification_outdir/<genome_id>
# we'll split the big result list by genome and write to subdirs
dir.create(config$out_base, recursive = TRUE, showWarnings = FALSE)

# derive genome IDs from filenames like "cds_around_cf_<genome>.tsv.gz"
all_genomes <- files %>%
  sub("^cds_around_cf_", "", .) %>%
  sub("\\.tsv\\.gz$", "", .) %>%
  unique()

for (g in all_genomes) {
  out_dir_g <- file.path(config$out_base, g)
  dir.create(out_dir_g, recursive = TRUE, showWarnings = FALSE)

  write_tsv(
    results$cf_classification %>% filter(genome_id == g),
    file.path(out_dir_g, paste0("cf_classification_min_classify_score_fraction_", config$cf_classify_min_score_fraction, ".tsv"))
  )

  write_tsv(
    results$classifications_by_gene %>% filter(genome_id == g),
    file.path(out_dir_g, paste0("classifications_by_gene_min_classify_score_fraction_", config$cf_classify_min_score_fraction, ".tsv"))
  )

  write_tsv(
    results$best_matching_gene_any_cf %>% filter(genome_id == g),
    file.path(out_dir_g, "best_matching_gene_any_cf.tsv")
  )
}

# Optional seq output (unchanged logic except using out_base)
if (isTRUE(config$output_gene_seqs) || isTRUE(config$output_operon_seqs)) {
  library(plyranges)
  library(BSgenome)
  library(Rsamtools)
  library(Biostrings)

  load_fa_seqs <- function(fa_file){
    Biostrings::getSeq(Rsamtools::FaFile(fa_file))
  }

  make_bakta_contignames <- function(dnastringset){
    names(dnastringset) <- paste0("contig_", seq_along(dnastringset))
    dnastringset
  }

  # You'll need these in YAML if you want seq export:
  # genome_seqs_dir: <dir with original contig FASTAs>
  # fa_suffix: ".fa.gz" (or similar)
  if (isTRUE(config$output_gene_seqs)) {
    for (g in all_genomes) {
      gene_seq_out_dir <- file.path(config$out_base, g, "seqs", "genes")
      dir.create(gene_seq_out_dir, recursive = TRUE, showWarnings = FALSE)

      gr_tmp <- results$classifications_by_gene %>%
        filter(genome_id == g) %>%
        as_granges() # requires plyranges

      contig_seqs <- load_fa_seqs(
        file.path(config$genome_seqs_dir, paste0(g, config$fa_suffix))
      ) %>%
        make_bakta_contignames()

      nucl_seqs <- BSgenome::getSeq(contig_seqs, gr_tmp)
      names(nucl_seqs) <- paste0(
        gr_tmp$operon_id, ",", gr_tmp$my_cds_id, ",",
        gr_tmp$best_operon_cf_match, ",", gr_tmp$gene
      )

      writeXStringSet(
        nucl_seqs,
        file.path(gene_seq_out_dir, paste0(g, ".fa.gz")),
        compress = TRUE
      )
    }
  }

  if (isTRUE(config$output_operon_seqs)) {
    for (g in all_genomes) {
      operon_seq_out_dir <- file.path(config$out_base, g, "seqs", "operons")
      dir.create(operon_seq_out_dir, recursive = TRUE, showWarnings = FALSE)

      operons.df <- results$classifications_by_gene %>%
        filter(genome_id == g) %>%
        group_by(seqnames, strand = operon_strand, genome_id, operon_id, best_operon_cf_match) %>% 
        summarise(start = min(start), end = max(end), .groups = "drop")

      gr_tmp <- operons.df %>%
        filter(genome_id == g) %>%
        as_granges()

      contig_seqs <- load_fa_seqs(
        file.path(config$genome_seqs_dir, paste0(g, config$fa_suffix))
      ) %>%
        make_bakta_contignames()

      nucl_seqs <- BSgenome::getSeq(contig_seqs, gr_tmp)
      names(nucl_seqs) <- paste0(gr_tmp$operon_id, ",", gr_tmp$best_operon_cf_match)

      writeXStringSet(
        nucl_seqs,
        file.path(operon_seq_out_dir, paste0(g, ".fa.gz")),
        compress = TRUE
      )
    }
  }
}