#!/usr/bin/env bash
# Run the CF cluster classification pipeline
set -euo pipefail

CONFIG=${1:? "Usage: find_gene_clusters.sh <config.yaml>"}

# Resolve this script's directory (so relative paths work no matter where you run it)
SCRIPT_DIR="$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"

# --- Read needed config values (requires shyaml; it's in the conda env) ---
clustered_genes_reference_db_aa_seq=$(cat "$CONFIG" | shyaml get-value clustered_genes_reference_db_aa_seq)
cds_features_dir=$(cat "$CONFIG" | shyaml get-value cds_features_dir)
gene_search_dir=$(cat "$CONFIG" | shyaml get-value gene_search_dir)
cores=$(cat "$CONFIG" | shyaml get-value cores)

cds_amino_seqs_dir="${cds_features_dir}/cds_amino_seqs"
diamond_out_dir="${gene_search_dir}/diamond_aligned_amino_cf_db"

# Create output dirs
mkdir -p "$diamond_out_dir"

# --- 1) DIAMOND: query = CDS AA, subject = CF AA reference DB ---
"${SCRIPT_DIR}/subscripts/parallel_diamond_cds_amino_cf_ref.sh" \
  "$clustered_genes_reference_db_aa_seq" \
  "$cds_amino_seqs_dir" \
  "$diamond_out_dir" \
  "$cores"

# --- 2) Filter DIAMOND hits ---
Rscript "${SCRIPT_DIR}/subscripts/filter_diamond_cds_cf_related.R" "$CONFIG"

# --- 3) Classify operons / CF profile ---
Rscript "${SCRIPT_DIR}/subscripts/classify_cf_profile.R" "$CONFIG"

echo "Done."