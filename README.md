# EcGc Pipeline

A Snakemake pipeline for detecting and classifying colonization factor (CF) gene clusters from bacterial genome assemblies.

The pipeline supports:

- Running **Bakta annotation automatically**
- Or using **pre-annotated GFF3 files**
- CDS extraction
- DIAMOND-based gene matching
- Operon-level CF classification

---

# Overview

For each genome, the pipeline:

1. Links the input FASTA file into a working directory
2. Runs Bakta (if no GFF3 file is provided)
3. Extracts CDS coordinates and sequences
4. Runs DIAMOND against a CF reference database
5. Classifies genes and operons based on predefined CF rules

The pipeline is fully configurable via `config/run.yaml`.

---

# Installation

Two conda environments are provided in the `env/` directory:

- `env/conda_bakta.yaml`
- `env/conda_cluster_classification.yaml`

Create them manually (recommended):

```bash
mamba env create -f env/conda_bakta.yaml -n bakta_env
mamba env create -f env/conda_cluster_classification.yaml -n ecgc
```

(`conda` may be used instead of `mamba`, but mamba is recommended.)

Activate the classification environment before running Snakemake:

```bash
conda activate ecgc
```

---

# Important: No --use-conda

This pipeline does **not** use Snakemake's `--use-conda`.

Instead, it expects required tools to be:

- Available in your PATH  
  **or**
- Explicitly defined in `run.yaml`

The following tools are configurable:

```yaml
bakta_bin: bakta
rscript_bin: Rscript
diamond_bin: diamond
```

If the tools are not in your PATH, replace them with full paths, e.g.:

```yaml
bakta_bin: /full/path/to/bakta
rscript_bin: /full/path/to/Rscript
diamond_bin: /full/path/to/diamond
```

---

# Input Files

## 1. samples.tsv

This file defines which genomes to analyse.

It must contain at least:

- `sample`
- `fasta`

Optionally:

- `gff3`

---

## Example 1 — Run Bakta Automatically

```
sample      fasta
GVE0064     /data/genomes/GVE0064.fa.gz
GVE0065     /data/genomes/GVE0065.fa.gz
```

In this case:

- Bakta runs automatically.
- GFF3 files are written to `bakta_outdir`.

---

## Example 2 — Use Existing GFF3 Files

```
sample      fasta                              gff3
GVE0064     /data/genomes/GVE0064.fa.gz        /data/gff/GVE0064.gff3
GVE0065     /data/genomes/GVE0065.fa.gz        /data/gff/GVE0065.gff3
```

In this case:

- Bakta is NOT run.
- The provided GFF3 files are used directly.

---

# Configuration (run.yaml)

Key parameters:

```yaml
samples_tsv: config/samples.tsv

# Working directory for symlinked FASTA files
genome_link_dir: results/work/genomes

# Output directories
bakta_outdir: results/bakta
cds_features_outdir: results/cds_features
cf_classification_outdir: results/cf_classification
tmp_dir: results/tmp

# Reference databases
cf_reference_nt_db: references/your_reference.fasta
cf_reference_aa_db: references/your_reference_translated.fasta
cf_reference_dmnd: references/your_reference_translated.dmnd
cf_rules_table: references/cf_rules.tsv
cf_lookup_table: references/cf_lookup.tsv

# Bakta
bakta_db: /path/to/bakta_db
bakta_bin: bakta

# Classification thresholds
min_cov: 80
min_pident: 80
max_dist: 10000
cf_classify_min_score_fraction: 0.95
```

All paths may be relative or absolute.

---

# Running the Pipeline

From the repository root:

```bash
snakemake -p --cores 8
```

Options explained:

- `--cores 8` → maximum number of cores Snakemake may use
- `-p` → print executed shell commands

---

# Output Structure

```
results/
├── bakta/
├── cds_features/
├── tmp/
└── cf_classification/
```

---

# Gene-Level Output Files

### best_matching_gene_any_cf.tsv

Best DIAMOND hit per gene (ignores CF consistency).

### classifications_by_gene_min_classify_score_fraction_*.tsv

Best gene assignment with CF-aware prioritization.

Important columns:

- `seqnames`
- `start`, `end`
- `strand`
- `operon_id`
- `my_cds_id`
- `gene`
- `cf`
- `score`
- `score_fraction`
- `pid`

---

# Operon-Level Output File

### cf_classification_min_classify_score_fraction_*.tsv

One row per detected operon.

Important columns:

- `cf_classification`
- `gene_profile`
- `mean_score_fraction`
- `operon_gene_count`
- `best_operon_cf_match`
- `best_operon_cf_match_all_genes_present`
- `best_operon_cf_match_genes_in_order`

---

# Reference Database Preparation

The DIAMOND database must be generated from the amino acid reference FASTA:

```bash
diamond makedb \
  --in references/baktamod_CFs_and_kappa_reference_db_translated.fasta \
  --db references/baktamod_CFs_and_kappa_reference_db_translated.dmnd
```

If using a different reference set, generate a new `.dmnd` file using the same command.

---

# Optional FASTA Output

If enabled:

```yaml
output_gene_seqs: true
output_operon_seqs: true
```

Sequences are written to:

```
results/cf_classification/<sample>/seqs/
├── genes/
└── operons/
```

---

# Notes

- CDS extraction assumes non-overlapping CDS features.
- Primarily tested with Bakta GFF3 output.
- Designed for large genome collections (1000+ genomes).
- Fully configurable via `run.yaml`.