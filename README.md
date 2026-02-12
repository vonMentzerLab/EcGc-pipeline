# EcGc: Colonization Factor Detection Pipeline

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