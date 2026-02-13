# Snakefile (v0): link inputs + run Bakta

import os, sys, yaml, re
from snakemake.shell import shell

RUNCFG = "config/run.yaml"

with open(RUNCFG) as f:
    cfg = yaml.safe_load(f)

SAMPLES_TSV = cfg.get("samples_tsv")

if not SAMPLES_TSV:
    sys.exit("ERROR: samples_tsv is missing in config/run.yaml")

############################################
# parse samples.tsv (tabs OR spaces, 2 or 3 cols)
############################################
with open(SAMPLES_TSV) as fh:
    raw_lines = [l.strip() for l in fh if l.strip()]

# split header + rows
parts = [re.split(r"\s+", l) for l in raw_lines]
header = parts[0]
rows = parts[1:]

# normalize header names (handle Sample vs sample, etc.)
header_lc = [h.lower() for h in header]

if "sample" not in header_lc or "fasta" not in header_lc:
    sys.exit("ERROR: samples.tsv must have at least 'sample' and 'fasta' columns.")

has_gff = "gff3" in header_lc

sample_idx = header_lc.index("sample")
fasta_idx  = header_lc.index("fasta")
gff_idx    = header_lc.index("gff3") if has_gff else None

SAMPLES = []
FA = {} # sample -> fasta path
GFF = {} # sample -> gff path (if provided)

for row in rows:
    sid = row[sample_idx]
    fa_path = row[fasta_idx]

    SAMPLES.append(sid)
    FA[sid] = fa_path

    if has_gff:
        gff_path = row[gff_idx]
        GFF[sid] = gff_path

HAS_GFF = has_gff  # expose this for downstream logic
SAMPLE_LIST_STR = " ".join(SAMPLES)

############################################
# pull config values
############################################
GENOME_LINK_DIR = cfg.get("genome_link_dir")  # /work/genomes

if not GENOME_LINK_DIR:
    sys.exit("ERROR: genome_link_dir is missing in config/run.yaml")

# --- Bakta settings ---
BAKTA_PATH = cfg.get("bakta_path", "")
BAKTA_DB = cfg.get("bakta_db")
BAKTA_OUT = cfg.get("bakta_outdir", "results/bakta")
BAKTA_THREADS = int(cfg.get("bakta_threads", 8))
BAKTA_GENUS = cfg.get("bakta_genus", "Escherichia")
BAKTA_GRAM = cfg.get("bakta_gram", "-")
BAKTA_BIN = cfg.get("bakta_bin", "bakta")

# --- Tool binaries ---
RSCRIPT_BIN = cfg["tools"]["rscript"]
DIAMOND_BIN = cfg["tools"]["diamond"]

# --- CF classification settings ---
CF_AA_DB = cfg.get("cf_reference_aa_db")
CF_DMND_DB = cfg.get("cf_reference_dmnd")
CF_RULES = cfg.get("cf_rules_table")
CF_LOOKUP = cfg.get("cf_lookup_table")

CLASS_THREADS = int(cfg.get("classification_threads", 8))
CLASS_OUT = cfg.get("cf_classification_outdir", "results/cf_classification")
TMP_DIR = cfg.get("tmp_dir", "results/tmp")

# CDS feature output base
CDS_BASE = cfg.get("cds_features_outdir", "results/cds_features")

############################################
# expected per-sample outputs
############################################
# Expected CDS output for all samples
LINKED_FASTAS = [f"{GENOME_LINK_DIR}/{s}.fa" for s in SAMPLES]

# Bakta (or user-provided) GFF at standard location
BAKTA_GFFS    = [f"{BAKTA_OUT}/{s}/{s}.gff3" for s in SAMPLES]

# CDS extraction outputs per sample
CDS_COORDS = [f"{CDS_BASE}/cds_coords/cds_features_{s}.tsv.gz" for s in SAMPLES]
CDS_NUCL = [f"{CDS_BASE}/cds_nucl_seqs/cds_nucl_seqs_{s}.fa.gz" for s in SAMPLES]
CDS_AA = [f"{CDS_BASE}/cds_amino_seqs/cds_amino_seqs_{s}.fa.gz" for s in SAMPLES]
ALL_CDS_OUTPUTS = CDS_COORDS + CDS_NUCL + CDS_AA

# DIAMOND output per sample
DIAMOND_OUT = [
    f"{TMP_DIR}/diamond_aligned_amino_cf_db/diamond_align_{s}.tsv"
    for s in SAMPLES
]

# final classification table per sample
FINAL_CLASS_TSV = [
    f"{CLASS_OUT}/{s}/cf_classification_min_classify_score_fraction_0.95.tsv"
    for s in SAMPLES
]

RUN_SUMMARY_TSV = f"{CLASS_OUT}/ecgc_run_summary.tsv"

############################################
# rule all: what a full run produces
############################################
rule all:
    input:
        LINKED_FASTAS,
        BAKTA_GFFS,
        ALL_CDS_OUTPUTS,
        DIAMOND_OUT,
        FINAL_CLASS_TSV,
        RUN_SUMMARY_TSV


############################################
# 1) link_inputs:
# make per-sample symlink from raw fasta -> {GENOME_LINK_DIR}/{sample}.fa
############################################
rule link_inputs:
    input:
        src_fa = lambda wc: FA[wc.sample]
    output:
        fa = f"{GENOME_LINK_DIR}" + "/{sample}.fa"
    run:
        os.makedirs(GENOME_LINK_DIR, exist_ok=True)
        src_fa = os.path.abspath(input.src_fa)
        dst = output.fa
        if os.path.islink(dst) or os.path.exists(dst):
            os.remove(dst)
        os.symlink(src_fa, dst)
        print(f"[link_inputs] linked {wildcards.sample} -> {dst}")

############################################
# 2) prepare_gff:
# if samples.tsv provided gff3, just link it into {BAKTA_OUT}/{sample}/{sample}.gff3
# otherwise run Bakta to generate it
############################################
rule prepare_gff:
    input:
        fa = rules.link_inputs.output.fa
    output:
        gff = f"{BAKTA_OUT}" + "/{sample}/{sample}.gff3"
    threads: BAKTA_THREADS
    run:
        os.makedirs(f"{BAKTA_OUT}/{wildcards.sample}", exist_ok=True)

        if HAS_GFF:
            # User provided GFF in samples.tsv, link it
            src_gff_abs = os.path.abspath(GFF[wildcards.sample])
            dst_gff = output.gff

            if os.path.islink(dst_gff) or os.path.exists(dst_gff):
                os.remove(dst_gff)

            os.symlink(src_gff_abs, dst_gff)
        else:
            # Run Bakta to generate GFF3
            shell(r"""
                set -euo pipefail

                if [ -n "{BAKTA_PATH}" ]; then
                    export PATH="{BAKTA_PATH}:$PATH"
                fi

                "{BAKTA_BIN}" \
                  --db {BAKTA_DB} \
                  --genus '{BAKTA_GENUS}' \
                  --gram '{BAKTA_GRAM}' \
                  --output {BAKTA_OUT}/{wildcards.sample} \
                  --prefix {wildcards.sample} \
                  --threads {threads} \
                  --force \
                  {input.fa}
                # sanity check
                test -s {output.gff}
            """)
            

############################################
# 3) cds_extract_all:
# run one R script that:
#   - reads config/run.yaml
#   - reads ALL Bakta GFFs + genome FASTAs
#   - writes:
#       {CDS_BASE}/cds_coords/cds_features_<sample>.tsv.gz
#       {CDS_BASE}/cds_nucl_seqs/cds_nucl_seqs_<sample>.fa.gz
#       {CDS_BASE}/cds_amino_seqs/cds_amino_seqs_<sample>.fa.gz
############################################
rule cds_extract_all:
    input:
        cfg  = RUNCFG,
        gffs = BAKTA_GFFS,
        fas  = LINKED_FASTAS
    output:
        coords = CDS_COORDS,
        nucl   = CDS_NUCL,
        aa     = CDS_AA
    threads: 1
    shell:
        r"""
        set -euo pipefail
        {RSCRIPT_BIN} \
          code/cds_extraction/get_bakta_cds_defs_and_seqs.R \
          {input.cfg}
        """

############################################
# 4) diamond_search_all:
# blastp each sample's CDS amino acid FASTA vs CF ref db
############################################ 
rule diamond_search_all:
    input:
        aa_fas = CDS_AA,
        dmnd_db = CF_DMND_DB
    output:
        DIAMOND_OUT
    threads: CLASS_THREADS
    shell:
        r"""
        set -euo pipefail
        mkdir -p {TMP_DIR}/diamond_aligned_amino_cf_db

        # loop over samples and run diamond blastp for each
        for S in {SAMPLE_LIST_STR}; do
            QUERY="{CDS_BASE}/cds_amino_seqs/cds_amino_seqs_${{S}}.fa.gz"
            OUT="{TMP_DIR}/diamond_aligned_amino_cf_db/diamond_align_${{S}}.tsv"

            # only run if missing or empty (so reruns are fast)
            if [ ! -s "$OUT" ]; then
                {DIAMOND_BIN} blastp \
                  --query "$QUERY" \
                  --db {CF_DMND_DB} \
                  --out "$OUT" \
                  --outfmt 6 qseqid sseqid pident length slen qlen gaps evalue \
                  --max-target-seqs 200 \
                  --threads {threads}
            fi
        done

        # sanity check: make sure all expected outputs exist & non-empty
        for S in {SAMPLE_LIST_STR}; do
            test -f "{TMP_DIR}/diamond_aligned_amino_cf_db/diamond_align_${{S}}.tsv"
        done
        """

############################################
# 5) classify_cf:
# run R classification logic that:
#   - filters hits
#   - calls CFs per operon
#   - writes per-sample classification tables
############################################
rule classify_cf:
    input:
        cfg = RUNCFG,
        cds_out = ALL_CDS_OUTPUTS,
        diamond = DIAMOND_OUT,
        rules_tsv = CF_RULES,
        lookup_tsv = CF_LOOKUP
    output:
        FINAL_CLASS_TSV,
        RUN_SUMMARY_TSV
    threads: CLASS_THREADS
    shell:
        r"""
        set -euo pipefail
        mkdir -p {CLASS_OUT}

        # 1. Filter diamond hits down to CF-related matches
        {RSCRIPT_BIN} \
        code/cluster_classification/subscripts/filter_diamond_cds_cf_related.R \
        {RUNCFG}

        # 2. Assign CF / fimbrial cluster classifications based on rules
        {RSCRIPT_BIN} \
        code/cluster_classification/subscripts/classify_cf_profile.R \
        {RUNCFG}

        # 3. Sanity check: make sure all per-sample classification tables were created
        for S in {SAMPLE_LIST_STR}; do
            test -f "{CLASS_OUT}/${{S}}/cf_classification_min_classify_score_fraction_0.95.tsv"
        done

        # 4. Run-level summary across all genomes
        {RSCRIPT_BIN} \
        code/cluster_classification/subscripts/make_ecgc_run_summary.R \
        {RUNCFG}

        # sanity check: summary exists
        test -s "{CLASS_OUT}/ecgc_run_summary.tsv"
        """

    