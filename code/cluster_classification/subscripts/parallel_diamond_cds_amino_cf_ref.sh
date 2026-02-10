#!/usr/bin/bash

# arg1: fasta file (e.g. CF_reference_db_translated.fasta)
# arg2: cds amino seqs dir
# arg3: diamond output dir
# arg4: number of cores

mkdir -p "$3"

db_fa_path="$1"
db_fa_path_wo_ext="${db_fa_path%%.*}"
db_fa_file_wo_ext=$(basename -- "$db_fa_path_wo_ext")

diamond makedb --in "$db_fa_path" --db "$db_fa_path_wo_ext"

parallel -j "$4" diamond blastp --outfmt "6 qseqid sseqid pident length slen qlen gaps evalue" -d "$db_fa_path_wo_ext".dmnd -q {} -p 1 -o "$3"/diamond_align_"$db_fa_file_wo_ext"_{/.}.tsv ::: "$2"/cds_amino_seqs*.fa*
