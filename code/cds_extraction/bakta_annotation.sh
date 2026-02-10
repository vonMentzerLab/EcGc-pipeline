#!/usr/bin/bash

fa_dir="$1" # Contig fasta files directory.
bakta_db="$2" # Bakta database
no_of_cores="$3"


# Run bakta on all found genomes and put the resulting files in the same directories as the .fa files. ( {//} is the directory of a file )
find "${fa_dir}" -maxdepth 1 -name "*fa.gz" -print0 | parallel -0 -j "${no_of_cores}" --resume-failed --joblog "${fa_dir}/bakta_parallel.log" bakta --db "${bakta_db}" --threads 1 --genus 'Escherichia' --gram '-' --output {//} {}
