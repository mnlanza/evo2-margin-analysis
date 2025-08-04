#!/bin/bash

# Exit on error, undefined var, or pipeline fail
set -euo pipefail

# === Usage Check ===
if [ $# -lt 6 ]; then
  echo "Usage: $0 POS GENE_START GENE_END SEQ_ID AID INPUT_FASTA TGT_CODON [left_margin=2000] [right_margin=1000]"
  exit 1
fi

# === Positional and Optional Arguments ===
pos="$1"
gene_start="$2"
gene_end="$3"
seq_id="$4"
aid="$5"
input_fasta="$6"
tgt_codon="$7"
left_margin="${8:-2000}"
right_margin="${9:-1000}"

# Debug output
echo "Debug: Checking parameters..."
echo "pos: $pos"
echo "gene_start: $gene_start"
echo "gene_end: $gene_end"
echo "seq_id: $seq_id"
echo "aid: $aid"
echo "input_fasta: $input_fasta"
echo "tgt_codon: $tgt_codon"
echo "left_margin: $left_margin"
echo "right_margin: $right_margin"

# Check required files exist
if [ ! -f "$input_fasta" ]; then
  echo "Error: Input FASTA file $input_fasta not found"
  exit 1
fi

if [ ! -f "input/codon_table" ]; then
  echo "Error: input/codon_table not found"
  exit 1
fi

# Create required directories
mkdir -p output jobs figures

# Define job ID first (moved up)
# Convert to lowercase and replace underscores with hyphens
job_id="$(echo "$aid" | tr '[:upper:]' '[:lower:]')"
job_version="$(echo "${seq_id//_/-}-${pos}-margin-${left_margin}" | tr '[:upper:]' '[:lower:]' | tr '_' '-')"

# Validate job names match GCP requirements (lowercase letters, numbers, and hyphens only)
if [[ ! "$job_id" =~ ^[a-z][a-z0-9-]*[a-z0-9]$ ]]; then
    echo "Error: Invalid job_id format: $job_id"
    exit 1
fi

if [[ ! "$job_version" =~ ^[a-z0-9][a-z0-9-]*[a-z0-9]$ ]]; then
    echo "Error: Invalid job_version format: $job_version"
    exit 1
fi

# Create aid-specific directories
mkdir -p "output/${job_id}" "figures/${job_id}"

# File naming
output_dir="output/${job_id}"
fasta_out="${output_dir}/margins_${seq_id}_${pos}.fasta"
fasta_info_table="${output_dir}/margins_${seq_id}_${pos}.tab"
job_dir="jobs/${job_id}-${job_version}"
logits_dir="$(pwd)/${job_dir}/output"  # Add logits_dir for clarity

# Create job directory before running anything
mkdir -p "$job_dir"

# Remove old files if they exist
echo "Removing old files if they exist..."
rm -f "$fasta_info_table"  # Changed from plot_info_table to fasta_info_table
rm -f "$fasta_out"

# generate all codon variants at $pos with different left margins
# change seq-id for each population
python3 scripts/gen_for_cod_var_marg.py \
  --fasta "$input_fasta" \
  --codon-table input/codon_table \
  --aa-coord "$pos" \
  --seq-id "$seq_id" \
  --output-fasta "$fasta_out" \
  --plot_info_table "$fasta_info_table" \
  --left-margin "$left_margin" \
  --right-margin "$right_margin" \
  --gene-start "$gene_start" \
  --gene-end "$gene_end" \
  --aid "$aid" \
  --mut-codon "$tgt_codon"

# # submit job
# evo_gcp submit --job "$job_id" \
#   --input_fasta "$(pwd)/$fasta_out" \
#   --job_version "$job_version" \
#   --output_type logits \
#   --wait

# download results
evo_gcp download --job "$job_id" \
  --job_version "$job_version" \
  --jobs_dir "$(pwd)/jobs"

echo "Processing files:"
echo "  FASTA: $fasta_out"
echo "  Info: $fasta_info_table"
echo "  Logits: $logits_dir"

echo "Running format_marg_plot_data.py..."
python3 scripts/format_marg_plot_data.py \
  --layer_data_tsv "$fasta_info_table" \
  --fasta_margins "$fasta_out" \
  --output_dir "$output_dir" \
  --logits_dir "$logits_dir"

if [ $? -ne 0 ]; then
  echo "Error: format_marg_plot_data.py failed"
  exit 1
fi

echo "Data processing complete. Output files in: $output_dir"


# # Create strand comparison table
# Rscript -e "
# source('scripts/create_strand_table.r')
# create_strand_table(
#   ifn_tsv='${job_dir}/output/input_summary.txt',
#   ofn='$compare_out')
# "

# # Plot strand comparison
# Rscript -e "
# source('scripts/plot_strand_scatter.r')
# plot_strand_scatter(
#   ifn_tab='$compare_out',
#   ifn_codon="$plot_info_table",
#   title='${aid}_pos_${pos}',
#   fdir='figures/${job_id}')
# "
