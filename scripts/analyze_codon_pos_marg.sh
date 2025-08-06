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

# Create aid-specific directories
mkdir -p "output/${job_id}" "figures/${job_id}"

# File naming
output_dir="output/${job_id}"
fasta_out="${output_dir}/margins_${seq_id}_${pos}.fasta"
seq_info_table="${output_dir}/margins_${seq_id}_${pos}.tab"
job_dir="jobs/${job_id}-${job_version}"
logits_dir="$(pwd)/${job_dir}/output"  # Add logits_dir for clarity

# Create job directory before running anything
mkdir -p "$job_dir"

# Define query table path (same as seq_info_table but with _query suffix)
query_table="${seq_info_table%.tab}_query.tab"

# Remove old files if they exist
echo "Removing old files if they exist..."
rm -f "$seq_info_table"  # Changed from plot_info_table to seq_info_table
rm -f "$fasta_out"
rm -f "$query_table"

# generate all codon variants at $pos with different left margins
# change seq-id for each population
python3 scripts/gen_for_cod_var_marg.py \
  --fasta "$input_fasta" \
  --codon-table input/codon_table \
  --aa-coord "$pos" \
  --seq-id "$seq_id" \
  --output-fasta "$fasta_out" \
  --plot_info_table "$seq_info_table" \
  --query-table "$query_table" \
  --left-margin "$left_margin" \
  --right-margin "$right_margin" \
  --gene-start "$gene_start" \
  --gene-end "$gene_end" \
  --aid "$aid" \
  --mut-codon "$tgt_codon"

# # Submit job and wait for completion
# echo "Submitting job to GCP..."
# evo_gcp submit --job "$job_id" \
#   --input_fasta "$(pwd)/$fasta_out" \
#   --query_table "$(pwd)/$query_table" \
#   --job_version "$job_version" \
#   --output_type logits \
#   --wait

# # Download results
# echo "Downloading results from GCP..."
# evo_gcp download --job "$job_id" \
#   --job_version "$job_version" \
#   --jobs_dir "$(pwd)/jobs"

Rscript scripts/plot_margins.R \
  --seq_data "$seq_info_table" \
  --fasta_file "$fasta_out" \
  --aid "$aid" \
  --job_dir "$logits_dir"  # Use logits_dir instead of job_dir

if [ $? -ne 0 ]; then
  echo "Error: plot_margins.R failed"
  exit 1
fi
