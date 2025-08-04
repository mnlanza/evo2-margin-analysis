#!/bin/bash
set -euo pipefail

# Arguments: aid and optional position
aid_filter="${1:-}"     # First argument is aid (optional)
pos_filter="${2:-}"     # Second argument is specific position (optional)

if [ -n "$pos_filter" ] && [ -z "$aid_filter" ]; then
    echo "Error: Cannot specify position without specifying aid"
    echo "Usage: $0 [AID] [POSITION]"
    echo "Examples:"
    echo "  $0          # Run all aids and positions"
    echo "  $0 BAA      # Run all positions for BAA"
    echo "  $0 BAA 154  # Run only position 154 for BAA"
    exit 1
fi

# Create all required directories upfront
mkdir -p input output jobs figures

# Create aid-specific directories from TSV
while IFS=$'\t' read -r aid _ _ _ _ _ _ _ _ _ || [ -n "$aid" ]; do
    if [ "$aid" = "aid" ]; then continue; fi
    aid_lower="$(echo "$aid" | tr '[:upper:]' '[:lower:]')"
    mkdir -p "output/${aid_lower}" "figures/${aid_lower}"
done < "input/updated_data.tsv"

# Shared parameters (using no margin for now and then adding different margins in gen_for_cod_var_marg.py)
left_margin=20000
right_margin=0
input_fasta="input/human_contigs_src.fasta"

# Process mutations
while IFS=$'\t' read -r aid gene contig start end strand flipped src_codon tgt_codon mut_pos || [ -n "$contig" ]; do
    # Skip header line
    if [ "$aid" = "aid" ]; then continue; fi
    
    # Skip if aid filter is set and doesn't match
    if [ -n "$aid_filter" ] && [ "$aid" != "$aid_filter" ]; then continue; fi

    # Skip if position filter is set and doesn't match
    if [ -n "$pos_filter" ] && [ "$mut_pos" != "$pos_filter" ]; then continue; fi

    echo "Running $contig at position $mut_pos (gene range $start-$end, aid: $aid)"
    ./scripts/analyze_codon_pos_marg.sh "$mut_pos" "$start" "$end" "$contig" "$aid" "$input_fasta" "$tgt_codon" "$left_margin" "$right_margin"
done < "input/updated_data.tsv"