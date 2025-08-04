from utils import initialize_dict
import numpy as np
import pandas as pd
import os
import sys
from Bio import SeqIO

def format_marg_plot_data(layer_data_tsv, fasta_margins, output_dir, logits_dir):
    print("\nDEBUG: Starting format_marg_plot_data with:")
    print(f"layer_data_tsv: {layer_data_tsv}")
    print(f"fasta_margins: {fasta_margins}")
    print(f"output_dir: {output_dir}")
    print(f"logits_dir: {logits_dir}")
    
    # Make output_dir absolute if it isn't already
    if not os.path.isabs(output_dir):
        output_dir = os.path.abspath(output_dir)
    
    if not os.path.exists(logits_dir):
        print(f"Error: Logits directory not found: {logits_dir}")
        sys.exit(1)
    
    if not os.path.exists(fasta_margins):
        print(f"Error: FASTA file not found: {fasta_margins}")
        sys.exit(1)
    
    # Create output directory if it doesn't exist
    os.makedirs(output_dir, exist_ok=True)

    # Read TSV file (header and one line)
    with open(layer_data_tsv, 'r') as f:
        header = f.readline().strip().split('\t')
        values = f.readline().strip().split('\t')
        data = dict(zip(header, values))

    # Convert numeric values
    aa_coord = int(data['coord'])
    left_margin = int(data['left_margin'])
    start = int(data['start'])
    end = int(data['end'])

    # String values
    codon = data['codon']
    amino_acid = data['aa']
    aid = data['aid']
    mut_codon = data['mut_codon']
    mut_aa = data['mut_aa']
    syn_codon = data['synon_cod']

    # Debug: print values
    print("\nValues from TSV:")
    print(f"aa_coord: {aa_coord}")
    print(f"amino_acid: {amino_acid}")
    print(f"codon: {codon}")
    print(f"left_margin: {left_margin}")
    print(f"mut_codon: {mut_codon}")
    print(f"mut_aa: {mut_aa}")

    # Debug: print full paths being attempted
    print("\nAttempting to load files:")

    # Define the variants we want to analyze
    variants = {
        'src': {
            'aa': amino_acid,
            'codon': codon,
            'logits': os.path.join(logits_dir, f"input_{aa_coord}_{amino_acid}_{codon}_{left_margin}_logits.npy")
        },
        'stop': {
            'aa': 'Z',
            'codon': 'TAG',
            'logits': os.path.join(logits_dir, f"input_{aa_coord}_Z_TAG_{left_margin}_logits.npy")
        },
        'mut': {
            'aa': mut_aa,
            'codon': mut_codon,
            'logits': os.path.join(logits_dir, f"input_{aa_coord}_{mut_aa}_{mut_codon}_{left_margin}_logits.npy")
        },
        'syn': {
            'aa': amino_acid,
            'codon': syn_codon,
            'logits': os.path.join(logits_dir, f"input_{aa_coord}_{amino_acid}_{syn_codon}_{left_margin}_logits.npy")
        }
    }

    # Read sequences from the combined FASTA file
    print(f"\nReading sequences from: {fasta_margins}")
    all_variants = {}
    
    if not os.path.exists(fasta_margins):
        print(f"Error: Combined FASTA file not found: {fasta_margins}")
        sys.exit(1)
        
    for record in SeqIO.parse(fasta_margins, "fasta"):
        # Store all sequences - we'll pick the ones we need later
        all_variants[record.id] = str(record.seq)
    
    # Group sequences by variant type and margin size      
    variant_sequences = {}
    for header, sequence in all_variants.items():
        # Parse header parts (format: aa_coord_aa_codon_marginseq)
        try:
            # Split by underscore, but only first 3 times to keep margin sequence intact
            parts = header.split('_', 3)
            if len(parts) != 4:
                print(f"Warning: Unexpected header format: {header}")
                continue
                
            seq_coord, seq_aa, seq_codon, margin_seq = parts
            
            # Validate coord matches
            if int(seq_coord) != aa_coord:
                print(f"Warning: Coordinate mismatch in {header}")
                continue
            
            # Create variant key
            variant_key = f"{seq_coord}_{seq_aa}_{seq_codon}"
            
            # Calculate margin size from margin sequence
            margin_size = len(margin_seq)  # This is the actual margin sequence from the header
            
            # Store sequence with its margin size
            if variant_key not in variant_sequences:
                variant_sequences[variant_key] = {}
            variant_sequences[variant_key][margin_size] = sequence
            
            print(f"Processed variant {variant_key} with margin size {margin_size}")
            
        except ValueError as e:
            print(f"Error processing header {header}: {e}")
            continue

    # Initialize results dictionary
    results = {
        'variants': {},
        'metadata': {
            'aa_coord': aa_coord,
            'left_margin': left_margin,
            'start': start,
            'end': end,
            'aid': aid
        }
    }

    # Initialize data storage for both metrics
    entropy_data = []
    loglik_data = []
    gene_length = end - start + 1

    # Process each variant
    for var_type, var_info in variants.items():
        base_key = f"{aa_coord}_{var_info['aa']}_{var_info['codon']}"
        if base_key not in variant_sequences:
            print(f"Warning: No sequences found for {var_type} ({base_key})")
            continue

        print(f"\nProcessing {var_type} variant...")

        # Process each margin size
        for margin_size, sequence in sorted(variant_sequences[base_key].items()):
            seq_dict = {base_key: sequence}
            result_dict = initialize_dict(var_info['logits'], base_key, seq_dict)
            
            print(f"Processing margin size {margin_size} for {var_type}")
            
            # Get values only for gene region
            # For sequences with margins, we need to slice off the margin
            # For no-margin sequences (margin_size == 0), we use the whole sequence
            if margin_size == 0:
                gene_entropy = result_dict['entropy']
                gene_loglik = result_dict['log_likelihood']
            else:
                gene_entropy = result_dict['entropy'][margin_size:]
                gene_loglik = result_dict['log_likelihood'][margin_size:]
            
            # Debug prints
            print(f"\nSequence lengths for {var_type} with margin {margin_size}:")
            print(f"Full sequence length: {len(sequence)}")
            print(f"Gene entropy length: {len(gene_entropy)}")
            print(f"Gene loglik length: {len(gene_loglik)}")
            print(f"Expected gene length: {gene_length}")
            
            # Store data
            entropy_data.append({
                'margin_size': margin_size,
                'variant': var_type,
                'values': gene_entropy
            })
            loglik_data.append({
                'margin_size': margin_size,
                'variant': var_type,
                'values': gene_loglik
            })

    # Create output files
    print("\nDEBUG: Creating output files...")
    print(f"Output directory: {output_dir}")
    
    # 1. Entropy TSV
    entropy_file = os.path.join(output_dir, "entropy.tsv")
    print(f"Writing entropy file to: {entropy_file}")
    with open(entropy_file, 'w') as f:
        # Write metadata header as comments
        f.write(f"# Analysis ID: {aid}\n")
        f.write(f"# Gene coordinates: {start}-{end}\n")
        f.write(f"# Amino acid coordinate: {aa_coord}\n")
        f.write(f"# Gene length: {gene_length} bp\n")
        # Write column headers
        f.write("margin_size\tvariant\t" + "\t".join(f"pos_{i+1}" for i in range(gene_length)) + "\n")
        # Write data
        for entry in entropy_data:
            f.write(f"{entry['margin_size']}\t{entry['variant']}\t")
            f.write("\t".join(f"{v:.6f}" for v in entry['values']))
            f.write("\n")

    # 2. Log-likelihood TSV
    loglik_file = os.path.join(output_dir, "loglik.tsv")
    with open(loglik_file, 'w') as f:
        # Write metadata header as comments
        f.write(f"# Analysis ID: {aid}\n")
        f.write(f"# Gene coordinates: {start}-{end}\n")
        f.write(f"# Amino acid coordinate: {aa_coord}\n")
        f.write(f"# Gene length: {gene_length} bp\n")
        # Write column headers
        f.write("margin_size\tvariant\t" + "\t".join(f"pos_{i+1}" for i in range(gene_length)) + "\n")
        # Write data
        for entry in loglik_data:
            f.write(f"{entry['margin_size']}\t{entry['variant']}\t")
            f.write("\t".join(f"{v:.6f}" for v in entry['values']))
            f.write("\n")

    print(f"\nOutput files created:")
    print(f"Entropy data: {entropy_file}")
    print(f"Log-likelihood data: {loglik_file}")

    # Return metadata for reference
    return {
        'metadata': {
            'aa_coord': aa_coord,
            'left_margin': left_margin,
            'start': start,
            'end': end,
            'aid': aid,
            'files': {
                'entropy': entropy_file,
                'loglik': loglik_file
            }
        }
    }

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('--layer_data_tsv', required=True, help='Path to layer data TSV file')
    parser.add_argument('--fasta_margins', required=True, help='Path to FASTA margins file')
    parser.add_argument('--output_dir', required=True, help='Output directory')
    parser.add_argument('--logits_dir', required=True, help='Directory containing logits files')
    
    args = parser.parse_args()
    
    format_marg_plot_data(
        layer_data_tsv=args.layer_data_tsv,
        fasta_margins=args.fasta_margins,
        output_dir=args.output_dir,
        logits_dir=args.logits_dir
    )
