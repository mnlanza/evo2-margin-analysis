#!/usr/bin/env python3
import sys
import csv
import argparse
from collections import defaultdict
from Bio import SeqIO
from Bio.Seq import Seq

def read_codon_table(codon_table_file):
    """read codon table and return codon->aa mapping"""
    with open(codon_table_file, 'r') as f:
        reader = csv.DictReader(f, delimiter='\t')
        return {row['codon']: row['aa'] for row in reader}

def read_fasta(fasta_file, target_id):
    """read sequence from fasta file by identifier"""
    sequences = SeqIO.parse(fasta_file, "fasta")
    
    for record in sequences:
        if record.id == target_id:
            return record.id, str(record.seq)
    
    # sequence not found
    print(f"error: sequence with identifier '{target_id}' not found in fasta file", file=sys.stderr)
    sys.exit(1)

def generate_all_codons():
    """generate all 64 possible codons"""
    bases = ['A', 'T', 'G', 'C']
    codons = []
    for b1 in bases:
        for b2 in bases:
            for b3 in bases:
                codons.append(b1 + b2 + b3)
    return codons

def find_synonymous_codon(codon_to_aa, original_codon):
    """Find a different codon that codes for the same amino acid"""
    original_aa = codon_to_aa.get(original_codon, 'X')
    # Get all codons that code for the same amino acid
    synonymous_codons = [cod for cod, aa in codon_to_aa.items() if aa == original_aa and cod != original_codon]
    # Return first synonymous codon if exists, otherwise return original
    return synonymous_codons[0] if synonymous_codons else original_codon

def main():
    parser = argparse.ArgumentParser(description='Generate codon variants for a specific position')
    parser.add_argument('--fasta', '-f', required=True, help='Input FASTA file')
    parser.add_argument('--codon-table', '-c', required=True, help='Codon table file')
    parser.add_argument('--aa-coord', '-a', type=int, required=True, help='Amino acid coordinate (1-based)')
    parser.add_argument('--seq-id', '-s', required=True, help='Sequence identifier')
    parser.add_argument('--output-fasta', '-o', required=True, help='Output FASTA file')
    parser.add_argument('--plot_info_table', '-v', required=True, help='Output codon file (original codon)')
    parser.add_argument('--left-margin', '-l', type=int, default=2000, help='Left margin')
    parser.add_argument('--right-margin', '-r', type=int, default=1000, help='Right margin')
    parser.add_argument('--gene-start', '-g', type=int, default=None, help='Gene start nucleotide position(1-based) from table')
    parser.add_argument('--gene-end', '-e', type=int, default=None, help='Gene end nucleotide position(1-based) from table')
    parser.add_argument('--aid', '-i', required=True, help='AID')
    parser.add_argument('--mut-codon', '-t', required=True, help='Target codon')

    args = parser.parse_args()
    
    # read codon table
    codon_to_aa = read_codon_table(args.codon_table)
    
    # read input fasta
    header, contig_seq = read_fasta(args.fasta, args.seq_id)

    print(f'contig length: {len(contig_seq)}')
    # making gene seq and checking gene nuc & aalengths
    gene_seq = contig_seq[args.gene_start-1:args.gene_end]
    print(f"sequence length: {len(gene_seq)} nt, {len(gene_seq) // 3} aa")

    # calculate nucleotide position (1-based aa coord to 0-based nucleotide)
    nt_start = (args.aa_coord - 1) * 3
    nt_end = nt_start + 3
    
    # check bounds
    if nt_end > len(gene_seq):
        print(f"error: position {args.aa_coord} is beyond sequence length", file=sys.stderr)
        sys.exit(1)

    # extract and print original codon as sanity check
    original_codon = gene_seq[(nt_start):(nt_end)]
    print(f"original codon: {original_codon}")
    original_aa = codon_to_aa.get(original_codon, 'X')
    mut_aa = codon_to_aa.get(args.mut_codon, 'X')
    
    # Find a synonymous codon
    synon_cod = find_synonymous_codon(codon_to_aa, original_codon)
    print(f"synonymous codon found: {synon_cod} (also codes for {original_aa})")

    print(f"original codon at position {args.aa_coord}: {original_codon} -> {original_aa}")

    # write original codon to output file
    print(f"writing original codon to {args.plot_info_table}")
    with open(args.plot_info_table, 'w') as out:
        out.write(f"coord\tcodon\taa\tleft_margin\tstart\tend\taid\tmut_codon\tmut_aa\tsynon_cod\n")
        out.write(f"{args.aa_coord}\t{original_codon}\t{original_aa}\t{args.left_margin}\t{args.gene_start}\t{args.gene_end}\t{args.aid}\t{args.mut_codon}\t{mut_aa}\t{synon_cod}\n")
    
    # get all possible codons
    all_codons = generate_all_codons()

    # # not implementing hard coded margins
    # left_margin = contig_seq[max(0, args.gene_start - 1 - args.left_margin):args.gene_start - 1]
    # right_margin = contig_seq[args.gene_end:min(len(contig_seq), args.gene_end + args.right_margin)]

    #### TESTING ####
    # variant_seq = gene_seq[:nt_start] + "XXX" + gene_seq[nt_end:]
    # print(f'here is the variant seq:{variant_seq} with length: {len(variant_seq)}')
    # print(f"gene seq: {gene_seq}")
    # print(f"Original seq: {contig_seq[:args.gene_start-1]}{gene_seq}{contig_seq[args.gene_end:]} and seq length {len(contig_seq[:args.gene_start-1]+gene_seq+contig_seq[args.gene_end:])}")
    # print(f"Margin seq: {left_margin}{gene_seq}{right_margin} and margin length: {len(left_margin+gene_seq+right_margin)}") 

    
    print(f"Generating file: {args.output_fasta}")
    # generate variants
    with open(args.output_fasta, 'w') as out:
        for codon in all_codons:
            aa = codon_to_aa.get(codon, 'X')  # X for unknown

            # replace stop codon symbol (*) with Z
            aa = 'Z' if aa == '*' else aa
            variant_seq = gene_seq[:nt_start] + codon + gene_seq[nt_end:]
            
            # First write the gene sequence without margins
            out.write(f">{args.aa_coord}_{aa}_{codon}_0\n")
            out.write(f"{variant_seq}\n")
            
            # Then write sequences with different margins
            for n in range(1, 6):  # Start from 1 since we already wrote 0 margin
                current_margin_size = n * args.left_margin
                # Get margin from contig sequence
                left_margin = contig_seq[max(0, args.gene_start - 1 - current_margin_size):args.gene_start - 1]
                margin_var_seq = left_margin + variant_seq
                out.write(f">{args.aa_coord}_{aa}_{codon}_{current_margin_size}\n")
                out.write(f"{margin_var_seq}\n")



if __name__ == "__main__":
    main()