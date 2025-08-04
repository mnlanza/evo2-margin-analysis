# utils.py
import numpy as np
from Bio import SeqIO
import math

def softmax(logits_vec):
    """Takes a vector of logits (real numbers) and transforms into probability distribution"""
    # Subtracting max for numerical stability
    exp_logits = np.exp(logits_vec - np.max(logits_vec))
    return exp_logits / np.sum(exp_logits)

def get_logits(filename):
    """Load and process logits from numpy file"""
    # Load numpy file
    base_np = np.load(filename)
    
    # Take first slice (equivalent to R's [1, , ])
    mat = base_np[0]
    
    # Apply softmax row-wise
    all_probs = np.apply_along_axis(softmax, 1, mat)
    
    # Extract nucleotide probabilities
    # Note: ASCII values for 'A', 'C', 'G', 'T' are 65, 67, 71, 84
    nuc_prob = {
        'A': all_probs[:, ord('A')],
        'C': all_probs[:, ord('C')],
        'G': all_probs[:, ord('G')],
        'T': all_probs[:, ord('T')]
    }
    
    # Create empty row for position 1 (uniform distribution)
    empty_row = {'A': 0.25, 'C': 0.25, 'G': 0.25, 'T': 0.25}
    
    # Convert to numpy array and add empty row at start
    nuc_array = np.array([nuc_prob['A'], nuc_prob['C'], nuc_prob['G'], nuc_prob['T']]).T
    empty_row_array = np.array([[0.25, 0.25, 0.25, 0.25]])
    nuc_array = np.vstack([empty_row_array, nuc_array])
    
    # Remove last row and return
    return nuc_array[:-1]

def initialize_dict(EVO2_npy_file, sequence_name, all_variants):
    """Initialize dataframe with probabilities, entropy, and log-likelihood"""
    # Load and normalize logits
    prob_matrix = get_logits(EVO2_npy_file)
    
    # Calculate entropy
    entropy = np.sum(-prob_matrix * np.log2(prob_matrix + 1e-10), axis=1)
    
    # Get sequence and convert to nucleotide list
    desired_seq = all_variants[sequence_name]
    seq_vec = list(desired_seq)
    
    # Calculate log-likelihood
    nuc_to_idx = {'A': 0, 'C': 1, 'G': 2, 'T': 3}
    log_likelihood = []
    for i, base in enumerate(seq_vec):
        base_idx = nuc_to_idx[base]
        log_likelihood.append(math.log2(prob_matrix[i, base_idx]))
    
    # Create result dictionary
    result = {
        'prob_matrix': prob_matrix,
        'entropy': entropy,
        'log_likelihood': np.array(log_likelihood)
    }
    
    return result