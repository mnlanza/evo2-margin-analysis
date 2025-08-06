#!/usr/bin/env Rscript

# Load required libraries
library(optparse)
library(reticulate)
library(Biostrings)
library(ggplot2)

source("scripts/utils.r")


# Command line argument parsing
option_list = list(
  make_option("--seq_data", type="character", help="Path to sequence data"),
  make_option("--fasta_file", type="character", help="Path to FASTA file with variants"),
  make_option("--aid", type="character", help="Sequence ID"),
  make_option("--job_dir", type="character", help="Directory containing the logits files")
)

# Parse command line arguments
opt <- parse_args(OptionParser(option_list=option_list))

# Print all command line arguments
print("Command line arguments received:")
print(opt)

# Read sequence data from TSV, forcing all columns to be character type
seq_data <- read.delim(opt$seq_data, header=TRUE, stringsAsFactors=FALSE, 
                      colClasses=c(coord="integer", 
                                 codon="character",
                                 aa="character",
                                 left_margin="integer",
                                 start="integer",
                                 end="integer",
                                 aid="character",
                                 mut_codon="character",
                                 mut_aa="character",
                                 synon_cod="character"))
print("Sequence data loaded:")
print(seq_data)

# Extract key parameters from sequence data
aa_coord <- seq_data$coord
src_codon <- seq_data$codon
mut_codon <- seq_data$mut_codon
syn_codon <- seq_data$synon_cod
left_margin <- seq_data$left_margin
aa_val <- seq_data$aa

# First read the FASTA file to get sequences
all_variants <- readDNAStringSet(opt$fasta_file)

# Convert amino acid position to nucleotide position
mut_pos_nt <- (aa_coord - 1) * 3 + 1
print(paste("Using amino acid position:", aa_coord))
print(paste("Converted to nucleotide position:", mut_pos_nt))

# Define variants with meaningful names
# Create variants with actual codons
variants <- c()
variants[seq_data$codon] <- seq_data$aa        # source codon
variants[seq_data$synon_cod] <- seq_data$aa    # synonymous codon  
variants[seq_data$mut_codon] <- seq_data$mut_aa # mutant codon
variants["TAA"] <- "Z"                         # stop codon

# Initialize list to store dataframes
variant_dfs <- list()

# Print debug info
print("Job directory:")
print(opt$job_dir)
print("Variants:")
print(variants)


# Define margins and ensure they're formatted as regular numbers
margins <- sprintf("%d", c(0, 20000, 40000, 60000, 80000, 100000))
for (margin in margins) {
  print(paste("Processing margin:", margin))
  # Iterate through variants and store results
  for (i in seq_along(variants)) {
    codon <- names(variants)[i]
    aa <- variants[i]
    
    print(paste("Processing variant:", i))
    print(paste("Codon:", codon))
    print(paste("AA:", aa))
    
    # Create name for this variant
    variant_name <- paste0(aa_coord, "_", aa, "_", codon, "_", margin)
    print(paste("Variant name:", variant_name))
    
    # Check if job_dir is provided
    if (is.null(opt$job_dir)) {
      stop("Error: --job_dir argument is required")
    }
    
    # Create file path
    evo2_npy_file <- file.path(opt$job_dir, paste0("input_", variant_name, "_logits.npy"))
    print(paste("Looking for file:", evo2_npy_file))
    
    # Initialize dataframe and store in list
    if (file.exists(evo2_npy_file)) {
      print("File found, initializing dataframe...")
      variant_dfs[[variant_name]] <- initialize_df(evo2_npy_file, variant_name, all_variants)
    } else {
      print(paste("Warning: Logits file not found:", evo2_npy_file))
    }
  }
}

# nucleotide positions to plot
nuc_plot_range <- 1:500

# Print total log-likelihoods for each variant
print("\nTotal log-likelihoods for each variant:")
for (variant_name in names(variant_dfs)) {
  total_ll <- sum(variant_dfs[[variant_name]]$log_likelihood)
  print(sprintf("%s: %.2f", variant_name, total_ll))
}
print("\n")

# Create base figures directory
base_figures_dir <- file.path("figures", opt$aid)

# Process each margin size
for (margin in margins) {
  # Create margin-specific directory
  margin_dir <- file.path(base_figures_dir, paste0("margin_", margin), paste0("full"))
  dir.create(margin_dir, recursive = TRUE, showWarnings = FALSE)
  
  # Process each variant type (src, mut, syn, stop)
  for (i in seq_along(variants)) {
    codon <- names(variants)[i]
    aa <- variants[i]
    variant_name <- paste0(aa_coord, "_", aa, "_", codon, "_", margin)
    
    # Skip if we don't have data for this variant
    if (!variant_name %in% names(variant_dfs)) {
      print(paste("Skipping", variant_name, "- no data available"))
      next
    }
    
    variant_data <- variant_dfs[[variant_name]]
    
    # Create entropy plot
    entropy_plot <- plot_entropy(
      gene_name = variant_name,
      gene_df = variant_data,
      index_rows = nuc_plot_range,
      aa_pos = aa_coord
    )
    ggsave(
      file.path(margin_dir, paste0("entropy_", aa, "_", codon, ".pdf")),
      entropy_plot,
      width = 8,
      height = 6
    )
    
    # Create log-likelihood plot
    loglik_plot <- plot_log_likelihood(
      gene_name = variant_name,
      gene_df = variant_data,
      index_rows = nuc_plot_range,
      aa_pos = aa_coord
    )
    ggsave(
      file.path(margin_dir, paste0("loglik_", aa, "_", codon, ".pdf")),
      loglik_plot,
      width = 8,
      height = 6
    )
    prob <- probability_correct_base(variant_name, variant_data, all_variants)
    print(paste("Probability of correct base:", prob))
    print(paste("Created plots for", variant_name))
  }
  print(paste("Completed margin size:", margin))
}

# Create output directories
stacked_dir <- file.path(base_figures_dir, "stacked")
compare_dir <- file.path(base_figures_dir, "compare")
total_ll_dir <- file.path(base_figures_dir, "total_loglik")

# Remove old directories if they exist
unlink(stacked_dir, recursive = TRUE)
unlink(compare_dir, recursive = TRUE)
unlink(total_ll_dir, recursive = TRUE)

# Create fresh directories
dir.create(stacked_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(compare_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(total_ll_dir, recursive = TRUE, showWarnings = FALSE)

print(paste("Created fresh directories:", stacked_dir, "and", compare_dir, "and", total_ll_dir))

# Create total log-likelihood plots for each variant
print("\nCreating total log-likelihood plots...")
# Get unique variant bases
variant_bases <- unique(sapply(names(variant_dfs), function(name) {
  parts <- strsplit(name, "_")[[1]]
  paste(parts[1:3], collapse="_")  # Keep coord, aa, and codon
}))

# Create single stacked plot with all variants
print("Creating stacked total log-likelihood plot for all variants")
total_ll_plot <- plot_total_loglik_margins(variant_dfs, variant_bases)

# Save plot
ggsave(
  file.path(total_ll_dir, "total_loglik_all_variants.pdf"),
  total_ll_plot,
  width = 10,  # Wider to accommodate labels
  height = 6 # Taller to accommodate stacked plots
)

# Create stacked plots for each variant
for (i in seq_along(variants)) {
  codon <- names(variants)[i]
  aa <- variants[i]
  variant_base <- paste0(aa_coord, "_", aa, "_", codon)
  
  # Create entropy stacked plot
  entropy_stacked <- plot_stacked_margins(
    variant_name = variant_base,
    variant_dfs = variant_dfs,
    aa_pos = aa_coord,
    index_rows = 1:3234,
    metric = "entropy"
  )
  ggsave(file.path(stacked_dir, paste0("entropy_", aa, "_", codon, ".pdf")), 
         entropy_stacked, width = 8, height = 6)
  
  # Create log-likelihood stacked plot
  loglik_stacked <- plot_stacked_margins(
    variant_name = variant_base,
    variant_dfs = variant_dfs,
    aa_pos = aa_coord,
    index_rows = 1:3234,
    metric = "log_likelihood"
  )
  ggsave(file.path(stacked_dir, paste0("loglik_", aa, "_", codon, ".pdf")), 
         loglik_stacked, width = 8, height = 6)
}

# Debug variant_dfs contents
print("Available variants in variant_dfs:")
print(names(variant_dfs))

# Create margin comparison plots for each variant
print("Creating entropy comparison plots...")
entropy_plots <- plot_margin_comparison(
  variant_dfs = variant_dfs,
  aa_pos = aa_coord,
  index_rows = nuc_plot_range,
  metric = "entropy"
)

print(paste("Got entropy plots for variants:", paste(names(entropy_plots), collapse=", ")))

# Save entropy plots
for (variant in names(entropy_plots)) {
  variant_parts <- strsplit(variant, "_")[[1]]
  aa <- variant_parts[2]
  codon <- variant_parts[3]
  ggsave(
    file.path(compare_dir, paste0("entropy_", aa, "_", codon, ".pdf")),
    entropy_plots[[variant]],
    width = 10,
    height = 6
  )
}

# Create and save log-likelihood plots
print("Creating log-likelihood comparison plots...")
loglik_plots <- plot_margin_comparison(
  variant_dfs = variant_dfs,
  aa_pos = aa_coord,
  index_rows = nuc_plot_range,
  metric = "log_likelihood"
)

print(paste("Got log-likelihood plots for variants:", paste(names(loglik_plots), collapse=", ")))

# Save log-likelihood plots
for (variant in names(loglik_plots)) {
  variant_parts <- strsplit(variant, "_")[[1]]
  aa <- variant_parts[2]
  codon <- variant_parts[3]
  ggsave(
    file.path(compare_dir, paste0("loglik_", aa, "_", codon, ".pdf")),
    loglik_plots[[variant]],
    width = 10,
    height = 6
  )
}
