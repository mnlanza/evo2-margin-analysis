suppressPackageStartupMessages({
library(reticulate)
library(Biostrings)
library(ggplot2)
np <- import("numpy")
})

# This script normalizes EVO2 input and adds a row to account for
# EVO2 not having prediction data for the first position (nuc)

#' Create stacked bar plots showing total log-likelihood for each margin size
#' 
#' @param variant_dfs List of dataframes containing variant data
#' @param variant_bases List of variant base names
#' @return ggplot object
plot_total_loglik_margins <- function(variant_dfs, variant_bases) {
  # Initialize data frame for plotting
  plot_data <- data.frame()
  
  # Get data for all variants and margins
  margins <- c(0, 20000, 40000, 60000, 80000, 100000)
  for (variant_base in variant_bases) {
    variant_parts <- strsplit(variant_base, "_")[[1]]
    aa <- variant_parts[2]
    codon <- variant_parts[3]
    variant_label <- paste0(aa, "_", codon)
    
    for (margin in margins) {
      df_name <- paste0(variant_base, "_", sprintf("%d", margin))
      if (df_name %in% names(variant_dfs)) {
        total_ll <- sum(variant_dfs[[df_name]]$log_likelihood)
        plot_data <- rbind(plot_data, 
                          data.frame(margin = paste0(margin/1000, "kb"),
                                   total_ll = total_ll,
                                   variant = variant_label))
      }
      print(paste("Checking", df_name, ":", df_name %in% names(variant_dfs)))
    }
  }
  
  # Convert margin to factor with specific order
  plot_data$margin <- factor(plot_data$margin, 
                           levels = paste0(c(0, 20, 40, 60, 80, 100), "kb"))
  
  # Convert variant to factor to maintain order
  plot_data$variant <- factor(plot_data$variant, levels = unique(plot_data$variant))
  
  # For each variant, calculate range to set y-axis limits
  variant_ranges <- by(plot_data, plot_data$variant, function(d) {
    min_val <- min(d$total_ll)
    max_val <- max(d$total_ll)
    range <- max_val - min_val
    # Extend range by 10% on each end
    list(
      min = min_val - range * 0.1,
      max = max_val + range * 0.1
    )
  })

  # Create stacked bar plot with custom y-axis ranges for each variant
  variants <- levels(plot_data$variant)
  plots <- lapply(seq_along(variants), function(i) {
    v <- variants[i]
    is_last <- i == length(variants)  # Check if this is the last plot
    
    variant_data <- plot_data[plot_data$variant == v,]
    y_min <- min(variant_data$total_ll)
    y_max <- max(variant_data$total_ll)
    y_range <- y_max - y_min
    y_min <- y_min - y_range * 0.1  # Add 10% padding
    y_max <- y_max + y_range * 0.1
    
    ggplot(variant_data, aes(x = margin, y = total_ll)) +
      geom_bar(stat = "identity", fill = "steelblue", alpha = 0.7) +
      coord_cartesian(ylim = c(y_min, y_max)) +  # Set custom y-axis range
      labs(
        title = v,
        x = if (is_last) "margin" else ""  # Only show x label on last plot
      ) +
      theme_minimal() +
      theme(
        # X-axis settings
        axis.text.x = if (is_last) element_text(angle = 45, hjust = 1, size = 8) else element_blank(),
        axis.title.x = if (is_last) element_text() else element_blank(),
        # Y-axis settings
        axis.text.y = element_text(size = 7),
        # Margins - only bottom plot gets margin for x-axis labels
        plot.margin = if (is_last) {
          margin(b = 20, t = 0, l = 5, r = 5)
        } else {
          margin(b = 0, t = 0, l = 5, r = 5)
        }
      )
  })
  
  # Combine plots using patchwork
  library(patchwork)
  p <- wrap_plots(plots, ncol = 1, heights = rep(1, length(plots))) +
    plot_layout(heights = unit(rep(1, length(plots)), "null")) +
    plot_annotation(
      title = "Total log-likelihood across margins for all variants",
      caption = "Note: Y-axes do not start at zero to better show differences",
      theme = theme(
        plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
        plot.caption = element_text(hjust = 0, face = "italic")
      )
    ) &
    theme(axis.title.x = element_text(size = 12),
          axis.title.y = element_text(size = 12),
          axis.text = element_text(size = 10))
    labs(x = "Margin size",
         y = "Total log-likelihood",
         title = "Total log-likelihood across margins for all variants",
         caption = "Note: Y-axes do not start at zero to better show differences") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          plot.title = element_text(size=14, face="bold", hjust = 0.5),
          axis.title = element_text(size=12),
          axis.text = element_text(size=10),
          strip.text.y = element_text(angle = 0),  # Make variant labels horizontal
          panel.spacing = unit(1, "lines"),  # Add space between plots
          plot.caption = element_text(hjust = 0, face = "italic"))
  
  return(p)
}

# takes a vector of logits (real #'s) and transforms into prob distr
softmax <- function(logits_vec) {
  # subtracting each element from its max (for numerical stability)
  exp_logits_vec <- exp(logits_vec - max(logits_vec))
  exp_logits_vec / sum(exp_logits_vec) # normalizing vector to sum to 1
}

#
get_logits <- function(filename) {
  base_np <- np$load(filename) # loading the function and putting into base_np

  # convert numpy file to Rx
  mat <- py_to_r(base_np)[1, , ] # python to r array (takes first slice)

  # turn into probability matrix using softmax (N samples x M tokens)
  # applies 3rd arg to 1st arg (1 = row wise) & t is transpose
  all <- t(apply(mat, 1, softmax))
  # focus only on nucleotides 
  nuc_prob <- data.frame(
    A = all[, utf8ToInt("A") + 1], # utf8ToInt ASCII + of "A" and then adding 1
    C = all[, utf8ToInt("C") + 1],
    G = all[, utf8ToInt("G") + 1],
    T = all[, utf8ToInt("T") + 1]
  )

  # shifting EVO2 data: EVO2 data starts at pos 2 (FASTA data starts at pos 1)
  empty_row <- data.frame(A = 0.25, C = 0.25, G = 0.25, T = 0.25) # at pos 1
  nuc_prob <- rbind(empty_row, nuc_prob)
  nuc_prob[-nrow(nuc_prob), ]
}

initialize_df <- function(EVO2_npy_file, sequence_name, all_variants) { 
  # Load and normalize logits into nucleotide probabilities
  prob_matrix <- get_logits(EVO2_npy_file) 

  # Calculate entropy and add it as a column
  prob_matrix$entropy <- apply(prob_matrix[,1:4], 1, 
    function(position_prob) -sum(position_prob * log2(position_prob)))
  
  # Get sequences and find the no-margin version for comparison
  seqs <- as.character(all_variants)
  # Get the base sequence name (without margin) by removing the last part after last underscore
  base_seq_name <- sub("_[^_]+$", "_0", sequence_name)
  # Get the no-margin sequence for comparison
  base_seq <- seqs[base_seq_name]
  base_seq_vec <- strsplit(base_seq, split = "")[[1]]
  
  # Calculate log-likelihood using the no-margin sequence
  log_likelihood <- mapply(function(base, i) {
      base_index <- match(base, c("A", "C", "G", "T"))
      log2(prob_matrix[i, base_index])}, base_seq_vec, seq_along(base_seq_vec))
  
  prob_matrix$log_likelihood <- log_likelihood
  return(prob_matrix)
}

build_plot_df <- function(gene_df, metric, index_rows = 200:300) {
  # Create dataframe with:
  # - pos: nucleotide positions we want to plot
  # - value: corresponding metric values (entropy or log-likelihood)
  # - codon_i: which position in codon (1,2,3) each nucleotide represents
  df <- data.frame(
    pos     = index_rows,
    value   = gene_df[[metric]][index_rows],
    codon_i = ((index_rows - 1) %% 3) + 1
  )
  df
}

#' Plot entropy across positions
#' 
#' @param gene_name Name of the gene variant
#' @param gene_df Dataframe containing entropy data
#' @param index_rows Range of positions to plot
#' @param highlighted Positions to highlight
#' @return ggplot object
plot_entropy <- function(gene_name, gene_df, index_rows = 200:300, aa_pos, 
                         color_by_codon = TRUE, point_color = "purple") {
  # Convert amino acid position to nucleotide positions
  mut_pos_nt <- (aa_pos - 1) * 3 + 1
  
  plot_df <- build_plot_df(gene_df, "entropy", index_rows)
  # Calculate plot limits from index_rows
  x_min <- min(index_rows)
  x_max <- max(index_rows)
  
  p <- ggplot(plot_df, aes(x = pos, y = value)) + 
    coord_cartesian(xlim = c(x_min, x_max))
    
  # Add points with or without codon coloring
  if (color_by_codon) {
    p <- p + geom_point(aes(color = factor(codon_i)), size = 1.5, alpha = 0.7) +
      scale_color_manual(name = "Codon position",
                        values = c("1" = "red", "2" = "orange", "3" = "green"))
  } else {
    p <- p + geom_point(color = point_color, size = 1.5, alpha = 0.7)
  }
  
  p <- p + 
    # Add vertical lines for mutation position
    # Add vertical lines for mutation position (using purple for better visibility)
    geom_vline(xintercept = mut_pos_nt, linetype="dashed", color="purple", linewidth=1, alpha=0.3) +
    geom_vline(xintercept = mut_pos_nt + 1, linetype="dashed", color="purple", linewidth=1, alpha=0.3) +
    geom_vline(xintercept = mut_pos_nt + 2, linetype="dashed", color="purple", linewidth=1, alpha=0.3) +
    labs(x = "Nucleotide position", y = "Entropy",
         title = sprintf("Entropy across positions %d-%d (%s)", 
                        min(index_rows), max(index_rows), gene_name),
         subtitle = sprintf("Mutation at amino acid position %d (nucleotides %d-%d)",
                          aa_pos, mut_pos_nt, mut_pos_nt + 2)) +
    theme_minimal() + 
    theme(legend.position = "right", 
          legend.box = "vertical",
          plot.title = element_text(size=14, face="bold"),
          plot.subtitle = element_text(size=10),
          axis.title = element_text(size=12),
          axis.text = element_text(size=10))
}

#' Plot log-likelihood across positions
#' 
#' @param gene_name Name of the gene variant
#' @param gene_df Dataframe containing log-likelihood data
#' @param index_rows Range of positions to plot
#' @param highlighted Positions to highlight
#' @return ggplot object
plot_log_likelihood <- function(gene_name, gene_df, index_rows = 200:300, aa_pos,
                                color_by_codon = TRUE, point_color = "purple") {
  # Convert amino acid position to nucleotide positions
  mut_pos_nt <- (aa_pos - 1) * 3 + 1
  
  plot_df <- build_plot_df(gene_df, "log_likelihood", index_rows)
  # Calculate plot limits from index_rows
  x_min <- min(index_rows)
  x_max <- max(index_rows)
  
  # Initialize base plot
  p <- ggplot(plot_df, aes(x = pos, y = -value)) +
    coord_cartesian(xlim = c(x_min, x_max))
    
  # Add points with or without codon coloring
  if (color_by_codon) {
    p <- p + geom_point(aes(color = factor(codon_i)), size = 1.5, alpha = 0.7) +
      scale_color_manual(name = "Codon position",
                        values = c("1" = "red", "2" = "orange", "3" = "green"))
  } else {
    p <- p + geom_point(color = point_color, size = 1.5, alpha = 0.7)
  }
  
  p <- p +
    # Add vertical lines for mutation position
    # Add vertical lines for mutation position (using purple for better visibility)
    geom_vline(xintercept = mut_pos_nt, linetype="dashed", color="purple", linewidth=1, alpha=0.3) +
    geom_vline(xintercept = mut_pos_nt + 1, linetype="dashed", color="purple", linewidth=1, alpha=0.3) +
    geom_vline(xintercept = mut_pos_nt + 2, linetype="dashed", color="purple", linewidth=1, alpha=0.3) +
    labs(x = "Nucleotide position", y = "Log-likelihood",
         title = sprintf("-log-likelihood across positions %d-%d (%s)",
                        min(index_rows), max(index_rows), gene_name),
         subtitle = sprintf("Mutation at amino acid position %d (nucleotides %d-%d)",
                          aa_pos, mut_pos_nt, mut_pos_nt + 2)) +
    theme_minimal() + 
    theme(legend.position = "right", 
          legend.box = "vertical",
          plot.title = element_text(size=14, face="bold"),
          plot.subtitle = element_text(size=10),
          axis.title = element_text(size=12),
          axis.text = element_text(size=10))
}

#' Calculate probability of correct base prediction

#' @param gene Gene variant name
#' @param df_gene Dataframe with predictions
#' @param all_variants DNAStringSet with reference sequences
#' @return String with probability result
probability_correct_base <- function(gene, df_gene, all_variants) {
  seqs <- as.character(all_variants)
  # Get the base sequence name (without margin) by removing the last part after last underscore
  base_seq_name <- sub("_[^_]+$", "_0", gene)
  # Get the no-margin sequence for comparison
  base_seq <- seqs[base_seq_name]
  base_seq_vec <- strsplit(base_seq, split = "")[[1]]


  correct_count <- 0
  
  for (i in 1:length(base_seq_vec)) {
    max_column <- which.max(df_gene[i, 1:4])
    predicted_base <- c("A", "C", "G", "T")[max_column]
    if (predicted_base == base_seq_vec[i]) {
      correct_count <- correct_count + 1  
    }
  }
  answer <- correct_count/length(base_seq_vec)
  sprintf("The probability that EVO2 predicts the correct base across the gyrase gene in e-coli for %s is %.4f",
          gene, answer)
}

#' Create stacked plots for a variant at different margins
#' 
#' @param variant_name Name of the variant
#' @param variant_dfs List of dataframes containing variant data
#' @param aa_pos Amino acid position
#' @param index_rows Range of positions to plot
#' @param metric Which metric to plot ("entropy" or "log_likelihood")
#' @return ggplot object
plot_stacked_margins <- function(variant_name, variant_dfs, aa_pos, index_rows = 300:400,
                               metric = "entropy") {
  # Get data for margins 
  margins <- c(0, 4000, 8000, 20000)
  plot_data <- data.frame()
  
  for (margin in margins) {
    df_name <- paste0(variant_name, "_", margin)
    if (!df_name %in% names(variant_dfs)) {
      stop(paste("Data not found for", df_name))
    }
    
    df <- variant_dfs[[df_name]]
    subset <- build_plot_df(df, metric, index_rows)
    subset$margin <- paste0(margin/1000, "kb")  # Convert to kb for display
    plot_data <- rbind(plot_data, subset)
  }
  
  # Convert margin to factor with specific order
  plot_data$margin <- factor(plot_data$margin, 
                           levels = c("0kb", "4kb", "8kb", "20kb"))
  
  # Create base plot
  p <- ggplot(plot_data, aes(x = pos)) +
    coord_cartesian(xlim = c(min(index_rows), max(index_rows))) +
    facet_grid(margin ~ ., scales = "free_y")  # Stack plots vertically
    
  # Define colors for each margin level
  margin_colors <- c("0kb" = "#1f77b4",    # blue
                    "4kb" = "#2ca02c",    # green
                    "8kb" = "#ff7f0e",    # orange
                    "20kb" = "#d62728")   # red
  
  # Add points with different colors per margin
  p <- p + geom_point(aes(y = value, color = margin), size = 1, alpha = 0.7) +
       scale_color_manual(values = margin_colors)
  
  # Add mutation position lines
  mut_pos_nt <- (aa_pos - 1) * 3 + 1
  p <- p +
    geom_vline(xintercept = mut_pos_nt, linetype="solid", color="yellow", linewidth=1, alpha=0.5) +
    geom_vline(xintercept = mut_pos_nt + 1, linetype="solid", color="yellow", linewidth=1, alpha=0.5) +
    geom_vline(xintercept = mut_pos_nt + 2, linetype="solid", color="yellow", linewidth=1, alpha=0.5) +
    labs(x = "Nucleotide position",
         y = ifelse(metric == "entropy", "Entropy", "Log-likelihood"),
         title = sprintf("%s across positions %d-%d\nVariant: %s",
                        ifelse(metric == "entropy", "Entropy", "Log-likelihood"),
                        min(index_rows), max(index_rows), variant_name),
         subtitle = sprintf("Mutation at amino acid position %d (nucleotides %d-%d)",
                          aa_pos, mut_pos_nt, mut_pos_nt + 2)) +
    theme_minimal() +
    theme(legend.position = "right",
          strip.text.y = element_text(angle = 0),  # Margin labels horizontal
          plot.title = element_text(size=14, face="bold"),
          plot.subtitle = element_text(size=10),
          axis.title = element_text(size=12),
          axis.text = element_text(size=10))
  
  return(p)
}

#' Create line plots comparing margins for all variants
#' 
#' @param variant_dfs List of dataframes containing variant data
#' @param aa_pos Amino acid position
#' @param index_rows Range of positions to plot
#' @param metric Which metric to plot ("entropy" or "log_likelihood")
#' @return ggplot object
plot_margin_comparison <- function(variant_dfs, aa_pos, index_rows = 300:400,
                                 metric = "entropy") {
  # Get data for margins 0 and 100k
  margins <- c(0, 4000)
  plots <- list()  # Store plots for each variant
  
  # Extract variant base names (without margins)
  # Get unique variant bases by removing the margin part
  variant_bases <- unique(sapply(names(variant_dfs), function(name) {
    parts <- strsplit(name, "_")[[1]]
    paste(parts[1:3], collapse="_")  # Keep coord, aa, and codon
  }))
  
  print("Found variant bases:")
  print(variant_bases)
  
  # Create a plot for each variant
  for (variant in variant_bases) {
    print(paste("Processing variant:", variant))
    plot_data <- data.frame()
    
    # Get data for both margins
    margin_data <- list()
    for (margin in margins) {
      # Format margin as regular number, not scientific notation
      df_name <- paste0(variant, "_", sprintf("%d", margin))
      print(paste("  Checking for:", df_name))
      if (df_name %in% names(variant_dfs)) {
        print(paste("    Found data for:", df_name))
        df <- variant_dfs[[df_name]]
        subset <- build_plot_df(df, metric, index_rows)
        subset$margin <- paste0(margin/1000, "kb")
        margin_data[[as.character(margin)]] <- subset
      } else {
        print(paste("    Missing data for:", df_name))
      }
    }
    
    print(paste("  Found data for", length(margin_data), "margins"))
    # Only proceed if we have data for both margins
    if (length(margin_data) == 2) {
      plot_data <- do.call(rbind, margin_data)
      
      # Create line plot
      mut_pos_nt <- (aa_pos - 1) * 3 + 1
      
      # Extract variant info
      variant_parts <- strsplit(variant, "_")[[1]]
      aa <- variant_parts[2]
      codon <- variant_parts[3]
      
      # Create plot for this variant
      p <- ggplot(plot_data, aes(x = pos, y = value, color = margin)) +
        coord_cartesian(xlim = c(min(index_rows), max(index_rows))) +
        # Add lines with good visibility
        geom_line(linewidth = 0.5, alpha = 0.7) +
      geom_vline(xintercept = mut_pos_nt, linetype="solid", color="yellow", linewidth=1, alpha=0.3) +
      geom_vline(xintercept = mut_pos_nt + 1, linetype="solid", color="yellow", linewidth=1, alpha=0.3) +
      geom_vline(xintercept = mut_pos_nt + 2, linetype="solid", color="yellow", linewidth=1, alpha=0.3) +
      scale_color_manual(values = c("0kb" = "blue", "4kb" = "red")) +
      labs(x = "Nucleotide position",
           y = ifelse(metric == "entropy", "Entropy", "Log-likelihood"),
           title = sprintf("%s comparison (AA: %s, Codon: %s)",
                          ifelse(metric == "entropy", "Entropy", "Log-likelihood"),
                          aa, codon),
           subtitle = sprintf("Mutation at amino acid position %d (nucleotides %d-%d)",
                            aa_pos, mut_pos_nt, mut_pos_nt + 2)) +
      theme_minimal() +
      theme(legend.position = "right",
            plot.title = element_text(size=14, face="bold"),
            plot.subtitle = element_text(size=10),
            axis.title = element_text(size=12),
            axis.text = element_text(size=10))
      
      plots[[variant]] <- p
    }
  }
  
  return(plots)  # Return list of plots, one per variant
}

#' Plot metric differences between variants
#' 
#' @param ref_gene Reference gene name
#' @param target_gene Target gene name
#' @param ref_df Reference gene dataframe
#' @param target_df Target gene dataframe
#' @param metric Metric to compare ("entropy" or "log_likelihood")
#' @param index_rows Range of positions to plot
#' @param highlight Positions to highlight
#' @return ggplot object
plot_metric_diff <- function(ref_gene, target_gene, ref_df, target_df, 
                           metric = "entropy", index_rows = 200:300, aa_pos,
                           color_by_codon = TRUE, point_color = "purple") {
  # Convert amino acid position to nucleotide positions
  mut_pos_nt <- (aa_pos - 1) * 3 + 1
  
  diff_vec <- target_df[[metric]] - ref_df[[metric]]
  df <- data.frame(
    pos = index_rows,
    diff_vec = diff_vec[index_rows],
    codon_i = ((index_rows - 1) %% 3) + 1
  )
  
  # Calculate plot limits from index_rows
  x_min <- min(index_rows)
  x_max <- max(index_rows)
  
  # Initialize base plot
  p <- ggplot(df, aes(x = pos, y = diff_vec)) +
    coord_cartesian(xlim = c(x_min, x_max))
    
  # Add points with or without codon coloring
  if (color_by_codon) {
    p <- p + geom_point(aes(color = factor(codon_i)), size = 1.5, alpha = 0.7) +
      scale_color_manual(name = "Codon position",
                        values = c("1" = "red", "2" = "orange", "3" = "green"))
  } else {
    p <- p + geom_point(color = point_color, size = 1.5, alpha = 0.7)
  }
  
  p <- p +
    # Add vertical lines for mutation position
    # Add vertical lines for mutation position (using purple for better visibility)
    geom_vline(xintercept = mut_pos_nt, linetype="dashed", color="purple", linewidth=1, alpha=0.3) +
    geom_vline(xintercept = mut_pos_nt + 1, linetype="dashed", color="purple", linewidth=1, alpha=0.3) +
    geom_vline(xintercept = mut_pos_nt + 2, linetype="dashed", color="purple", linewidth=1, alpha=0.3) +
    scale_color_manual(name = "Codon position",
                      values = c("1" = "red", "2" = "orange", "3" = "green")) +
    labs(x = "Nucleotide position",
         y = sprintf("Δ %s", gsub("_", " ", metric, fixed = TRUE)),
         title = sprintf("%s difference: %s vs %s",
                        gsub("_", " ", metric, fixed = TRUE),
                        target_gene, ref_gene),
         subtitle = sprintf("Mutation at amino acid position %d (nucleotides %d-%d)",
                          aa_pos, mut_pos_nt, mut_pos_nt + 2)) +
    theme_minimal() + 
    theme(legend.position = "right", 
          legend.box = "vertical",
          plot.title = element_text(size=14, face="bold"),
          plot.subtitle = element_text(size=10),
          axis.title = element_text(size=12),
          axis.text = element_text(size=10))
}