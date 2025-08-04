#!/usr/bin/env Rscript
library(tidyverse)

# Read data
entropy_data <- read.delim("entropy.tsv", comment.char="#")
loglik_data <- read.delim("loglik.tsv", comment.char="#")

# Example plot
entropy_data %>%
  pivot_longer(starts_with("pos_"), names_to="position", values_to="entropy") %>%
  ggplot(aes(x=position, y=entropy, color=variant)) +
  geom_line() +
  facet_wrap(~margin_size)