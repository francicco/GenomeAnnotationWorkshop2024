library(ggplot2)
library(tidyverse)


# Parse command-line arguments
args <- commandArgs(trailingOnly = TRUE)

# Check if the correct number of arguments is provided
if (length(args) != 2) {
  stop("Usage: Rscript Analyze_Diamond_topHit_coverage.R w_pct_hit_length.input_file file.png")
}

# Extract input and output file paths from the arguments
input_file <- args[1]
output_file <- args[2]

AlnFrac <- read_tsv(input_file, col_names = TRUE) %>%
  select(1:15)

head(AlnFrac)

Density_plot <- ggplot(AlnFrac, aes(pct_hit_len_aligned, after_stat(count), color = Annotation, linetype = Annotation)) +
  geom_density(size=0.3) +
  scale_x_reverse() +
  labs(title=paste('UniProt comparison'),
       x='Percentage of AA length aligned', y='Count') 

ggsave(output_file, plot = Density_plot, width = 10, height = 6)
