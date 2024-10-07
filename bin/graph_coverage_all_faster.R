#!/usr/bin/env Rscript
library(tidyverse)

args <- commandArgs(trailingOnly = TRUE)

# Check if there are at least 5 arguments (file paths, group size, low depth, cutoff)
if (length(args) < 5) {
  stop("Insufficient number of arguments. Usage: script.R file1 file2 ... group_size low_depth cutoff")
}

files <- args[1:(length(args) - 3)]
group_size <- strtoi(args[length(args) - 2])
low_depth <- strtoi(args[length(args) - 1])
cutoff <- strtoi(args[length(args)])

# Read all files into a single data frame
all_data <- map_dfr(files, read_tsv, col_names = c('id', 'locus', 'depth')) %>%
  mutate(sample = gsub("_depth.tsv", "", basename(files)))

# Filter unique ids
unique_ids_subset <- all_data %>%
  distinct(id, .keep_all = TRUE)

# Compact data
compact_data <- all_data %>%
  group_by(id, sample) %>%
  mutate(group = ceiling(row_number() / group_size)) %>%
  group_by(id, sample, group) %>%
  summarise(locus = mean(locus * group_size), depth = mean(depth)) %>%
  ungroup() %>%
  mutate(depth = pmin(depth, cutoff))

# Filter low depth data
low_data <- compact_data %>%
  filter(depth < low_depth)

# Get unique IDs and samples
ids <- unique(compact_data$id)
samples <- unique(compact_data$sample)

# Plotting
for (id in ids) {
  for (sample in samples) {
    filtered_data <- compact_data %>%
      filter(id == !!id, sample == !!sample)

    p <- ggplot(filtered_data, aes(x = locus, y = depth)) +
      geom_point(data = low_data, color = 'red', size = 0.1) +
      geom_line(linewidth = 0.4) +
      ggtitle(paste("ID: ", id, ", Sample: ", sample)) +
      theme_minimal()

    outfile <- paste0(id, "_", sample, "_depth.png")
    ggsave(outfile, plot = p, width = 10, height = 8, dpi = 600, limitsize = FALSE)
  }
}
