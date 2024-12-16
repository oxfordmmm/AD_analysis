#!/usr/bin/env Rscript
library(tidyverse) 

args <- commandArgs(trailingOnly = TRUE)

# Check if there are at least 5 arguments (file paths, group size, low depth, cutoff)
if (length(args) < 5) {
  stop("Insufficient number of arguments. Usage: script.R file1 file2 ... group_size low_depth cutoff")
}

files <- args[1:(length(args)-3)]

group_size <- strtoi(args[length(args)-2])
low_depth <- strtoi(args[length(args)-1])
cutoff <- strtoi(args[length(args)])

cat("Input files:\n")
cat(head(files), sep = "\n")

cat("group_size:\n")
cat(group_size, sep = "\n")

cat("low_dept:\n")
cat(low_depth , sep = "\n")

cat("cutoff:\n")
cat(cutoff, sep = "\n")

all_data <- data.frame()

bacterial_chroms=c('CP085971.1','NZ_CP025371.1','NC_005043.1','NZ_LR214945.1')

for (file in files) {
  # Read data from the current file
  depth_df <- read_tsv(file, col_names = c('id', 'locus', 'depth'))
  # filter bacterial chromosomes
  depth_df <- depth_df %>% 
     filter(locus %in% bacterial_chroms)
  depth_df

  # Extract the sample name from the file name
  sample_name <- gsub("_depth.tsv", "", basename(file))

  # Add a "sample" column to df
  # get rid of the id2 - it is wrong
  depth_df <- depth_df %>% mutate(sample = sample_name)
  # ,
  #                                 id2 = case_when(
  #                                   str_detect(sample, "(A|B|C)") ~
  #                                   paste0("RV", str_extract(sample, "(A|B|C)"), "_", id),
  #                                   TRUE ~ id
  #                                   )
  #                                 )

  # Append data to the all_data data frame
  all_data <- bind_rows(all_data, depth_df)
}

cat("All data filtered by unique ids:\n")

unique_ids_subset <- all_data %>%
  distinct(id, .keep_all = TRUE)

head(unique_ids_subset, 5)
tail(unique_ids_subset, 5)

all_data <- all_data %>% 
  mutate(sample_id = paste(sample,id,sep = "_"))

# Get all samples
ids <- all_data %>% 
  select(id) %>% 
  distinct() %>% 
  pull(id)

# Manipulate df
compact_data <- all_data %>% 
group_by(sample_id) %>% 
mutate(group = ceiling(row_number() / group_size)) %>% 
group_by(sample_id, group) %>% 
summarise(locus = mean(group) * group_size, depth = mean(depth), sample = sample, id = id) %>% 
ungroup() %>% 
mutate(depth = pmin(depth, cutoff)) %>%
select(id, sample, group, locus, depth)


low_data <- all_data %>% 
filter(depth < low_depth)

nsamples <- length((unique(compact_data$sample)))
print("number of samples:")
print(nsamples)
head(unique(compact_data$sample))


# Plot a graph per id
for (id in ids) {
  filtered_data <- compact_data %>% filter(id == !!id)

  p <- filtered_data %>%
    ggplot(aes(x = locus, y = depth)) +
    geom_point(data = low_data, color = 'red', size = 0.1) +
    geom_line(linewidth = 0.3) +
    ggtitle(paste("ID: ", id)) +
    facet_wrap(vars(sample), scales = 'fixed') +
    theme(strip.text = element_text(size = 5))
  
  print("saving...")
  print(id)
  outfile <- paste0(id, "_depth.png")
  ggsave(outfile, plot = p, width = 10, height = 8, dpi = 600, limitsize = FALSE) 
}