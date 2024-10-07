#!/usr/bin/env Rscript
library(tidyverse)
library(vroom)
library(gtools)

args <- commandArgs(trailingOnly = TRUE)
file <- args[1]
group_size <- strtoi(args[2]) # e.g. 50
low_depth <- strtoi(args[3]) # e.g. 5
cutoff <- strtoi(args[4]) # e.g. 100
outfile <- args[5]
sample <- args[6]

df <- read_tsv(file, col_names = c('id', 'locus', 'depth'))
head(df)

# Reorder the levels of id to be alphabetical and numerical (uses gtools)
# Do this to the original df so all downstream subsets are ordered the same
levels <- mixedsort(as.character(unique(df$id)))
df$id <- factor(df$id, levels = levels)

ids <- df %>% 
  select(id) %>% 
  distinct() %>% 
  pull(id) 

print("IDs:")
head(ids)

# Compresses data by grouping by group_size and taking the mean depth over each group
compact_df <- df %>% 
  group_by(id) %>% 
  mutate(group = ceiling(row_number() / group_size)) %>% 
  group_by(id, group) %>% 
  summarise(locus = mean(group) * group_size, depth = mean(depth)) %>% 
  ungroup() %>% 
  mutate(depth = pmin(depth, cutoff))

print("Compact df:")
head(compact_df)

low_df <- df %>% 
  filter(depth < low_depth)

print("Low df:")
head(low_df)

#### Plot for all ####
p <- compact_df %>%  
  ggplot(aes(x=locus, y=depth)) +
  ggtitle(paste("Sample: ", sample)) +
  geom_line(linewidth=1) +
  geom_point(data=low_df, color='red', size=1) +
  facet_wrap(vars(id), scales='free') +
  theme(strip.text = element_text(size = 3))
#p

ggsave(outfile, plot=p, width=25, height=15, dpi=300, limitsize = F)

print("Plotted for all")

#### Filter and plot for controls ####``
control_ids <- c("MS2", "armored_rna", "murine_respirovirus", 
  "orthoreovirus_1", "orthoreovirus_2", "orthoreovirus_3", "orthoreovirus_4", 
  "orthoreovirus_5", "orthoreovirus_6", "orthoreovirus_7", "orthoreovirus_8", 
  "orthoreovirus_9", "orthoreovirus_10", "zika")

controls_df <- compact_df %>% 
  filter(id %in% control_ids)

print("Controls df:")
head(controls_df)

low_controls_df <- df %>% 
  filter(id %in% control_ids) %>%
  filter(depth < low_depth)

print("Low controls df:")
head(low_controls_df)

p_controls <- controls_df %>% 
  ggplot(aes(x=locus, y=depth)) +
  ggtitle(paste("Sample: ", sample)) +
  geom_point(data = low_controls_df, color = 'red', size = 1) +
  geom_line(linewidth=1) +
  facet_wrap(vars(id), scales='free') +
  theme(strip.text = element_text(size = 10))

#p_controls
ggsave(paste0("controls_",outfile), plot=p_controls, width=25, height=15, dpi=300, limitsize = F)

print("P controls saved")

#### Filter and plot of positive samples ####
positive_df <- compact_df %>% 
  filter(!id %in% control_ids) %>%
  group_by(id) %>%
  filter(sum(depth) != 0) %>%
  ungroup()


print("Positive segments df:")
head(positive_df)

low_positives_df <- df %>% 
  filter(!id %in% control_ids) %>%
  group_by(id) %>%
  filter(sum(depth) != 0) %>%
  ungroup() %>%
  filter(depth < low_depth)

print("Low positive segments df:")
head(low_positives_df)

# Check if positive_df is empty
if (nrow(positive_df) > 0) {
  # Create plot only if positive_df is not empty
  p_positive <- positive_df %>% 
    ggplot(aes(x=locus, y=depth)) +
    ggtitle(paste("Sample: ", sample)) +
    geom_point(data=low_positives_df, color='red', size=1) +
    geom_line(linewidth=1) +
    facet_wrap(vars(id), scales='free') +
    theme(strip.text = element_text(size = 10))
  
  # Save the plot
  ggsave(paste0("positives_",outfile), plot=p_positive, width=25, height=15, dpi=300, limitsize = F)
  
  print("P positive saved")
} else {
  # Print a message if there are no positive segments
  print("No positive segments found in the sample.")
}

#### Filter and plot of positives and spikes ####
all_positives_df <- compact_df %>% 
  group_by(id) %>%
  filter(sum(depth) != 0) %>%
  ungroup()

low_all_positives_df <- df %>% 
  group_by(id) %>%
  filter(sum(depth) != 0) %>%
  ungroup() %>%
  filter(depth < low_depth)

p_all_positives <- all_positives_df %>% 
  ggplot(aes(x=locus, y=depth)) +
  ggtitle(paste("Sample: ", sample)) +
  geom_point(data=low_all_positives_df, color='red', size=1) +
  geom_line(linewidth=1) +
  facet_wrap(vars(id), scales='free') +
  theme(strip.text = element_text(size = 10))

p_all_positives

ggsave(paste0("all_positives_",outfile), plot=p_all_positives, width=25, height=15, dpi=300, limitsize = F)

p <- compact_df %>%  
  ggplot(aes(x=locus, y=depth)) +
  ggtitle(paste("Sample: ", sample)) +
  geom_point(data=low_df, color='red', size=1) +
  geom_line(linewidth=1) +
  facet_wrap(vars(id), scales='free') +
  theme(strip.text = element_text(size = 10))

ggsave(outfile, plot=p, width=25, height=15, dpi=300, limitsize = F)
