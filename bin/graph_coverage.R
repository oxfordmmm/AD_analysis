#!/usr/bin/env Rscript
library(tidyverse)

args <- commandArgs(trailingOnly = TRUE)
file <- args[1]
group_size <- strtoi(args[2]) # e.g. 50
low_depth <- strtoi(args[3]) # e.g. 5
cutoff <- strtoi(args[4]) # e.g. 100
outfile <- args[5]
sample <- args[6]

df <- read_tsv(file, col_names = c('id', 'locus', 'depth'))
bacterial_chroms=c('CP085971.1','NZ_CP025371.1','NC_005043.1','NZ_LR214945.1')

# filter bacterial chromosomes
df <- df %>% 
  filter(locus %in% bacterial_chroms)
df

ids <- df %>% 
  select(id) %>% 
  distinct() %>% 
  pull(id)

# Compresses data by grouping by group_size and taking the mean depth over each group
compact_df <- df %>% 
  group_by(id) %>% 
  mutate(group = ceiling(row_number() / group_size)) %>% 
  group_by(id, group) %>% 
  summarise(locus = mean(group) * group_size, depth = mean(depth)) %>% 
  ungroup() %>% 
  mutate(depth = pmin(depth, cutoff))

low_df <- df %>% 
  filter(depth < low_depth)

p <- compact_df %>%  
  ggplot(aes(x=locus, y=depth)) +
  ggtitle(paste("Sample: ", sample)) +
  geom_point(data=low_df, color='red', size=1) +
  geom_line(linewidth=1) +
  facet_wrap(vars(id), scales='free') +
  theme(strip.text = element_text(size = 5))

ggsave(outfile, plot=p, width=25, height=15, dpi=300, limitsize = F, format='png')
#ggsave(outfile, plot=p, width=25, height=15, dpi=300, limitsize = F, format='pdf')

