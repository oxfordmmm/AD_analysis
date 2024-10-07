#!/usr/bin/env Rscript
library(tidyverse)

args <- commandArgs(trailingOnly = TRUE)
blast_file <- args[1]
kaiju_file <- args[2]
k2_file <- args[3]

murine_id <- 11191
orthoreovirus_id <- 538123
zika <- 64320

df <- read_tsv(blast_file) %>% 
  mutate(forward = sstart < send) %>% 
  mutate(low = pmin(sstart, send), high = pmax(sstart, send)) %>% 
  arrange(sseqid, low) %>% 
  group_by(sseqid) %>% 
  mutate(diff = low - lag(high) - 1) %>% 
  mutate(diff = pmax(0, diff)) %>% 
  ungroup() %>% 
  print(n=200)

kaiju_df <- read_tsv(kaiju_file, col_names = c('kaiju_code', 'sseqid', 'kaiju_taxon')) %>% 
  mutate(kaiju_taxon = str_remove(kaiju_taxon, "\t.*")) %>% 
  print()
kraken2_df <- read_tsv(k2_file, col_names = c('k2_code', 'sseqid', 'k2_taxon', 'read_length')) %>% 
  select(k2_code, sseqid, k2_taxon, read_length) %>% 
  print()
classification <- inner_join(kaiju_df, kraken2_df) %>% 
  mutate(taxon = pmax(k2_taxon, kaiju_taxon), code = case_when(
    taxon == 11191 ~ 'murine',
    taxon == 538123 ~ 'orthoreovirus',
    taxon == 64320 ~ 'zika',
    taxon > 1 ~ 'C',
    taxon == 1 ~ 'V',
    T ~ 'U'
  )) %>% 
  mutate(code = factor(code, levels=c('murine', 'orthoreovirus', 'zika', 'C', 'V', 'U'))) %>% 
  print(n=100) %>% 
  select(sseqid, code, taxon, read_length)

all <- classification %>% 
  left_join(df) %>% 
  mutate(primer_hit = !is.na(qseqid)) %>% 
  select(sseqid, code, taxon, read_length, primer_hit) %>% 
  distinct()

all %>% 
  ggplot(aes(x=code, fill=primer_hit)) +
  geom_bar()
ggsave('read_classification_split_by_primer_hit.png', width=8, height=8)

all %>% 
  ggplot(aes(x=code, y=read_length, fill=primer_hit)) +
  geom_violin(scale='count') +
  scale_y_log10()
ggsave('read_length_split_by_primer_hit.png', width=8, height=8)


# Graph hits vs hit length
hits <- df %>% 
  group_by(sseqid) %>% 
  summarise(hits = n(), slen = mean(slen), hit_length = sum(length), sample_length = slen - hit_length, avg_diff = mean(diff, na.rm=T)) %>% 
  left_join(classification) %>% 
  print()

hits %>% 
  filter(slen < 30000) %>% 
  ggplot(aes(x=hits, y=slen, color=code)) +
  geom_point(alpha=0.2) +
  facet_wrap(vars(code)) +
  scale_x_log10()
ggsave('hits_vs_read_length.png', width=15, height=8)


hits %>% 
  filter(slen < 40000) %>% 
  ggplot(aes(x=slen, y=hits, color=code)) +
  geom_point() +
  facet_wrap(vars(code))
ggsave('read_length_vs_hits.png', width=15, height=8)

# Graph gaps

hits %>% 
  filter(avg_diff < 100) %>% 
  ggplot(aes(x=avg_diff)) +
  geom_histogram() +
  scale_y_log10()
ggsave('avg_gap_100.png', width=10, height=10)

hits %>% 
  filter(avg_diff < 2000) %>% 
  ggplot(aes(x=avg_diff)) +
  geom_histogram() +
  scale_y_log10()
ggsave('avg_gap_2000.png', width=10, height=10)


hits %>% 
  group_by(avg_diff) %>% 
  summarise(count = n()) %>% 
  print(n=500)

gaps <- df %>% 
  filter(! is.na(diff))

gaps %>% 
  filter(diff <= 100) %>% 
  ggplot(aes(x=diff)) +
  geom_bar() +
  scale_y_log10()
ggsave('gap_counts.png', width=20, height=10)
