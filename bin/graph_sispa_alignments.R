#!/usr/bin/env Rscript

library(tidyverse)

args <- commandArgs(trailingOnly = TRUE)
read_stats_file <- args[1]
alignments_file <- args[2]
barcode <- args[3]

read_stats <- read_csv(read_stats_file) %>% 
  rename(query = id) %>% 
  replace(is.na(.), '') %>% 
  mutate(primer_density = if_else(has_primer, 100 * 31 * primer_count / read_length, 0))

# Make this more reproducable. Also add in upper/lower case
alns <- read_csv(alignments_file) %>% 
  select(query, ref, align_len) %>% 
  mutate(ref = if_else(str_detect(ref, 'orthoreovirus'), 'orthoreovirus', ref)) %>% 
  mutate(ref = if_else(str_detect(ref, 'Influenza_A_virus'), 'Influenza_A', ref)) %>%
  mutate(ref = if_else(str_detect(ref, 'Influenza_B_virus'), 'Influenza_B', ref)) %>%
  mutate(ref = if_else(str_detect(ref, 'Rhinovirus|rhinovirus'), 'Rhinovirus', ref)) %>%
  mutate(ref = if_else(str_detect(ref, 'entero|Entero'), 'Enterovirus', ref)) %>%
  mutate(ref = if_else(str_detect(ref, 'Respiratory_syncytial_virus'), 'RSV', ref)) %>%
  mutate(ref = if_else(str_detect(ref, 'Coxsackie'), 'Coxsackievirus', ref)) %>%
  mutate(ref = if_else(str_detect(ref, 'echovirus'), 'Echovirus', ref)) %>%
  mutate(ref = if_else(str_detect(ref, 'Human_parainfluenza'), 'Human_parainfluenza', ref)) %>%
  mutate(ref = if_else(str_detect(ref, 'Severe_acute_respiratory_syndrome_coronavirus_2|Human_coronavirus_OC43'), 'Coronavirus', ref)) %>% 
  print()


df <- read_stats %>% 
  left_join(alns) %>% 
  mutate(ref = replace_na(ref, 'none')) %>% 
  replace(is.na(.), 0) %>% 
  arrange(desc(align_len))
df
df$barcode <- barcode
df %>% write_csv('read_align_stats.csv')


qc_summary <- df %>% 
  summarise(
    have_primer = sum(has_primer),
    have_primer_and_map = sum(has_primer & ref != 'none'),
    failed_trimming = sum(failed_trimming),
    primer_in_middle = sum(fail_reason == 'primer_in_middle'),
    too_many_primers = sum(fail_reason == 'too_many_primers'),
    resulting_read_too_short = sum(fail_reason == 'resulting_read_too_short'),
  ) %>% 
  print()

reads_summary <- df %>% 
  group_by(ref) %>% 
  summarise(read_count = n()) %>% 
  pivot_wider(names_from = ref, values_from = read_count, names_prefix = "reads_" ) %>% 
  mutate(reads_all = rowSums(across(everything())), reads_mapped = reads_all - reads_none) %>% 
  select(reads_all, reads_mapped, everything(), -reads_none) %>% 
  print()

bases_summary <- df %>% 
  group_by(ref) %>% 
  summarise(base_count = sum(read_length)) %>% 
  pivot_wider(names_from = ref, values_from = base_count, names_prefix = "bases_" ) %>% 
  mutate(bases_all = rowSums(across(everything())), bases_mapped = bases_all - bases_none) %>% 
  select(bases_all, bases_mapped, everything(), -bases_none) %>% 
  print()

df %>% 
  ggplot(aes(x = read_length, fill = ref)) +
  geom_histogram(bins = 50) +
  scale_x_log10()
ggsave('read_length_plot.png', width=10, height=6)

read_length_quartiles <- df %>%
  summarise(
    read_length_Q1 = quantile(read_length, 0.25),
    read_length_Q2 = quantile(read_length, 0.5), # Median (50th percentile)
    read_length_Q3 = quantile(read_length, 0.75),
    read_length_mean = mean(read_length),
  ) %>% 
  print()
read_length_quartiles_mapped <- df %>%
  filter(ref != 'none') %>% 
  summarise(
    read_length_Q1 = quantile(read_length, 0.25),
    read_length_Q2 = quantile(read_length, 0.5), # Median (50th percentile)
    read_length_Q3 = quantile(read_length, 0.75),
    read_length_mean = mean(read_length),
  ) %>% 
  rename_with(~paste("mapped_", ., sep = ""), everything()) %>% 
  print()


barcodes <- tibble(barcode = barcode) %>% print()
all_summary_stats <- bind_cols(
  barcodes,
  reads_summary,
  qc_summary,
  bases_summary,
  read_length_quartiles,
  read_length_quartiles_mapped,
) %>% 
  write_csv('summary_stats.csv')# %>% 
  #pivot_longer(names_to = 'stat', cols=everything()) %>% 
  #print(n=40)

