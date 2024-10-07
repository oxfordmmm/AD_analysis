#!/usr/bin/env Rscript
library(tidyverse)

args <- commandArgs(trailingOnly = TRUE)
read_stats_file <- args[1]
alignments_file <- args[2]

read_stats <- read_csv(read_stats_file) %>% 
  rename(query = id) %>% 
  mutate(primer_density = if_else(has_primer, 100 * 31 * primer_count / read_length, 0))

alns <- read_csv(alignments_file) %>% 
  mutate(split = case_when(
    str_detect(query, 'split') ~ as.numeric(str_split_i(query, '_split', 2)),
    T ~ 0
  )) %>% 
  mutate(query = str_split_i(query, '_split', 1)) %>% 
  group_by(query) %>% 
  summarise(ref =names(which.max(table(ref))), num_alns = n(), align_sum = sum(align_len), align_mean = mean(align_len), align_median = median(align_len), aligned_pseudoreads = n_distinct(split)) %>% 
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
  arrange(desc(num_alns), query) %>% 
  print()


df <- read_stats %>% 
  left_join(alns) %>% 
  mutate(ref = replace_na(ref, 'none')) %>% 
  replace(is.na(.), 0) %>% 
  arrange(desc(align_sum))
df
df %>% write_csv('read_align_stats.csv')


qc_summary <- df %>% 
  summarise(
    have_primer = sum(has_primer),
    have_primer_and_map = sum(has_primer & ref != 'none')
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



df %>% 
  filter(pseudoread_count > 0) %>%
  mutate(pseudoread_count = pmin(pseudoread_count, 40)) %>% 
  ggplot(aes(x = pseudoread_count, fill = ref)) +
  geom_bar()
ggsave('pseudoread_count_plot.png', width=10, height=6)

pseudo_circle_stats <- function(df) {
  df %>% 
    filter(pseudoread_count > 2) %>% 
    summarise(
      pseudo_circles = n(),
      pseudo_circle_size_Q1 = quantile(pseudoread_median_length, 0.25),
      pseudo_circle_size_Q2 = quantile(pseudoread_median_length, 0.5),
      pseudo_circle_size_Q3 = quantile(pseudoread_median_length, 0.75),
      pseudo_circle_size_Q95 = quantile(pseudoread_median_length, 0.95),
      pseudo_circle_size_mean = mean(pseudoread_median_length),
    ) %>% 
    return()
}

pseudo_circles <- df %>% 
  pseudo_circle_stats() %>% 
  print()

mapped_pseudo_circles <- df %>% 
  filter(ref != 'none') %>% 
  pseudo_circle_stats() %>% 
  rename_with(~paste("mapped_", ., sep = ""), everything()) %>% 
  print()

align_circles <- df %>% 
  filter(num_alns > 2) %>% 
  summarise(
    align_circles = n(),
    align_circle_size_Q1 = quantile(align_median, 0.25),
    align_circle_size_Q2 = quantile(align_median, 0.5),
    align_circle_size_Q3 = quantile(align_median, 0.75),
    align_circle_size_Q95 = quantile(align_median, 0.95),
    align_circle_size_mean = mean(align_median),
  ) %>% 
  print()

df %>% 
  filter(pseudoread_count > 0) %>%
  ggplot(aes(x = pseudoread_median_length, fill = ref)) +
  geom_histogram() +
  scale_x_log10()
ggsave('pseudoread_median_length_plot.png', width=10, height=6)

df %>% 
  filter(num_alns > 0) %>%
  mutate(num_alns = pmin(num_alns, 10)) %>% 
  ggplot(aes(x = num_alns, fill = ref)) +
  geom_bar()
ggsave('align_counts_plot.png', width=10, height=6)

df %>% 
  filter(num_alns > 0) %>%
  ggplot(aes(x = align_median, fill = ref)) +
  geom_histogram() +
  scale_x_log10()
ggsave('align_median_lengths_plot.png', width=10, height=6)

all_summary_stats <- bind_cols(
  reads_summary,
  qc_summary,
  bases_summary,
  read_length_quartiles,
  read_length_quartiles_mapped,
  pseudo_circles,
  mapped_pseudo_circles,
  align_circles
) %>% 
  write_csv('summary_stats.csv') %>% 
  pivot_longer(names_to = 'stat', cols=everything()) %>% 
  print(n=40)


## example circles

best_circles <- df %>% 
  top_n(10) %>% 
  pull('query')

read_csv(alignments_file) %>% 
  mutate(split = case_when(
    str_detect(query, 'split') ~ as.numeric(str_split_i(query, '_split', 2)),
    T ~ 0
  )) %>% 
  mutate(query_base = str_split_i(query, '_split', 1)) %>%
  filter(query_base %in% best_circles) %>% 
  mutate(query_base = factor(query_base, levels = best_circles)) %>% 
  arrange(query_base, split, query_start) %>% 
  write_csv("best_circles.csv") %>% 
  print()



##### Old grouping method

# levels = c('none', 'one', 'some', 'many')
# df_levels <- df %>% 
#   mutate(primer_level = case_when(
#     primer_count == 0 ~ 'none',
#     primer_count == 1 ~ 'one',
#     primer_count <= 5 ~ 'some',
#     primer_count > 5 ~ 'many'
#   )) %>% 
#   mutate(pseudoreads_level = case_when(
#     pseudoread_count == 0 ~ 'none',
#     pseudoread_count == 1 ~ 'one',
#     pseudoread_count <= 5 ~ 'some',
#     pseudoread_count > 5 ~ 'many'
#   )) %>% 
#   mutate(alns_level = case_when(
#     num_alns == 0 ~ 'none',
#     num_alns == 1 ~ 'one',
#     num_alns <= 5 ~ 'some',
#     num_alns > 5 ~ 'many'
#   )) %>% 
#   mutate(
#     primer_level = factor(primer_level, levels = levels),
#     pseudoreads_level = factor(pseudoreads_level, levels = levels),
#     alns_level = factor(alns_level, levels = levels)
#   )
# 
# levels2 <- c('no_primer_no_map', 'no_primer_one_map', 'no_primer_multi_map',
#              'primer_no_maps', 'primer_one_map', 'primer_2_to_5_maps', 'primer_many_maps')
# df_levels2 <- df_levels %>% 
#   mutate(read_type = case_when(
#     primer_level == 'none' & alns_level == 'none' ~ 'no_primer_no_map',
#     primer_level == 'none' & alns_level == 'one' ~ 'no_primer_one_map',
#     primer_level == 'none' ~ 'no_primer_multi_map',
#     alns_level == 'none' ~ 'primer_no_maps',
#     alns_level == 'one' ~ 'primer_one_map',
#     alns_level == 'some' ~ 'primer_2_to_5_maps',
#     alns_level == 'many' ~ 'primer_many_maps',
#   )) %>% 
#   mutate(read_type = factor(read_type, levels = levels2)) %>% 
#   group_by(read_type) %>% 
#   summarise(count = n(), align_avg = mean(align_sum), align_sum = sum(align_sum),
#             mean_primer_density = mean(primer_density)) %>% 
#   print() %>% 
#   write_csv('alignment_summary.csv')
# 
# df_grouped <- df_levels %>% 
#   group_by(pseudoreads_level, alns_level) %>% 
#   summarise(count = n(), align_avg = mean(align_sum), align_sum = sum(align_sum),
#             mean_primer_density = mean(primer_density)) %>% 
#   mutate(desc = sprintf("n:%s, aln_avg:%s, aln_sum:%s, densiy:%s", count, align_avg, align_sum, mean_primer_density))
# 
# make_table <- function(df, var){
#   var <- ensym(var)
#   
#   df %>% 
#     select(pseudoreads_level, alns_level, !!var) %>% 
#     pivot_wider(names_from = alns_level, values_from = !!var) %>% 
#     return()
# }  
#   
# make_table(df_grouped, count) %>% 
#   write_csv('aligns_by_segments_count.csv')
#   
# make_table(df_grouped, align_avg) %>% 
#   write_csv('aligns_by_segments_align_avg.csv')
# 
# make_table(df_grouped, align_sum) %>% 
#   write_csv('aligns_by_segments_align_sums.csv')
# 
# make_table(df_grouped, mean_primer_density) %>% 
#   write_csv('aligns_by_segments_mean_primer_density.csv')
