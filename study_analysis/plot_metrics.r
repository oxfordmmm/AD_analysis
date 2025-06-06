library(tidyverse)

# Read in data from command line argument
args <- commandArgs(trailingOnly = TRUE)
df <- read_csv(args[1])

# Recode 'pass' column
df <- df %>%
  mutate(pass = if_else(pass == '0', 'False', pass),
         `full pass` = if_else(pass == 'True' & PCs_passed == TRUE, TRUE, FALSE))

print(unique(df$pass))
print(unique(df$`full pass`))

# Create classification result columns
df <- df %>%
  mutate(
    TP = if_else(gold_standard == TRUE & `full pass` == TRUE, 1, 0),
    FP = if_else(gold_standard == FALSE & `full pass` == TRUE, 1, 0),
    TN = if_else(gold_standard == FALSE & `full pass` == FALSE, 1, 0),
    FN = if_else(gold_standard == TRUE & `full pass` == FALSE, 1, 0)
  )

# Select relevant columns
cols <- c('Run', 'barcode', 'seq_name', 'reverse_transcription_control', '2 reads pass',
          'pathogen', 'TP', 'FP', 'TN', 'FN',
          'Sample_reads_percent_of_run', 'Sample_reads_percent_of_refs',
          'Sample_reads_percent_of_type_run', 'AuG_trunc10')
df <- df %>% select(all_of(cols))

# Melt (pivot_longer) classification columns
df_melted <- df %>%
  pivot_longer(cols = c(TP, FP, TN, FN), names_to = "Test result", values_to = "Test result Value") %>%
  filter(`Test result Value` == 1) %>%
  mutate(`AuG_trunc10/Sample_reads_percent_of_refs` = Sample_reads_percent_of_refs / AuG_trunc10)

# Print summary of the ratio
summary(df_melted$`AuG_trunc10/Sample_reads_percent_of_refs`)

# Plot using ggplot2
ggplot(df_melted, aes(x = AuG_trunc10, y = Sample_reads_percent_of_refs,
                      color = `Test result`, shape = reverse_transcription_control)) +
  geom_point() +
  geom_abline(intercept = 0, slope = 0.1, color = "red") +
  geom_hline(yintercept = 0.007, linetype = "dashed", color = "grey") +
  geom_vline(xintercept = 0.003, linetype = "dashed", color = "grey") +
  scale_x_log10() +
  scale_y_log10() +
  labs(title = "Validation set, test_result includes AND ratio") +
  theme_minimal()

# Save plot
ggsave("AuG_trunc10_vs_Sample_reads_percent_of_refs_ggplot.pdf", width = 8, height = 6)
