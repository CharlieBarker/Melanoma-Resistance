library(dplyr)
library(ggplot2)
library(readr)
library(ggpubr)

setwd("~/Desktop/Melanoma_Resistance/")

# --- Load data ---
mRNA_files <- list.files(path = "./results/lfc/mRNA", full.names = TRUE, recursive = TRUE)
protein_files <- list.files(path = "./results/lfc/protein", full.names = TRUE, recursive = TRUE)

# Helper to clean file names
clean_name <- function(x) gsub("_lfc\\.csv", "", gsub("(.*\\/)", "", x))

# Read and combine
mRNA_list <- lapply(mRNA_files, read.csv)
names(mRNA_list) <- clean_name(mRNA_files)
complete_mRNA <- bind_rows(mRNA_list, .id = "condition")

protein_list <- lapply(protein_files, read.csv)
names(protein_list) <- clean_name(protein_files)
complete_protein <- bind_rows(protein_list, .id = "condition")

# --- Prepare pairings: match WT vs ARID1A for each drug ---
# Expected condition naming: e.g. "vemurafenib_lfc" and "vemurafenib_arid1a_lfc"
# Adjust if your names differ slightly
pair_conditions <- tibble(
  drug = c("vemurafenib", "trametinib", "combination"),
  wt_pattern = c("vemurafenib", "trametinib", "combination"),
  ko_pattern = c("vemurafenib_arid1a", "trametinib_arid1a", "combination_arid1a")
)

# --- Function to compute correlations ---
compute_drug_cor <- function(df, drug, wt_pattern, ko_pattern, datatype = "mRNA") {
  df_wt <- df %>% filter(condition == wt_pattern) %>% select(Gene, logFC) %>% rename(wt_logFC = logFC)
  df_ko <- df %>% filter(condition == ko_pattern) %>% select(Gene, logFC) %>% rename(ko_logFC = logFC)
  merged <- inner_join(df_wt, df_ko, by = "Gene")
  
  cor_val <- cor(merged$wt_logFC, merged$ko_logFC, method = "pearson", use = "complete.obs")
  p_val <- cor.test(merged$wt_logFC, merged$ko_logFC)$p.value
  
  merged$drug <- drug
  merged$datatype <- datatype
  
  list(data = merged, stats = tibble(drug, datatype, correlation = cor_val, p_value = p_val))
}

# --- Compute for mRNA and protein ---
results_list <- list()

for (i in 1:nrow(pair_conditions)) {
  results_list[[paste0("mRNA_", pair_conditions$drug[i])]] <-
    compute_drug_cor(complete_mRNA, pair_conditions$drug[i],
                     pair_conditions$wt_pattern[i],
                     pair_conditions$ko_pattern[i],
                     "mRNA")
  
  results_list[[paste0("protein_", pair_conditions$drug[i])]] <-
    compute_drug_cor(complete_protein, pair_conditions$drug[i],
                     pair_conditions$wt_pattern[i],
                     pair_conditions$ko_pattern[i],
                     "protein")
}

# Combine data and stats
merged_data <- bind_rows(lapply(results_list, function(x) x$data))
cor_summary <- bind_rows(lapply(results_list, function(x) x$stats))
print(cor_summary)

# Replace drug and datatype names with prettier labels
merged_data <- merged_data %>%
  mutate(
    drug = factor(drug,
                  levels = c("vemurafenib", "trametinib", "combination"),
                  labels = c("Vemurafenib", "Trametinib", "Combination")),
    datatype = factor(datatype,
                      levels = c("mRNA", "protein"),
                      labels = c("Transcriptomics", "Proteomics"))
  )

pdf(file = 'paper/Supplementary_plots/arid1a_vs_wt_drug_response_supp.pdf', width = 14, height = 6)
# --- Plot correlations ---
ggplot(merged_data, aes(x = wt_logFC, y = ko_logFC)) +
  geom_point(alpha = 0.5) +
  geom_smooth(method = "lm", se = FALSE, color = "black") +
  facet_grid(datatype ~ drug, scales = 'free') +
  cowplot::theme_cowplot() +
  theme(
    plot.title = element_text(size=18),
    axis.text.x = element_text(angle = 70, hjust = 1, size = rel(1)),
  ) +
  xlab("") +
  stat_cor(aes(label = paste(..r.label.., ..p.label.., sep = "~`,`~")),
           method = "pearson", label.x.npc = "left", label.y.npc = "top") +
  labs(
    title = "Correlation of ARID1A-WT vs ARID1A-KO responses per drug",
    x = "ARID1A-WT log2 Fold Change",
    y = "ARID1A-KO log2 Fold Change"
  )   + grids(linetype = "dashed")

dev.off()
