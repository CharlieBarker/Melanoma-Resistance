library(ggplot2)
library(dplyr)
library(ggrepel)
# Load the volcanoPlus package
library(volcanoPlus)
library(ghibli)
library(ggpubr)

# 
# ██████╗░██╗░░██╗░█████╗░░██████╗██╗░░██╗
# ██╔══██╗██║░░██║██╔══██╗██╔════╝╚██╗██╔╝
# ██████╔╝███████║██║░░██║╚█████╗░░╚███╔╝░
# ██╔═══╝░██╔══██║██║░░██║░╚═══██╗░██╔██╗░
# ██║░░░░░██║░░██║╚█████╔╝██████╔╝██╔╝╚██╗
# ╚═╝░░░░░╚═╝░░╚═╝░╚════╝░╚═════╝░╚═╝░░╚═╝
# 
# Version 0.5.2
# Copyright (C) 2024 Alessandro Lussana
# Licence Apache 2.0
# 
# Command: /home/charlie/miniconda3/envs/phuego/bin/phosx -h
# 
# usage: phosx [-h] [-p PSSM] [-q PSSM_QUANTILES] [-n N_PERMUTATIONS] [-k N_TOP_KINASES] [-m MIN_N_HITS] [-c N_PROC] [--plot-figures] [-d OUTPUT_DIR] [-o OUTPUT_PATH] [-v] seqrnk
# 
# Kinase activity inference from phosphosproteomics data based on substrate sequence specificity
# 
# positional arguments:
#   seqrnk                Path to the seqrnk file.
# 
# options:
#   -h, --help            show this help message and exit
# -p PSSM, --pssm PSSM  Path to the h5 file storing custom PSSMs; defaults to built-in PSSMs
# -q PSSM_QUANTILES, --pssm-quantiles PSSM_QUANTILES
# Path to the h5 file storing custom PSSM score quantile distributions under the key 'pssm_scores'; defaults to built-in PSSM scores quantiles
# -n N_PERMUTATIONS, --n-permutations N_PERMUTATIONS
# Number of random permutations; default: 1000
# -k N_TOP_KINASES, --n-top-kinases N_TOP_KINASES
# Number of top-scoring kinases potentially associatiated to a given phosphosite; default: 8
# -m MIN_N_HITS, --min-n-hits MIN_N_HITS
# Minimum number of phosphosites associated with a kinase for the kinase to be considered in the analysis; default: 4
# -c N_PROC, --n-proc N_PROC
# Number of cores used for multithreading; default: 1
# --plot-figures        Save figures in pdf format; see also --output_dir
# -d OUTPUT_DIR, --output-dir OUTPUT_DIR
# Output files directory; only relevant if used with --plot_figures; defaults to 'phosx_output/'
# -o OUTPUT_PATH, --output-path OUTPUT_PATH
# Main output table; if not specified it will be printed in STDOUT
# -v, --version         Print package version and exit

#   Command: /home/charlie/miniconda3/envs/phuego/bin/phosx -m 20 WT_KO_comb_vs_untreated.seqrnk


# Read data
wt <- read.delim("results/lfc_phosX/wt_kinase_activities.tmp")
ko <- read.delim("results/lfc_phosX/kinase_activities.tmp")


df <- data.frame(
  logFC = wt$KS,
  adj.P.Val = wt$p.value,
  Gene = wt$X
)

# Calculate -log10 adjusted p-values
df$negLogP <- -log10(df$adj.P.Val)

# Define significance thresholds
logFC_cutoff <- 0.2
pval_cutoff <- 0.05

# Assign colors
df$color <- "NotSig"
df$color[df$logFC > logFC_cutoff & df$adj.P.Val < pval_cutoff] <- "Up"
df$color[df$logFC < -logFC_cutoff & df$adj.P.Val < pval_cutoff] <- "Down"


# Flag kinases of interest
genes_of_interest <- c("MAPK1", "MAP2K1", "BRAF")
df$highlight <- ifelse(df$Gene %in% genes_of_interest, "Yes", "No")

ggplot(df, aes(x = logFC, y = negLogP, color = color)) +
  geom_point(size = 3, alpha = 0.8) +
  # Highlight points of interest with larger shape
  geom_point(data = subset(df, highlight == "Yes"),
             aes(x = logFC, y = negLogP), 
             shape = 21, fill = "blue", color = "black", size = 4, stroke = 1.2) +
  # Label genes of interest
  ggrepel::geom_text_repel(data = subset(df, highlight == "Yes"),
                           aes(label = Gene), size = 4, color = "black", min.segment.length = 3) +
  scale_color_manual(values = c("Down" = "blue", "Up" = "red", "NotSig" = "grey")) +
  geom_vline(xintercept = c(-logFC_cutoff, logFC_cutoff), linetype = "dashed") +
  geom_hline(yintercept = -log10(pval_cutoff), linetype = "dashed") +
  cowplot::theme_cowplot() +
  theme(
    plot.title = element_text(size=18),
    axis.text.x = element_text(angle = 70, hjust = 1, size = rel(1)),
  ) +
  xlab("") +
  grids(linetype = "dashed")+
  labs(
    x = "Activity Score",
    y = "P value",
    title = "A375 Kinase activity after treatment with combination therapy"
  )

# Merge by kinase name (X)
merged <- wt %>%
  dplyr::select(X, Activity.Score) %>%
  rename(WT_Activity = Activity.Score) %>%
  inner_join(
    ko %>% dplyr::select(X, Activity.Score) %>% rename(KO_Activity = Activity.Score),
    by = "X"
  )

# Calculate distance from diagonal
merged <- merged %>%
  mutate(diff = abs(WT_Activity - KO_Activity))

# Keep top 10 most different for labeling
merged$label <- ''
merged$label[grep(merged$X, pattern = 'MAPK')] <- merged$X[grep(merged$X, pattern = 'MAPK')]

# Scatter plot with labels
ggplot(merged, aes(x = WT_Activity, y = KO_Activity)) +
  geom_point(color = "steelblue", size = 3, alpha = 0.7) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
  geom_text_repel(
    data = merged,
    aes(label = label),
    size = 4,
    max.overlaps = 20,
    box.padding = 0.3,
    point.padding = 0.2
  ) +
  theme_minimal(base_size = 14) +
  labs(
    title = "Kinase activity comparison",
    x = "WT Activity Score",
    y = "ARID1A KO Activity Score"
  )
