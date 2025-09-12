
library(data.table)
library(dplyr)
library(purrr)
library(stringr)
library(ggplot2)
library(ComplexHeatmap)

setwd("~/Desktop/Melanoma_Resistance//")


# Define the drugs and their corresponding colors
drug_colors <- c(
  "Untreated" = "#A2AEBB",
  "Vemurafenib" = "#FFBA08",
  "Trametinib" = "#D00000",
  "Combination" = "#3F88C5"
)
replacement_Vec<-c("Untreated","Vermurafenib_1uM","Trametinib_10nM","vemurafenib.trametinib")
names(replacement_Vec)<- c("Untreated", "Vemurafenib", "Trametinib", "Combination")


TNKS1BP1<-"Q9C0C2"
ATM<-"Q13315"
TP53BP<-"Q12888"  

split_exp_conditions<-function(df){
  df$drug<-unlist(map(str_split(df$variable, pattern = "__"),1))
  df$ko<-unlist(map(str_split(df$variable, pattern = "__"),2))
  df$ko<-gsub("\\..*","",df$ko)
  df$variable<-NULL
  df$drug<-names(replacement_Vec)[match(df$drug, unname(replacement_Vec))]
  # Define the desired order of levels
  desired_order <- c("Untreated", "Vemurafenib", "Trametinib", "Combination")
  # Convert my_column to a factor with the specified order
  df$drug <- factor(df$drug, levels = desired_order)
  
  df <- df %>%
    group_by(X, drug, ko) %>%
    summarise( 
      n=n(),
      mean=mean(value),
      sd=sd(value)
    ) %>%
    mutate( se=sd/sqrt(n))  %>%
    mutate( ic=se * qt((1-0.05)/2 + .5, n-1))
  
  return(df)
}

phospho=data.frame(read.csv("data/input_data/unfiltered_phosphosites.csv"))
phospho_lm=data.frame(read.csv("results/lfc/phospho/combination_lfc.csv"))
library(ggplot2)
library(ggrepel)

# Your protein list

# Add columns for plotting
phospho_lm$negLog10FDR <- -log10(phospho_lm$adj.P.Val)
phospho_lm$Protein <- sub(";.*", "", phospho_lm$Gene)  # extract protein before first ";"

proteins_of_interest <- unique(phospho_lm$Protein[phospho_lm$adj.P.Val < 0.001])

phospho_lm$highlight <- ifelse(phospho_lm$Protein %in% proteins_of_interest, phospho_lm$X, NA)


ggplots_list<-list()
for (protein_of_interest in c("LEO1", "RIF1", 'SIRT1', 'SRRM2', 'SSRP1', 'ULK1', 'ATM')) {
  together_plot<-reshape2::melt(phospho[grep(phospho$X, pattern = protein_of_interest),])
  phos_df<-split_exp_conditions(together_plot)
  phos_wt<-phos_df[phos_df$ko=="WT",]
  ggplots_list[[protein_of_interest]]<-phos_wt
}

phos_total<-bind_rows(ggplots_list, .id = "psite")
to_miss<-c("TNKS1BP1;S712;", "TNKS1BP1;S920;")

# Step 1: Call the pdf command to start the plot
pdf(file = "~/Desktop/Melanoma_Resistance/results/vis/Factor1/phos_vis.pdf",   # The directory you want to save the file in
    width = 12, # The width of the plot in inches
    height = 10) # The height of the plot in inches

# Volcano plot
ggplot(phospho_lm, aes(x = logFC, y = negLog10FDR)) +
  geom_point(aes(color = !is.na(highlight)), alpha = 0.7, size = 2) +
  scale_color_manual(values = c("grey60", "red")) +
  geom_text_repel(aes(label = highlight), size = 3, max.overlaps = 10) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "black") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +
  theme_minimal(base_size = 14) +
  labs(
    title = "Volcano Plot of Phosphosites",
    x = "Log2 Fold Change",
    y = "-Log10 Adjusted P-value"
  )
# Step 2: Create the plot with R code
ggplot(phos_total) +
  geom_bar(aes(x = drug, y = mean, fill = drug), stat = "identity", alpha = 0.5) +
  geom_errorbar(aes(x = drug, ymin = mean - ic, ymax = mean + ic), width = 0.4, colour = "orange", alpha = 0.9, size = 1.5) +
  cowplot::theme_cowplot() + 
  facet_wrap(~X, nrow = 3) + 
  scale_fill_manual(values = drug_colors) +
  labs(y = "Phosphosite abundance") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  geom_hline(yintercept = 0, linetype = 'dotted', col = 'red') +
  theme(axis.text.x = element_text(size = 12, hjust = 1),
        strip.text = element_text(size = 8),
        strip.background = element_rect(fill = "lightgray", color = "gray", linewidth = 0.5),
        panel.border = element_rect(color = "gray", linewidth = 0.5))



library(dplyr)

# Fisher’s method helper
fisher_method <- function(pvals) {
  pvals <- pvals[!is.na(pvals) & pvals > 0]  # drop NAs/zeros
  if (length(pvals) == 0) return(NA)
  chisq <- -2 * sum(log(pvals))
  pval <- pchisq(chisq, df = 2 * length(pvals), lower.tail = FALSE)
  return(pval)
}

# Apply per protein
protein_level <- phospho_lm %>%
  mutate(Protein = sub(";.*", "", Gene)) %>%      # strip protein name
  group_by(Protein) %>%
  summarise(
    n_sites = n(),
    fisher_p = fisher_method(P.Value),
    .groups = "drop"
  ) %>%
  mutate(fisher_FDR = p.adjust(fisher_p, method = "BH")) %>%
  arrange(fisher_FDR)
top_30_proteins<-head(protein_level$Protein, 50)

# compute theoretical curve for identical per-site p (e.g. 0.1)
p_site <- 0.01
max_k <- max(protein_level$n_sites, na.rm = TRUE)

theo <- data.frame(n_sites = 1:max_k)
theo$fisher_p <- sapply(theo$n_sites, function(k) {
  chisq <- -2 * k * log(p_site)
  pchisq(chisq, df = 2*k, lower.tail = FALSE)
})
theo$minuslog10_p <- -log10(theo$fisher_p)


# Join protein_level with theo by n_sites
protein_level_labeled <- protein_level %>%
  left_join(theo, by = "n_sites") %>%
  # Add a label column: TRUE if protein fisher_p < theo fisher_p, FALSE otherwise
  mutate(label = ifelse(fisher_p.x < fisher_p.y, TRUE, FALSE)) %>%
  # Clean up column names
  dplyr::select(Protein, n_sites, fisher_p = fisher_p.x, fisher_FDR, theo_fisher_p = fisher_p.y, label)
protein_level_labeled$highlight <- ''
protein_level_labeled$highlight[protein_level_labeled$label]<-protein_level_labeled$Protein[protein_level_labeled$label]

protein_level_labeled %>%
  ggplot(aes(x = n_sites, y = -log10(fisher_p))) +  # use raw fisher_p
  geom_point(aes(color = label), size = 2, alpha = 0.7) +
  scale_color_manual(values = c("grey60", "red")) +
  geom_text_repel(aes(label = highlight), size = 5, max.overlaps = 10) +
  geom_line(data = theo, aes(x = n_sites, y = minuslog10_p), color = "blue", linetype = "dashed", size = 1) +
  labs(
    title = "Protein-level Phosphorylation Enrichment (Fisher’s method)",
    x = "Number of phosphosites tested",
    y = "-log10(Fisher p-value)",
    color = "Highlighted proteins"
  )  + 
  cowplot::theme_cowplot() +
  theme(
    plot.title = element_text(size = 18, face = "bold"),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 1)
  ) +
  grids(linetype = "dashed") 

dev.off()

phos_total<-phos_total[!phos_total$X %in% to_miss,]
