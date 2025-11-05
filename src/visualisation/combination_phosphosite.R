
library(data.table)
library(dplyr)
library(purrr)
library(stringr)
library(ggplot2)
library(ComplexHeatmap)
library(wesanderson)

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

phospho_lm=data.frame(read.csv("results/lfc/phospho/combination_lfc.csv"))
phospho_lm$Protein <- sub(";.*", "", phospho_lm$Gene)  # extract protein before first ";"
# Remove everything before the first semicolon, including the semicolon
psite_only <- sub(".*?;", "", phospho_lm$Gene)

# Remove trailing semicolon if present
psite_only <- sub(";$", "", psite_only)
phospho_lm$psite <- psite_only

interesting_proteins<-c("LEO1", "RIF1", 'SIRT1', 'SRRM2', 'SSRP1', 'ULK1', 'ATM', 'TP53BP1')
interesting_phospho_lm <- phospho_lm[phospho_lm$Protein %in% interesting_proteins & phospho_lm$adj.P.Val < 0.1,]

# Define mapping of proteins to categories
protein_categories <- c(
  "ATM" = "DNA Repair",
  "TP53BP1" = "DNA Repair",
  "RIF1" = "DNA Repair",
  "SIRT1" = "EMT",
  "ULK1" = "EMT",
  "LEO1" = "Stemness",
  "SRRM2" = "Spliceosome",
  "SSRP1" = "Spliceosome"
)

# Add Category column to your data
interesting_phospho_lm <- interesting_phospho_lm %>%
  mutate(Category = protein_categories[Protein])
interesting_phospho_lm$sign <- interesting_phospho_lm$logFC < 0

interesting_phospho_lm <- interesting_phospho_lm %>%
  mutate(
    Category = protein_categories[Protein],
    sign = logFC < 0,
    sig_label = case_when(
      adj.P.Val < 0.001 ~ "***",
      adj.P.Val < 0.01  ~ "**",
      adj.P.Val < 0.05  ~ "*",
      TRUE ~ ""
    ),
    # position for the stars
    y_star = ifelse(logFC >= 0, logFC + 0.1*max(abs(logFC)),  # above bar
                    logFC - 0.1*max(abs(logFC)))              # below bar
  )

# Step 1: Call the pdf command to start the plot
pdf(file = "~/Desktop/Melanoma_Resistance/paper/Figures/combination_phosphosites.pdf",   # The directory you want to save the file in
    width = 25, # The width of the plot in inches
    height = 3.5) # The height of the plot in inches
write.csv(interesting_phospho_lm, file = 'paper/all_figure_data/3C.csv')
ggplot(interesting_phospho_lm,
       aes(x = psite, y = logFC, fill=sign)) +
  geom_col(width = 0.7) +
  geom_text(aes(y = y_star, label = sig_label), vjust = 0.5, size = 5) +
  facet_grid(~ Category+ Protein, scales = "free_x", space = "free_x") + 
  geom_hline(yintercept = 0, color = "black", linetype = "dashed") +  # <-- horizontal line
  cowplot::theme_cowplot(font_size = 16) +  # increase base font size
  theme(
    plot.title = element_text(size = 20),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 14),
    axis.text.y = element_text(size = 14),
    strip.text = element_text(size = 14),  # facet labels
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 1)
  ) +
  labs(x = "Phosphosite Quantified", y = "log2 Fold Change") +
  grids(linetype = "dashed") +
  scale_fill_manual(values = wes_palette("FrenchDispatch"))

dev.off()

