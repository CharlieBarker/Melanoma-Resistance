
library(data.table)
library(MOFA2)
library(dplyr)
library(purrr)
library(stringr)
library(ggplot2)
library(ComplexHeatmap)

setwd("~/phd/MelanomaProject/mofa/")

TNKS1BP1<-"Q9C0C2"
ATM<-"Q13315"
TP53BP<-"Q12888"  

phospho=data.frame(read.csv("input_data/phosphosites.csv"))
together_plot<-reshape2::melt(phospho[grep(phospho$X, pattern = "TP53BP"),])
#together_plot<-reshape2::dcast(together_plot, variable ~ X)
together_plot$Drug = unlist(map(str_split(together_plot$variable, pattern = "__"), 1))
together_plot$Genetic = gsub("\\..*","",unlist(map(str_split(together_plot$variable, pattern = "__"), 2)))
# Change point shapes, colors and sizes

# Define the drugs and their corresponding colors
drug_colors <- c(
  "Untreated" = "#A2AEBB",
  "Vemurafenib" = "#FFBA08",
  "Trametinib" = "#D00000",
  "Combination" = "#3F88C5"
)

library(dplyr)

# Define a mapping of old drug names to new names
drug_mapping <- c(
  "Untreated" = "Untreated",
  "Vermurafenib_1uM" = "Vemurafenib",
  "Trametinib_10nM" = "Trametinib",
  "vemurafenib.trametinib" = "Combination"
)

# Replace the drug names in the dataframe
together_plot <- together_plot %>%
  mutate(Drug = case_when(
    Drug %in% names(drug_mapping) ~ as.character(drug_mapping[Drug]),
    TRUE ~ as.character(Drug)
  ))

# Create the ggplot
# Save the modified ggplot as a variable
your_ggplot <- ggplot(together_plot, aes(x = Drug, y = value, fill = Drug)) +
  geom_boxplot() +
  facet_grid(Genetic ~ X, scales = "free") +
  cowplot::theme_cowplot(font_size = 25) +
  geom_jitter(color = "black", size = 0.4, alpha = 0.9) +
  scale_fill_manual(values = drug_colors) +
  theme(axis.text.x = element_text(size = 12, angle = 45, hjust = 1),
        strip.text = element_text(size = 8),
        strip.background = element_rect(fill = "lightgray", color = "gray", linewidth = 0.5),
        panel.border = element_rect(color = "gray", linewidth = 0.5))


ggsave("~/Desktop/Melanoma_Resistance/results/vis/Factor1/TP53BP_phosphosite.pdf", plot = your_ggplot, width = 8, height = 10)
