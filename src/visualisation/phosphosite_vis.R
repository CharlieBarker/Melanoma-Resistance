
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
  "Vermurafenib_1uM" = "#FFBA08",
  "Trametinib_10nM" = "#D00000",
  "vemurafenib.trametinib" = "#3F88C5"
)


# Create the ggplot
ggplot(together_plot, aes(x = Drug, y = value, fill = Drug)) +
  geom_boxplot() +
  facet_wrap(Genetic ~ X, scales = "free") +
  cowplot::theme_cowplot(font_size = 25) +
  geom_jitter(color = "black", size = 0.4, alpha = 0.9) +
  scale_fill_manual(values = drug_colors)

