home_dir<-"~/Desktop/Melanoma_Resistance/"
setwd(home_dir)

source("./src/functions/default_variables.R")

library(pathview)
library(MOFA2)
library(tidyr)
library(ggplot2)
library(dplyr)
library(ggpubr)


###enrichment of the weights #### 
MOFAobject.trained<-load_model(file = "./results/mofa/mofa_object.hdf5")

weights <- get_weights(MOFAobject.trained, 
                       views = "all", 
                       as.data.frame = TRUE 
)

list_of_inputs<-list(#signalome_view=data.frame(read.csv("input_data/signalome_by_sample.csv")),
  phospho=data.frame(read.csv("./data/input_data/phosphosites.csv"))
  ,protein=data.frame(read.csv("./data/input_data/proteins.csv"))
  ,mRNA=data.frame(read.csv("./data/input_data/rna_expression.csv"))
)

# Define the drugs and their corresponding colors
drug_colors <- c(
  "Untreated" = "#A2AEBB",
  "Vemurafenib" = "#FFBA08",
  "Trametinib" = "#D00000",
  "Combination" = "#3F88C5"
)

RTKs<-c("IGF1R", 
        "ERBB3", "EGFR", "EGF",
        "EPHA2", "EPHA7",
        "KIT",
        "NTRK2",
        "KDR",
        "FLT1",
        "FGFR2", "FGF1", "FGFR1", "FGF2",
        "CD44", 
        "INSR",
        "MMP1")
RTK_figure1 <- c("EGFR", "EGF", "KIT", "FGFR2", "FGF1", "FGFR1", "FGF2")
neg_feedback_dusp<-c("DUSP1", "DUSP14", "DUSP2",   "DUSP21",  "DUSP4",   "DUSP6")
neg_feedback_spry<-c("SPRED1", "SPRED2", "SPRY1", "SPRY2", "SPRY4")

#for interest
for_interest<-c("MMP1", "HLA.DRB1")

all_of_interest <- c(RTK_figure1, neg_feedback_dusp, neg_feedback_spry)
rtk_mrna<-reshape2::melt(list_of_inputs$mRNA[list_of_inputs$mRNA$X %in% all_of_interest,])
rtk_mrna$X <- factor(rtk_mrna$X, levels=all_of_interest)
rtk_mrna$variable<-sub("*\\.[0-9]", "", rtk_mrna$variable)
rtk_mrna <- rtk_mrna |>
  separate_wider_delim(variable, delim = "__", names = c("drug", "ko"))


rtk_mrna$drug<-names(replacement_Vec)[match(rtk_mrna$drug, unname(replacement_Vec))]
# Define the desired order of levels
desired_order <- c("Untreated", "Vemurafenib", "Trametinib", "Combination")
# Convert my_column to a factor with the specified order
rtk_mrna$drug <- factor(rtk_mrna$drug, levels = desired_order)

output_file<-"./paper/plots/phuego/factor1_rna_abundances.pdf"
pdf(# The directory you want to save the file in
  width = 14, # The width of the plot in inches
  height = 5,
  file = output_file)
rtk_mrna %>%
  filter(ko == "WT")  %>%
  ggplot( aes(x=drug, y=value, fill=drug)) +
  geom_boxplot() +
  geom_jitter(color="black", size=0.4, alpha=0.9) +
  scale_fill_manual(values = drug_colors) +
  cowplot::theme_cowplot() +
  theme(
    plot.title = element_text(size=11),
    axis.text.x = element_text(angle = 70, hjust = 1, size = rel(1)),
  ) +
  xlab("") +
  facet_wrap(~X, scales = "free_y", nrow = 2) + 
  grids(linetype = "dashed")+
  labs(
    x = "Drug treatment",
    y = "Voom-normalised RNA expression",
  )
dev.off()
  
