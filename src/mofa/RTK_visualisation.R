home_dir<-"~/Desktop/Melanoma_Resistance/"
setwd(home_dir)

library(pathview)
library(MOFA2)
library(tidyr)
library(ggplot2)

###enrichment of the weights #### 
MOFAobject.trained<-load_model(file = "./results/mofa/mofa_object.hdf5")

weights <- get_weights(MOFAobject.trained, 
                       views = "all", 
                       as.data.frame = TRUE 
)

replacement_Vec<-c("Untreated","Vermurafenib_1uM","Trametinib_10nM","vemurafenib_and_trametinib")
names(replacement_Vec)<- c("Untreated", "Vemurafenib", "Trametinib", "Combination")

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
neg_feedback<-c("DUSP1", "DUSP2", "DUSP4")
#for interest
for_interest<-c("MMP1", "HLA.DRB1")


rtk_mrna<-reshape2::melt(list_of_inputs$mRNA[list_of_inputs$mRNA$X %in% for_interest,])
rtk_mrna$variable<-sub("*\\.[0-9]", "", rtk_mrna$variable)
rtk_mrna <- rtk_mrna |>
  separate_wider_delim(variable, delim = "__", names = c("drug", "ko"))


rtk_mrna$drug<-names(replacement_Vec)[match(rtk_mrna$drug, unname(replacement_Vec))]
# Define the desired order of levels
desired_order <- c("Untreated", "Vemurafenib", "Trametinib", "Combination")
# Convert my_column to a factor with the specified order
rtk_mrna$drug <- factor(rtk_mrna$drug, levels = desired_order)

rtk_mrna %>%
  ggplot( aes(x=ko, y=value, fill=drug)) +
  geom_boxplot() +
  geom_jitter(color="black", size=0.4, alpha=0.9) +
  scale_fill_manual(values = drug_colors) +
  theme(
    legend.position="none",
    plot.title = element_text(size=11)
  ) +
  ggtitle("mRNA abundances of central RTKs") +
  xlab("") + cowplot::theme_cowplot() +
  facet_wrap(~X, scales = "free_y")
