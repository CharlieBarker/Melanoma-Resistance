home_dir<-"~/Desktop/Melanoma_Resistance/"
setwd(home_dir)

library(pathview)
library(MOFA2)
library(tidyr)
library(ggplot2)
library(dplyr)
library(ggpubr)
library(GGally)

###enrichment of the weights #### 
MOFAobject.trained<-load_model(file = "./results/mofa/mofa_object.hdf5")

weights <- get_weights(MOFAobject.trained, 
                       views = "all", 
                       as.data.frame = TRUE 
)

replacement_Vec<-c("Untreated","Vermurafenib_1uM","Trametinib_10nM","vemurafenib.trametinib")
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

# Define vector of proteins of interest
#Factor 1
prot_of_interest <- c("ERBB3", "BCAR3", "SPRED1", "ARHGEF2", "ABL2", "RRAGC")

# Initialize an empty list to store filtered data for each protein
filtered_data <- list()

# Loop over each protein of interest
for (prot in prot_of_interest) {
  # Use grep to filter the data based on the protein name
  filtered_data[[prot]] <- list_of_inputs$phospho[grep(list_of_inputs$phospho$X, pattern = prot),]
}
df <- do.call("rbind", filtered_data)

rownames(df)<-df$X
df$X<-NULL
rtk_mrna<-data.frame(t(df))

#rtk_mrna$X <- factor(rtk_mrna$X, levels=for_interest)
rtk_mrna$variable<-sub("*\\.[0-9]", "", rownames(rtk_mrna))
rtk_mrna <- rtk_mrna |>
  separate_wider_delim(variable, delim = "__", names = c("drug", "ko"))


rtk_mrna$drug<-names(replacement_Vec)[match(rtk_mrna$drug, unname(replacement_Vec))]
# Define the desired order of levels
desired_order <- c("Untreated", "Vemurafenib", "Trametinib", "Combination")
# Convert my_column to a factor with the specified order
rtk_mrna$drug <- factor(rtk_mrna$drug, levels = desired_order)


# nCOL<-ncol(rtk_mrna)-2
# ggpairs(rtk_mrna[rtk_mrna$ko=="WT",],                 # Data frame
#         columns = 1:nCOL,        # Columns
#         aes(color = drug,  # Color by group (cat. variable)
#             alpha = 0.5),
#         upper = list(continuous = "points")) + cowplot::theme_cowplot()+grids(linetype = "dashed") +  scale_fill_manual(values = drug_colors) +  scale_color_manual(values = drug_colors)
# 

molten_rtk_mrna<-reshape2::melt(rtk_mrna)
molten_rtk_mrna %>%
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
  facet_wrap(~variable, scales = "free_y", nrow = 2) + 
  grids(linetype = "dashed")+
  labs(
    x = "Drug treatment",
    y = "Psite abundance",
  )