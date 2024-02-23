
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

phospho=data.frame(read.csv("data/input_data/phosphosites.csv"))

ggplots_list<-list()
for (protein_of_interest in c("ATM", "TP53BP")) {
  together_plot<-reshape2::melt(phospho[grep(phospho$X, pattern = protein_of_interest),])
  phos_df<-split_exp_conditions(together_plot)
  phos_wt<-phos_df[phos_df$ko=="WT",]
  ggplots_list[[protein_of_interest]]<-phos_wt
}

phos_total<-bind_rows(ggplots_list, .id = "psite")
to_miss<-c("TNKS1BP1;S712;", "TNKS1BP1;S920;")
phos_total<-phos_total[!phos_total$X %in% to_miss,]

# Step 1: Call the pdf command to start the plot
pdf(file = "~/Desktop/Melanoma_Resistance/results/vis/Factor1/phos_vis.pdf",   # The directory you want to save the file in
    width = 6, # The width of the plot in inches
    height = 6) # The height of the plot in inches

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

# Step 3: Run dev.off() to create the file!
dev.off()


kinase_module<-read.csv(file = "./results/vis/Networks/Factor2_up_kinase_module.csv")
