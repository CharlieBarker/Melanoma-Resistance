setwd("~/Desktop/Melanoma_Resistance/")

library(ggplot2)
library(reshape2)
library(stringr)
library(purrr)
library(viridis)
library(hrbrthemes)

all_abundace<-list(
  protein=read.csv(file = "./data/input_data/proteins.csv"),
  rna=read.csv(file = "./data/input_data/rna_expression.csv")
)

proteins_of_interest<-list(
  CD44="P16070"
)

# Define the drugs and their corresponding colors
drug_colors <- c(
  "Untreated" = "#A2AEBB",
  "Vemurafenib" = "#FFBA08",
  "Trametinib" = "#D00000",
  "Combination" = "#3F88C5"
)
proteins_df<-melt(all_abundace$protein[all_abundace$protein$X %in% unlist(unname(proteins_of_interest)),])
rna_df<-melt(all_abundace$rna[all_abundace$rna$X %in% names(proteins_of_interest),])
replacement_Vec<-c("Untreated","Vermurafenib_1uM","Trametinib_10nM","vemurafenib_and_trametinib")
names(replacement_Vec)<- c("Untreated", "Vemurafenib", "Trametinib", "Combination")

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

proteins_df<-split_exp_conditions(proteins_df)
rna_df<-split_exp_conditions(rna_df)

proteins_df<-rna_df#[rna_df$ko=="WT",]


# Step 1: Call the pdf command to start the plot
pdf(file = "~/Desktop/Melanoma_Resistance/results/vis/jnk.pdf",   # The directory you want to save the file in
    width = 6, # The width of the plot in inches
    height = 6) # The height of the plot in inches

# Step 2: Create the plot with R code
ggplot(proteins_df) +
  geom_bar(aes(x = drug, y = mean, fill = drug), stat = "identity", alpha = 0.5) +
  geom_errorbar(aes(x = drug, ymin = mean - ic, ymax = mean + ic), width = 0.4, colour = "orange", alpha = 0.9, size = 1.5) +
  cowplot::theme_cowplot() + 
  facet_wrap(~X+ko) + 
  scale_fill_manual(values = drug_colors) +
  labs(y = "RNA abundance") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  geom_hline(yintercept = 0, linetype = 'dotted', col = 'red') +
  theme(axis.text.x = element_text(size = 12, hjust = 1),
        strip.text = element_text(size = 8),
        strip.background = element_rect(fill = "lightgray", color = "gray", linewidth = 0.5),
        panel.border = element_rect(color = "gray", linewidth = 0.5))
# Step 3: Run dev.off() to create the file!
dev.off()


