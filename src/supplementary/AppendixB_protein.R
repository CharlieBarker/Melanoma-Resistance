# Analysis of IOVS mouse lens data (Supplemental Table S01):
# Khan, Shahid Y., et al. "Proteome Profiling of Developing Murine Lens Through Mass Spectrometry."
# Investigative Ophthalmology & Visual Science 59.1 (2018): 100-107.

# load libraries (this gets ggplot2 and related libraries)
library(tidyverse)
library(reshape2)
# these are from Bioconductor
library(limma) 
library(edgeR) 
library(sva)
library(psych)
source("./src/functions/default_variables.R")


#find correct environment 
packLib="/usr/lib/R"
if (file.exists(packLib)) {
  setwd(dir = "~/Desktop/Melanoma_Resistance/")
}else {
  setwd(dir = "~/Desktop/Melanoma_Resistance/")
}

design<-read.csv(file = "./data/proteomic/design.csv",header = T)

list_of_inputs<-list(#signalome_view=data.frame(read.csv("input_data/signalome_by_sample.csv")),
  phospho=data.frame(read.csv("./data/proteomic/processed/combat_peptide.csv"))
  ,protein=data.frame(read.csv("./data/input_data/proteins.csv"))
  ,mRNA=data.frame(read.csv("./data/input_data/rna_expression.csv"))
)
colnames(list_of_inputs$phospho)<-design$total[match(colnames(list_of_inputs$phospho), design$code)]
# function computes CVs per sample
make_CVs <- function(df, design) {
  list_dfs<-list()
  exp_conditions<-colnames(df)
  for (variable in names(table(exp_conditions)[table(exp_conditions) != 1])) {
    cond_df<-data.frame(df[,exp_conditions==variable])
    cond_df$ave <- rowMeans(cond_df)
    cond_df$sd <- apply(data.frame(cond_df[,colnames(cond_df)!="ave"]), 1, sd)
    cond_df$cv <- 100 * cond_df$sd / cond_df$ave
    list_dfs[[variable]]<-cond_df
  }
  # save results in 3 data frames and put into a list
  ave_df <- data.frame(lapply(list_dfs, function(x){x$ave}))
  sd_df <- data.frame(lapply(list_dfs, function(x){x$sd}))
  cv_df <- data.frame(lapply(list_dfs, function(x){x$cv}))
  
  cv_molt<-reshape2::melt(cv_df)
  colnames(cv_molt)<-c("Sample", "CV")
  return(cv_molt)
}
apply_cv_on_data <-function(counts, design){
  # Use X column to make row names
  rownames(counts) <- counts$X
  counts <- counts[, -1]  # Remove the 'X' column
  # Remove the number after the full stop in the column names
  colnames(counts) <- gsub("\\.\\d+", "", colnames(counts))
  
  # Extract drug names and genetic elements from column names
  col_names <- colnames(counts)
  drug_names <- sub("^(.*?)__.*", "\\1", col_names)
  gene_names <- sub(".*?__(.*)", "\\1", col_names)
  # Replace drug names with the names in replacement vector
  new_drug_names <- names(replacement_Vec)[match(drug_names, unname(replacement_Vec))]
  
  # Combine drug names and genetic elements to form new column names
  snames <- paste(new_drug_names, gene_names, sep = ".")
  
  # Assign new column names
  colnames(counts) <- snames
  
  # get CVs and averages
  list_protein <- make_CVs(counts, design)
  return(list_protein)
}

mRNA_csv<-apply_cv_on_data(list_of_inputs$mRNA, design)
protein_csv<-apply_cv_on_data(list_of_inputs$protein, design)
phospho_csv<-apply_cv_on_data(list_of_inputs$phospho, design)


#PLOTS

colour_scheme<-c("#F94144", "#F3722C", "#F8961E", "#F9844A", "#F9C74F", "#90BE6D", "#43AA8B")
sample_colour_scheme<-c("#001219","055363", "#005F73", "#0A9396", "#94D2BD",
                        "#E9D8A6", "#EE9B00", "#CA6702", "#BB3E03",
                        "#AE2012", "#9B2226", "#C61E24", "#5727A1")
fontSize = 10

# Change box plot colors by groups
cv1<-ggplot(mRNA_csv, aes(y=CV, x=Sample, fill=Sample))+
  geom_boxplot()+cowplot::theme_cowplot(font_size = fontSize) +
  ylim(0, 100) + ggtitle("A.   Coefficient Of Variation for each mRNA transcript per sample (CVs)") + 
  labs(y = "CV (%)", x = "Samples") + 
  theme(axis.text.x = element_blank())+
  scale_fill_manual(values=sample_colour_scheme)
cv2<-ggplot(protein_csv, aes(y=CV, x=Sample, fill=Sample))+
  geom_boxplot()+cowplot::theme_cowplot(font_size = fontSize) +
  ylim(0, 100) + ggtitle("B.   Coefficient Of Variation for each Protein per sample (CVs)") + 
  labs(y = "CV (%)", x = "Samples") + 
  theme(axis.text.x = element_blank())+
  scale_fill_manual(values=sample_colour_scheme)
cv3<-ggplot(phospho_csv, aes(y=CV, x=Sample, fill=Sample))+
  geom_boxplot()+cowplot::theme_cowplot(font_size = fontSize) +
  ylim(0, 100) + ggtitle("C.   Coefficient Of Variation for each Phospho-peptide per sample (CVs)") + 
  labs(y = "CV (%)", x = "Samples") + 
  theme(axis.text.x = element_blank())+
  scale_fill_manual(values=sample_colour_scheme)

#WRITE PDFS
pdf(file = "./paper/Supplementary_plots/COV_plots.pdf",width=8.27,height=10)
ggpubr::ggarrange(cv1, cv2, cv3, nrow = 3) #CVs 
dev.off()

