setwd("~/Desktop/Melanoma_Resistance/")

library(ggplot2)
library(reshape2)
library(stringr)
library(purrr)
library(viridis)
library(hrbrthemes)
library(ggpubr)


all_abundace<-list(
  protein=read.csv(file = "./data/input_data/proteins.csv"),
  rna=read.csv(file = "./data/input_data/rna_expression.csv")
)

# Create a list of proteins of interest to ARID1A with their respective UniProt IDs
sorger_feedback <- list(
  # RAS family proteins
  HRAS="P01112",   # HRas protein
  KRAS="P01116",   # KRas protein
  NRAS="P01111",   # NRas protein

  # MAPK pathway components
  ARAF="P10398",   # A-Raf proto-oncogene serine/threonine-protein kinase
  BRAF="P15056",   # B-Raf proto-oncogene serine/threonine-protein kinase
  RAF1="P04049",   # Raf-1 proto-oncogene serine/threonine-protein kinase
  MAP2K1="Q02750", # Mitogen-activated protein kinase kinase 1 (MEK1)
  MAP2K2="P36507", # Mitogen-activated protein kinase kinase 2 (MEK2)
  MAPK3="P27361",  # Mitogen-activated protein kinase 3 (ERK1)
  MAPK1="P28482",  # Mitogen-activated protein kinase 1 (ERK2)

  # Negative feedback regulators of MAPK
  DUSP4="Q13115",  # Dual specificity protein phosphatase 4
  DUSP6="Q16828",  # Dual specificity protein phosphatase 6

  # Negative feedback regulators of EGFR
  ERRFI1="Q9UJM3", # ERBB receptor feedback inhibitor 1
  SPRY4="Q9C004",  # Protein sprouty homolog 4
  SPRY2="O43597",  # Protein sprouty homolog 2

  # EGFR pathway components
  EGFR="P00533",   # Epidermal growth factor receptor
  GRB2="P62993",   # Growth factor receptor-bound protein 2
  CBL="P22681",    # E3 ubiquitin-protein ligase CBL
  SHC1="P29353",   # SHC-transforming protein 1
  PTPN11="Q06124", # Tyrosine-protein phosphatase non-receptor type 11
  SOS1="Q07889"    # Son of sevenless homolog 1
)

# Define the drugs and their corresponding colors
drug_colors <- c(
  "Untreated" = "#A2AEBB",
  "Vemurafenib" = "#FFBA08",
  "Trametinib" = "#D00000",
  "Combination" = "#3F88C5"
)

proteins_df<-melt(all_abundace$protein[all_abundace$protein$X %in% unname(sorger_feedback),])
rna_df<-melt(all_abundace$rna[all_abundace$rna$X %in% names(sorger_feedback),])
replacement_Vec<-c("Untreated","Vermurafenib_1uM","Trametinib_10nM","vemurafenib_and_trametinib")
names(replacement_Vec)<- c("Untreated", "Vemurafenib", "Trametinib", "Combination")

split_exp_conditions<-function(df, error=F){
  df$drug<-unlist(map(str_split(df$variable, pattern = "__"),1))
  df$ko<-unlist(map(str_split(df$variable, pattern = "__"),2))
  df$ko<-gsub("\\..*","",df$ko)
  df$variable<-NULL
  df$drug<-names(replacement_Vec)[match(df$drug, unname(replacement_Vec))]
  # Define the desired order of levels
  desired_order <- c("Untreated", "Vemurafenib", "Trametinib", "Combination")
  # Convert my_column to a factor with the specified order
  df$drug <- factor(df$drug, levels = desired_order)
  if (error) {
    df <- df %>%
      group_by(X, drug, ko) %>%
      summarise(
        n=n(),
        mean=mean(value),
        sd=sd(value)
      ) %>%
      mutate( se=sd/sqrt(n))  %>%
      mutate( ic=se * qt((1-0.05)/2 + .5, n-1))

  }

  return(df)
}


proteins_df<-split_exp_conditions(proteins_df)
rna_df<-split_exp_conditions(rna_df)
rna_df$X <- factor(rna_df$X, levels = names(sorger_feedback))

proteins_df$X<- names(sorger_feedback)[match(proteins_df$X, unname(sorger_feedback))]

complete_df<-rbind(data.frame(proteins_df, data="Protein abundace"),
                   data.frame(rna_df, data="mRNA abundace"))


# Step 1: Call the pdf command to start the plot
pdf(file = "~/Desktop/Melanoma_Resistance/paper/Supplementary_plots/gerosa_negative_feedback.pdf",   # The directory you want to save the file in
    width = 14,  # The width of the plot in inches
    height = 18) # The height of the plot in inches

complete_df %>%
  ggplot(aes(x = value, y = X, fill = drug)) +
  geom_boxplot() +
  geom_jitter(color = "black", size = 0.4, alpha = 0.9) +
  scale_fill_manual(values = drug_colors) +
  ggtitle("mRNA abundances of central RTKs") +
  xlab("Abundance") +  # Change x-axis label to "Abundance"
  ylab("Gene/Protein") +  # Change y-axis label to "Gene/Protein"
  cowplot::theme_cowplot() +
  facet_wrap(data ~ ko, scales = "free") +
  grids(linetype = "dashed") +
  theme(
    plot.title = element_text(size = 15, face = "bold"),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 1)  # Use linewidth instead of size
  )

dev.off()




