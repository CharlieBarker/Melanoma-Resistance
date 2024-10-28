
#find correct environment
packLib="/usr/lib/R"
if (file.exists(packLib)) {
  reticulate::use_condaenv("py37", required = T)
  setwd(dir = "~/Desktop/Melanoma_Resistance/")
}else {
  setwd(dir = "~/Desktop/Melanoma_Resistance/")
}

library(readxl)
library(ggplot2)
library(dplyr)
library(ggpubr)
library(diffusr)
library(ggrepel)
library(wesanderson)
library(igraph)
library(purrr)
library(stringr)
library(cowplot)
library(dplyr)
library(tidyr)
library(EnsDb.Hsapiens.v86)
library(ComplexHeatmap)
library(reshape2)

# Function to simplify the column names
simplify_column_names <- function(colname) {
  colname %>%
    # Remove the initial path "./data/kinomics///"
    gsub("^.*/kinomics///", "", .) %>%
    # Remove the ".xlsx" extension
    gsub("\\.xlsx", "", .) %>%
    # Replace "__" with " - " for clarity
    gsub("__", " - ", .)
}


all_lfcs<-list(
  rna=read.csv(file = "./results/lfc/mRNA//arid1a_lfc.csv"),
  tcga=read.csv(file = "./results/tcga_lfc/arid1a_mutant_vs_wt.csv")
)


# Define the categories for each node
nodes_to_see_proteins <- c(
  # Cadherins
  CDH11 = "Cadherins",
  CDH19 = "Cadherins",
  CDH4 =  "Cadherins",
  CDH6 =  "Cadherins",


  # Collagens
  COL1A2 = "Collagens",
  COL1A1 = "Collagens",
  COL3A1 = "Collagens",
  COL11A2 = "Collagens",
  COL5A3 = "Collagens",
  COL5A1 = "Collagens",

  # ECM
  TGFBI = "ECM",
  ANK1 = "ECM",
  ITGA8 = "ECM",
  LAMC3 = "ECM",
  LARP1 = "ECM",
  NRG1 = "ECM",
  RELN = "ECM",
  SPP1 = "ECM",

  # HLA Class II
  CD74 = "HLA Class II",
  CIITA = "HLA Class II",
  HLA.DMB = "HLA Class II",
  HLA.DOA = "HLA Class II",
  HLA.DQA1 = "HLA Class II",
  HLA.DRA = "HLA Class II",

  # L1CAM interactions
  L1CAM = "L1CAM interactions",
  NRP2 = "L1CAM interactions",
  CNTN1 = "L1CAM interactions",
  NCAM1 = "L1CAM interactions",

  # MMPs
  MMP1 = "MMPs",
  MMP3 = "MMPs",
  MMP10 = "MMPs",
  MMP16 = "MMPs",
  MMP11 = "MMPs"


)

all_lfcs$tcga$X<-gsub(x = all_lfcs$tcga$X, pattern = "-", replacement = ".")
all_lfcs_df<-bind_rows(all_lfcs, .id = "data_type")
all_lfcs_df<-all_lfcs_df[all_lfcs_df$X %in% names(nodes_to_see_proteins),]


# Create logFC matrix
logFC_matrix <- all_lfcs_df %>%
  dplyr::select(data_type, X, logFC) %>%
  pivot_wider(names_from = X, values_from = logFC)

# Create the significance matrix
significance_matrix <- all_lfcs_df %>%
  dplyr::select(X, adj.P.Val, data_type) %>%
  mutate(Significance = case_when(
    adj.P.Val < 0.001 ~ "***",
    adj.P.Val < 0.01 ~ "**",
    adj.P.Val < 0.1 ~ "*",
    TRUE ~ ""
  )) %>%
  dplyr::select(-adj.P.Val) %>%
  # Ensure we have unique values for each gene
  group_by(X, data_type) %>%
  summarise(Significance = max(Significance), .groups = 'drop') %>%
  pivot_wider(names_from = X, values_from = Significance)
significance_matrix[is.na(significance_matrix)] <- ""
data_matrix <- data.matrix(logFC_matrix[, -1])
rownames(data_matrix) <- logFC_matrix$data_type

library(ComplexHeatmap)
library(grid)

# Create the heatmap
significance_matrix_to_plot<-data.frame(significance_matrix[-1])
rownames(significance_matrix_to_plot) <- significance_matrix[1]$data_type

data_matrix<-data_matrix[, names(nodes_to_see_proteins)]
significance_matrix_to_plot<-significance_matrix_to_plot[, names(nodes_to_see_proteins)]


tcga_heatmap<-Heatmap(t(data_matrix["tcga",]), cluster_rows = FALSE, cluster_columns = FALSE,
                         cell_fun = function(j, i, x, y, width, height, fill) {
                           # Use the significance_matrix directly to add significance labels
                           grid.text(as.character(t(significance_matrix_to_plot["tcga",colnames(data_matrix)][i, j])), x, y, gp = gpar(fontsize = 10))
                         },
                         name = "Protein abundance (LFC)",
                         na_col = "darkgrey",
                         rect_gp = gpar(col = "white", lwd = 5),
                         column_gap=unit(.02, "npc"), row_title_gp = gpar(fontsize = 10),
                         column_title_gp = gpar(fontsize = 10), column_split = unname(nodes_to_see_proteins),
                         row_labels = "ARID1A mutant patients vs Unmutated (Protein Abundance)", heatmap_width = unit(28, "cm"))

rna_heatmap<-Heatmap(t(data_matrix["rna",]), cluster_rows = FALSE, cluster_columns = FALSE,
                         cell_fun = function(j, i, x, y, width, height, fill) {
                           # Use the significance_matrix directly to add significance labels
                           grid.text(as.character(t(significance_matrix_to_plot["rna",colnames(data_matrix)][i, j])), x, y, gp = gpar(fontsize = 10))
                         },
                         name = "Protein abundance (LFC)",
                         na_col = "darkgrey",
                         rect_gp = gpar(col = "white", lwd = 5),
                         column_gap=unit(.02, "npc"), row_title_gp = gpar(fontsize = 10),
                         column_title_gp = gpar(fontsize = 10), column_split = unname(nodes_to_see_proteins),
                         row_labels = "ARID1A KO vs Parental (RNA Abundance)", heatmap_width = unit(28, "cm"))

# Combine the main heatmap and the significance annotation
combined_heatmap <-  rna_heatmap  %v% tcga_heatmap


# Step 2: Call the pdf command to start the plot
pdf(file = "~/Desktop/Melanoma_Resistance/paper/Figures/ARID1A_tcga_heatmap.pdf",
    width = 18,  # The calculated width of the plot in inches
    height = 2) # The calculated height of the plot in inches

# Draw the combined heatmap
draw(combined_heatmap, merge_legend = TRUE)

# Step 3: Close the PDF device
dev.off()
