
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
library(tidyr)
library(EnsDb.Hsapiens.v86)
library(ComplexHeatmap)
library(reshape2)
library(grid)


kinomics_files<-list.files(path = "./data/kinomics//", full.names = T, recursive = T)

results_list<-list()
for (UKA_file in kinomics_files[grep(kinomics_files, pattern = "vs")]) {
  kinomics_results<-read_xlsx(path =UKA_file)
  results_list[[UKA_file]]<-kinomics_results
}
# Apply the conversion to each data frame in the list
results_list <- lapply(results_list, function(x) {
  x$`SD Kinase Statitistic` <- as.numeric(x$`SD Kinase Statitistic`)
  return(x)
})

# Bind the rows
complete_results <- bind_rows(results_list, .id = "file")
complete_results$experiment<-unlist(map(str_split(complete_results$file, pattern = "/"), last))
complete_results$background <- unlist(map(str_split(complete_results$file, pattern = "/"), 6))

# Assuming `complete_results` is your data frame containing the required columns

# Reorder Kinase Name by Mean Kinase Statistic
complete_results$`Kinase Name` <- reorder(complete_results$`Kinase Name`, complete_results$`Mean Kinase Statistic`)
complete_results$arid1a_status<-unlist(map(str_split(complete_results$background, pattern = "__"),2))
complete_results$kinase_type<-unlist(map(str_split(complete_results$background, pattern = "__"),1))
complete_results$file<- unlist(map(str_split(complete_results$file, "/"), last))
# Transform the dataset
experiment_matrix <- complete_results %>%
  # Spread the data into a matrix format
  dplyr::select(`Kinase Uniprot ID`, file, `Median Kinase Statistic`) %>%
  pivot_wider(names_from = file, values_from = `Median Kinase Statistic`)
genename_df <- AnnotationDbi::select(EnsDb.Hsapiens.v86, keys = experiment_matrix$`Kinase Uniprot ID`, keytype = "UNIPROTID", columns = "GENENAME");
experiment_matrix$`Kinase Name` <- genename_df$GENENAME[match(experiment_matrix$`Kinase Uniprot ID`, genename_df$UNIPROTID)]

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

# Apply the function to column names
simplified_colnames <- sapply(colnames(experiment_matrix), simplify_column_names)

# Set the new column names
colnames(experiment_matrix) <- simplified_colnames
#get list of receptors

# Define the categories for each node
nodes_to_see <- c(
  # JNKs
  MAPK8 = "JNK pathway",
  MAPK9 = "JNK pathway",
  MAPK10 = "JNK pathway",

  #MAPKs
  MAPK1 = "MAPK pathway",
  MAPK3 = "MAPK pathway",
  MAPK7 = "MAPK pathway",

  # JNKs
  FYN = "Other Kinases",
  ABL2 = "Other Kinases",
  GSK3B = "Other Kinases",
  PRKD1 = "Other Kinases",
  PTK2 = "Other Kinases",
  TXK = "Other Kinases",

  # RTKs
  EGFR = "RTKs",
  FGFR1 = "RTKs",
  FGFR2 = "RTKs",
  FLT1 = "RTKs",
  FLT4 = "RTKs",
  IGF1R = "RTKs",
  INSR = "RTKs",
  KIT = "RTKs",
  MET = "RTKs",
  NTRK2 = "RTKs",
  TEK = "RTKs",

  # Phosphatases
  SPRED1 = "Phosphatases",
  SPRED2 = "Phosphatases",
  SPRY1 = "Phosphatases",
  SPRY2 = "Phosphatases",
  SPRY4 = "Phosphatases",
  DUSP1 = "Phosphatases",
  DUSP2 = "Phosphatases",
  DUSP4 = "Phosphatases",
  DUSP6 = "Phosphatases"
)

experiments_to_see <- c("Trametinib.vs.Untreated", "Combined.vs.Untreated", "Vemurafenib.vs.Untreated")
heatmap_data <- data.frame(experiment_matrix[match(names(nodes_to_see)[nodes_to_see != "Phosphatases"], experiment_matrix$`Kinase Name`),])
abundance_plot<-heatmap_data$`Kinase Uniprot ID`
names(abundance_plot)<-heatmap_data$`Kinase Name`
rownames(heatmap_data) <- heatmap_data$`Kinase Name`


# Extract the matrix for heatmap
data_matrix <- data.matrix(heatmap_data)
rownames(data_matrix) <- heatmap_data$Kinase.Name
colnames(data_matrix) <- colnames(heatmap_data)

data_matrix<-t(data_matrix[,!colnames(data_matrix) %in% c('Kinase.Name', 'Kinase.Uniprot.ID')])
kinase_activity<-Heatmap(data.matrix(data_matrix[experiments_to_see,]),
                         name = "Median Kinase Activity (LFC)",
                         na_col = "darkgrey",
                         rect_gp = gpar(col = "white", lwd = 5), cluster_columns = F, cluster_rows = F,
                         column_gap=unit(.02, "npc"), row_title_gp = gpar(fontsize = 10),
                         column_title_gp = gpar(fontsize = 10),
                         column_split = unname(nodes_to_see)[nodes_to_see != "Phosphatases"],
                         heatmap_width = unit(25, "cm")
)


mRNA_LFC<-list(
  `Combination vs Untreated`=read.csv(file = "./results/lfc/mRNA//combination_lfc.csv"),
  `Trametinib vs Untreated`=read.csv(file = "./results/lfc/mRNA//trametinib_lfc.csv"),
  `Vemurafenib vs Untreated`=read.csv(file = "./results/lfc/mRNA//vemurafenib_lfc.csv")
)
protein_LFC<-list(
  `Combination vs Untreated`=read.csv(file = "./results/lfc/protein/combination_lfc.csv"),
  `Trametinib vs Untreated`=read.csv(file = "./results/lfc/protein/trametinib_lfc.csv"),
  `Vemurafenib vs Untreated`=read.csv(file = "./results/lfc/protein/vemurafenib_lfc.csv")
)

process_LFC_data <- function(LFC_list, nodes_to_see) {
  # Bind rows into a single dataframe
  all_lfcs_df <- bind_rows(LFC_list, .id = "data_type")

  # Filter for rows where `X` is in `nodes_to_see`
  all_lfcs_df <- all_lfcs_df[all_lfcs_df$X %in% names(nodes_to_see), ]

  # Ensure `X` has consistent levels based on `nodes_to_see`
  all_lfcs_df$X <- factor(all_lfcs_df$X, levels = names(nodes_to_see))

  # Create logFC matrix
  logFC_matrix <- all_lfcs_df %>%
    dplyr::select(data_type, X, logFC) %>%
    pivot_wider(names_from = X, values_from = logFC)

  # Add missing columns as `NA` to ensure all `nodes_to_see` are included
  missing_cols <- setdiff(names(nodes_to_see), colnames(logFC_matrix))
  logFC_matrix[missing_cols] <- NA

  # Reorder columns based on `nodes_to_see`
  logFC_matrix <- logFC_matrix %>%
    dplyr::select(data_type, !!!syms(names(nodes_to_see)))

  # Create the significance matrix
  significance_matrix <- all_lfcs_df %>%
    dplyr::select(Gene, adj.P.Val, data_type) %>%
    mutate(Significance = case_when(
      adj.P.Val < 0.001 ~ "***",
      adj.P.Val < 0.01 ~ "**",
      adj.P.Val < 0.1 ~ "*",
      TRUE ~ ""
    )) %>%
    dplyr::select(-adj.P.Val) %>%
    # Ensure we have unique values for each gene
    group_by(Gene, data_type) %>%
    summarise(Significance = max(Significance), .groups = 'drop') %>%
    pivot_wider(names_from = Gene, values_from = Significance)

  # Add missing columns to the significance matrix as `NA`
  missing_significance_cols <- setdiff(names(nodes_to_see), colnames(significance_matrix))
  significance_matrix[missing_significance_cols] <- NA

  # Replace NA values with an empty string
  significance_matrix[is.na(significance_matrix)] <- ""

  # Reorder significance matrix columns
  significance_matrix <- significance_matrix %>%
    dplyr::select(data_type, !!!syms(names(nodes_to_see)))

  # Convert logFC matrix to data matrix
  data_matrix <- data.matrix(logFC_matrix[, -1])
  rownames(data_matrix) <- logFC_matrix$data_type

  # Adjust significance matrix for plotting
  significance_matrix_to_plot <- significance_matrix[-1]

  # Return a list containing data matrix and significance matrix
  return(list(
    data_matrix = data_matrix,
    significance_matrix_to_plot = significance_matrix_to_plot
  ))
}


mRNA_to_plot<-process_LFC_data(mRNA_LFC, nodes_to_see)
protein_to_plot<-process_LFC_data(protein_LFC, nodes_to_see)

# Reorder data matrices to ensure phosphatases come first
mRNA_to_plot$data_matrix <- mRNA_to_plot$data_matrix[,order(colnames(mRNA_to_plot$data_matrix) == "Phosphatases", decreasing = F)]
protein_to_plot$data_matrix <- protein_to_plot$data_matrix[,order(colnames(protein_to_plot$data_matrix) == "Phosphatases", decreasing = F)]
mRNA_to_plot$significance_matrix_to_plot <- mRNA_to_plot$significance_matrix_to_plot[,order(colnames(mRNA_to_plot$significance_matrix_to_plot) == "Phosphatases", decreasing = F)]
protein_to_plot$significance_matrix_to_plot <- protein_to_plot$significance_matrix_to_plot[, order(colnames(protein_to_plot$significance_matrix_to_plot) == "Phosphatases", decreasing = F)]

# Create the mRNA heatmap
mRNA_heatmap <- Heatmap(
  mRNA_to_plot$data_matrix,
  cluster_rows = FALSE, cluster_columns = FALSE,
  cell_fun = function(j, i, x, y, width, height, fill) {
    grid.text(as.character(mRNA_to_plot$significance_matrix_to_plot[i, j]), x, y, gp = gpar(fontsize = 10))
  },
  name = "mRNA abundance (LFC)",
  na_col = "darkgrey",
  rect_gp = gpar(col = "white", lwd = 5),
  column_gap = unit(.02, "npc"),
  row_title_gp = gpar(fontsize = 10),
  column_title_gp = gpar(fontsize = 10),
  column_split = unname(nodes_to_see)[match(colnames(mRNA_to_plot$data_matrix), names(nodes_to_see))],
  heatmap_width = unit(28, "cm")
)

# Create the protein heatmap
protein_heatmap <- Heatmap(
  protein_to_plot$data_matrix,
  cluster_rows = FALSE, cluster_columns = FALSE,
  cell_fun = function(j, i, x, y, width, height, fill) {
    grid.text(as.character(protein_to_plot$significance_matrix_to_plot[i, j]), x, y, gp = gpar(fontsize = 10))
  },
  name = "Protein abundance (LFC)",
  na_col = "darkgrey",
  rect_gp = gpar(col = "white", lwd = 5),
  column_gap = unit(.02, "npc"),
  row_title_gp = gpar(fontsize = 10),
  column_title_gp = gpar(fontsize = 10),
  column_split = unname(nodes_to_see)[match(colnames(protein_to_plot$data_matrix), names(nodes_to_see))],
  heatmap_width = unit(28, "cm")
)

# Combine the heatmaps with phosphatases first
combined_heatmap <- mRNA_heatmap %v% protein_heatmap

# Generate the PDF
pdf(file = "~/Desktop/Melanoma_Resistance/paper/Figures/drug_agnostic_heatmap.pdf",
    width = 20,  # The calculated width of the plot in inches
    height = 3)  # The calculated height of the plot in inches

# Draw the combined heatmap
draw(combined_heatmap, merge_legend = TRUE)
dev.off()

pdf(file = "~/Desktop/Melanoma_Resistance/paper/Figures/drug_agnostic_kinases.pdf",
    width = 20,  # The calculated width of the plot in inches
    height = 2.2)  # The calculated height of the plot in inches

kinase_activity
# Close the PDF device
dev.off()
