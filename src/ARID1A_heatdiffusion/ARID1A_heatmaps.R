
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
  PRKD1 = "Other Kinases",
  PTK2 = "Other Kinases",
  CDK2 = "Other Kinases",

  # RTKs
  EGFR = "RTKs",
  FGFR1 = "RTKs",
  FGFR2 = "RTKs",
  FLT1 = "RTKs",
  KDR = "RTKs",

  # Ephrins
  EPHB1 = "Ephrin receptors",
  EPHA7 = "Ephrin receptors",
  EPHA8 = "Ephrin receptors",
  EPHA4 = "Ephrin receptors",
  EPHA2 = "Ephrin receptors",
  EPHA3 = "Ephrin receptors"
)

experiments_to_see <- c("ARID1A Combined vs Untreated", "Combined vs Untreated", "Untreated ARID1A vs WT")
heatmap_data <- experiment_matrix[match(names(nodes_to_see), experiment_matrix$`Kinase Name`),]
abundance_plot<-heatmap_data$`Kinase Uniprot ID`
names(abundance_plot)<-heatmap_data$`Kinase Name`
rownames(heatmap_data) <- heatmap_data$`Kinase Name`


# Extract the matrix for heatmap
data_matrix <- data.matrix(heatmap_data[, -1])
rownames(data_matrix) <- heatmap_data$`Kinase Name`

data_matrix<-t(data_matrix[,!colnames(data_matrix) %in% 'Kinase Name'])
kinase_activity<-Heatmap(data.matrix(data_matrix[experiments_to_see,]),
                         name = "Median Kinase Activity (LFC)",
                         na_col = "darkgrey",
                         rect_gp = gpar(col = "white", lwd = 5), cluster_columns = F, cluster_rows = F,
                         column_gap=unit(.02, "npc"), row_title_gp = gpar(fontsize = 10),
                         column_title_gp = gpar(fontsize = 10), column_split = unname(nodes_to_see), heatmap_width = unit(25, "cm")
)


all_lfcs<-list(
  protein=read.csv(file = "./results/lfc/protein/arid1a_lfc.csv"),
  rna=read.csv(file = "./results/lfc/mRNA//arid1a_lfc.csv")
)


# Define the categories for each node
nodes_to_see_proteins <- c(
  # JNKs
  MAPK8 = "JNK pathway",
  MAPK9 = "JNK pathway",
  MAPK10 = "JNK pathway",

  # MAPKs
  MAPK1 = "MAPK pathway",
  MAPK3 = "MAPK pathway",
  MAPK7 = "MAPK pathway",

  # RTKs
  EGFR = "RTK",
  FGFR1 = "RTK",
  FGFR2 = "RTK",
  KDR = "RTK",
  FLT1 = "RTK",
  ROS1 = "RTK",

  # Other Kinases
  PRKD1 = "Other kinase",
  PTK2 = "Other kinase",
  CDK2 = "Other kinase",
  FYN = "Other kinase",

  # Ephrins
  EPHA2 = "Ephrin receptors",
  EPHA3 = "Ephrin receptors",
  EPHA4 = "Ephrin receptors",
  EPHA7 = "Ephrin receptors",
  EPHA8 = "Ephrin receptors",
  EPHB1 = "Ephrin receptors",

  # Transcription Factors
  JUN = "Transcription factor",
  MYC = "Transcription factor",

  # Others
  CD44 = "Other",
  ITGA4 = "Other"
)


all_lfcs_df<-bind_rows(all_lfcs, .id = "data_type")
all_lfcs_df<-all_lfcs_df[all_lfcs_df$X %in% names(nodes_to_see_proteins),]


# Create logFC matrix
logFC_matrix <- all_lfcs_df %>%
  dplyr::select(data_type, Gene, logFC) %>%
  pivot_wider(names_from = Gene, values_from = logFC)

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
significance_matrix[is.na(significance_matrix)] <- ""
data_matrix <- data.matrix(logFC_matrix[, -1])
rownames(data_matrix) <- logFC_matrix$data_type

library(ComplexHeatmap)
library(grid)

# Create the heatmap
significance_matrix_to_plot<-significance_matrix[-1]

data_matrix<-data_matrix[, names(nodes_to_see_proteins)]
significance_matrix_to_plot<-significance_matrix_to_plot[, names(nodes_to_see_proteins)]

protein_heatmap<-Heatmap(t(data_matrix[2,]), cluster_rows = FALSE, cluster_columns = FALSE,
        cell_fun = function(j, i, x, y, width, height, fill) {
          # Use the significance_matrix directly to add significance labels
          grid.text(as.character(t(significance_matrix_to_plot[2,colnames(data_matrix)][i, j])), x, y, gp = gpar(fontsize = 10))
        },
        name = "Protein abundance (LFC)",
        na_col = "darkgrey",
        rect_gp = gpar(col = "white", lwd = 5),
        column_gap=unit(.02, "npc"), row_title_gp = gpar(fontsize = 10),
        column_title_gp = gpar(fontsize = 10), column_split = unname(nodes_to_see_proteins),
        row_labels = "ARID1A KO vs Parental (Protein Abundance)", heatmap_width = unit(28, "cm"))
rna_heatmap<-Heatmap(t(data_matrix[1,]), cluster_rows = FALSE, cluster_columns = FALSE,
                         cell_fun = function(j, i, x, y, width, height, fill) {
                           # Use the significance_matrix directly to add significance labels
                           grid.text(as.character(t(significance_matrix_to_plot[1,colnames(data_matrix)][i, j])), x, y, gp = gpar(fontsize = 10))
                         },
                         name = "mRNA abundance (LFC)",
                         na_col = "darkgrey",
                         rect_gp = gpar(col = "white", lwd = 5),
                         column_gap=unit(.02, "npc"), row_title_gp = gpar(fontsize = 10),
                         column_title_gp = gpar(fontsize = 10), column_split = unname(nodes_to_see_proteins),
                         row_labels = "ARID1A KO vs Parental (mRNA Abundance)", heatmap_width = unit(28, "cm"))


# Combine the main heatmap and the significance annotation
combined_heatmap <- protein_heatmap %v% rna_heatmap


# Step 2: Call the pdf command to start the plot
pdf(file = "~/Desktop/Melanoma_Resistance/paper/Figures/ARID1A_heatmap.pdf",
    width = 18,  # The calculated width of the plot in inches
    height = 2) # The calculated height of the plot in inches

# Draw the combined heatmap
draw(combined_heatmap, merge_legend = TRUE)
kinase_activity
# Step 3: Close the PDF device
dev.off()
