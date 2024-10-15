
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
library(ggrepel)
library(wesanderson)
library(igraph)
library(purrr)
library(stringr)
library(cowplot)
library(dplyr)
library(tidyr)
library(EnsDb.Hsapiens.v86)

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

#get phuego graphs
KDE <- "0.5"
factorS <- c("Factor1", "Factor2", "Factor3")
results_dir <- "./results/phuego/results/"
factor_graphs<-list()

#provisionally we need to find a way to make the phuego networks directed.

for (factor in factorS) {
  file_graphml_up<-paste0(results_dir, factor, "/increased/KDE_", KDE, "/networks/KDE.graphml")
  file_graphml_down<-paste0(results_dir, factor, "/decreased/KDE_", KDE, "/networks/KDE.graphml")
  factor_graphs[[factor]][["up"]]<-read_graph(file=file_graphml_up,format = "graphml")
  factor_graphs[[factor]][["down"]]<-read_graph(file=file_graphml_down,format = "graphml")
  E(factor_graphs[[factor]][["up"]])$source <- "phuego"
  E(factor_graphs[[factor]][["down"]])$source <- "phuego"
}

get_nodes<-function(igraph_object, conv=F){
  nodes<-V(igraph_object)$name
  weights<-igraph::page_rank(igraph_object)

  if (conv) {
    genename_df <- AnnotationDbi::select(EnsDb.Hsapiens.v86, keys = nodes,
                                         keytype = "UNIPROTID",
                                         columns = "GENENAME")
    return(sort(genename_df$GENENAME))
  }
  else {
    return(sort(weights$vector))
  }
}

#get nodes from all graphs
factor_nodes <- lapply(factor_graphs, function(sublist) {
  sapply(sublist, get_nodes, T)
})


complete_results$arid1a_status<-unlist(map(str_split(complete_results$background, pattern = "__"),2))
complete_results$kinase_type<-unlist(map(str_split(complete_results$background, pattern = "__"),1))


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
experiment_matrix$in_network <- experiment_matrix$`Kinase Uniprot ID` %in% unlist(lapply(kinase_list_all_factors, rownames))
library(ggplot2)
library(ggrepel)
library(cowplot)
library(rlang)
library(ghibli)

# Define a function for creating the plots to avoid redundancy
create_plot <- function(data, kinase_type) {
  x_col <- sym(paste0(kinase_type, " - ARID1A_KO/ARID1A Combined vs Untreated"))
  y_col <- sym(paste0(kinase_type, " - WT/Combined vs Untreated"))

  ggplot(data, aes(x = !!x_col, y = !!y_col)) +
    geom_point(colour="lightblue") +
    cowplot::theme_cowplot() +
    theme(
      plot.title = element_text(size = 15, face = "bold"),
      panel.border = element_rect(colour = "black", fill = NA, linewidth = 1),
    ) +
    grids(linetype = "dashed") +
    xlab(paste("Kinase Activity in ARID1A (Combined vs Untreated)")) +
    ylab("Kinase Activity in WT (Combined vs Untreated)") +
    ggtitle(paste(kinase_type, "Kinase Activity Comparison")) +

    # Add axis lines
    geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
    geom_vline(xintercept = 0, linetype = "dashed", color = "black") +

    # Add text labels for points where values exceed thresholds
    geom_text_repel(
      aes(label = ifelse(
        (!!x_col > 1 | !!x_col < -1 | !!y_col > 1 | !!y_col < -1),
        `Kinase Name`,
        NA)),  # Labels are shown only for selected points
      size = 3,  # Adjust text size as needed
      box.padding = 0.35,  # Padding around the text label to avoid overlap
      point.padding = 0.5,  # Space between the label and the point
      max.overlaps = Inf,  # Avoid limit on label overlaps
      min.segment.length = 0  # Draw segment lines for all labels
    )+
    # ghibli stuff
    scale_colour_ghibli_d("LaputaMedium", direction = -1)
}

# Create plots for PTK and STK
ptk_plot <- create_plot(experiment_matrix, "PTK")
stk_plot <- create_plot(experiment_matrix, "STK")
# Step 1: Call the pdf command to start the plot
pdf(file = "~/Desktop/Melanoma_Resistance/paper/Supplementary_plots/kinomics_arid1a.pdf",   # The directory you want to save the file in
    width = 10,  # The width of the plot in inches
    height = 6) # The height of the plot in inches

# Display the plots
plot_grid(ptk_plot, stk_plot, labels = c("A", "B"))

dev.off()
