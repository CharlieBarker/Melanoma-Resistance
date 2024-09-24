
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
  sapply(sublist, get_nodes)
})


# Assuming kinase_list is a vector of kinase names you want to mark with stars
kinase_list_all_factors <- list(Factor1=stack(factor_nodes$Factor1), #factor 1 describes changes in all drugs
                                Factor2=stack(factor_nodes$Factor2), #factor 2 describes changes specific to combination
                                Factor3=stack(factor_nodes$Factor3) #factor 3 describes ARID1A
)
#this one works, we know at least

complete_results$arid1a_status<-unlist(map(str_split(complete_results$background, pattern = "__"),2))
complete_results$kinase_type<-unlist(map(str_split(complete_results$background, pattern = "__"),1))

genename_df <- AnnotationDbi::select(EnsDb.Hsapiens.v86, keys = complete_results$`Kinase Uniprot ID`, keytype = "UNIPROTID", columns = "GENENAME");
complete_results$Gene_name <- genename_df$GENENAME[match(complete_results$`Kinase Uniprot ID`, genename_df$UNIPROTID)]
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
simplified_colnames <- sapply(complete_results$file, simplify_column_names)

# Set the new column names
complete_results$experiment <- unlist(map(str_split(simplified_colnames, pattern = " - "),2))

kinases_of_interest<-c("MAPK1", "MAPK3", "PRKD1", "FYN", "MAPK8", "MAPK9", "MAPK10")
drug_targets<-complete_results[complete_results$Gene_name %in% kinases_of_interest,]
pal <- wes_palette("Zissou1", 100, type = "continuous")

drug_targets <- drug_targets[drug_targets$experiment %in% c("WT/Combined vs Untreated", "ARID1A_KO/ARID1A Combined vs Untreated"),]
drug_targets$Gene_name <- factor(drug_targets$Gene_name, levels=kinases_of_interest)

pdf(file = "~/Desktop/Melanoma_Resistance/results/kinomics_microarray/prkd1_story.pdf",   # The directory you want to save the file in
    width = 6,  # The width of the plot in inches
    height = 8) # The height of the plot in inches
drug_targets$experiment[drug_targets$experiment=="ARID1A_KO/ARID1A Combined vs Untreated"] <- "ARID1A KO/Combined vs Untreated"
ggplot(drug_targets[!is.na(drug_targets$experiment),], aes(x=Gene_name, y=`Median Kinase Statistic`,
                                                           colour=`Median Kinase Statistic`)) +
  scale_colour_gradientn(colours = pal) +
  geom_hline(yintercept = 0, color = "black") + # Add horizontal line at y=0
  geom_segment(aes(x=Gene_name, xend=Gene_name, y=0, yend=`Median Kinase Statistic`), color="grey") +
  geom_point(aes(size = `Mean Specificity Score`)) +  # Adjust the size of points here
  cowplot::theme_cowplot() +
  theme(
    plot.title = element_text(size = 15, face = "bold"),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 1),
    axis.text.x = element_text(angle = 90, hjust = 1) # Rotate x-axis labels
  ) +
  grids(linetype = "dashed") +
  facet_wrap(~experiment, scales = "free_x", nrow = 2) +
  xlab("") +
  ylab("Median Kinase Statistic")


dev.off()
