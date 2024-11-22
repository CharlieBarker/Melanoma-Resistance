

# Load required libraries
library(readxl)
library(ggplot2)
library(dplyr)
library(ggpubr)
library(ggrepel)
library(wesanderson)
library(ggraph)
library(igraph)
library(cowplot)
library(purrr)
library(stringr)
library(ghibli)
library(EnsDb.Hsapiens.v86)
library(tidyr)
library(MOFA2)
library(ggConvexHull)


knitr::opts_chunk$set(dev = "ragg_png")

# Set environment and working directory
packLib <- "/usr/lib/R"
setwd(dir = "~/Desktop/Melanoma_Resistance/")
if (file.exists(packLib)) {
  reticulate::use_condaenv("py37", required = TRUE)
}

source("./src/functions/default_variables.R")
pal <- wes_palette("Zissou1", 100, type = "continuous")

# Read kinomics data files
kinomics_files <- list.files(path = "./data/kinomics/", full.names = TRUE, recursive = TRUE)

# Read and process kinomics data files
results_list <- kinomics_files %>%
  keep(~ grepl("vs", .)) %>%
  map(~ read_xlsx(.)) %>%
  map(~ mutate(.x, `SD Kinase Statitistic` = as.numeric(`SD Kinase Statitistic`)))

# Bind the rows
complete_results <- bind_rows(results_list, .id = "file")

# Extract experiment and background information
complete_results <- complete_results %>%
  mutate(
    experiment = map_chr(file, ~ tail(str_split(.x, pattern = "/")[[1]], 1)),
    background = map_chr(file, ~ str_split(.x, pattern = "/")[[1]][6])
  )

# Reorder Kinase Name by Mean Kinase Statistic
complete_results$`Kinase Name` <- reorder(complete_results$`Kinase Name`, complete_results$`Mean Kinase Statistic`)

# Get phuego graphs
KDE <- "0.5"
factorS <- c("Factor1", "Factor2", "Factor3")
results_dir <- "./results/phuego/results/"
factor_graphs <- list()
factor_genes_names <- list()

factor_to_vis<-"Factor1"

# Process each factor
for (factor in factorS) {
  file_graphml_up <- file.path(results_dir, factor, "increased", paste0("KDE_", KDE), "networks", "KDE.graphml")
  file_graphml_down <- file.path(results_dir, factor, "decreased", paste0("KDE_", KDE), "networks", "KDE.graphml")

  factor_graphs[[factor]][["up"]] <- read_graph(file = file_graphml_up, format = "graphml")
  factor_graphs[[factor]][["down"]] <- read_graph(file = file_graphml_down, format = "graphml")

  for (direction in c("up", "down")) {
    V(factor_graphs[[factor]][[direction]])$source <- "phuego"
    V(factor_graphs[[factor]][[direction]])$direction <- direction
    V(factor_graphs[[factor]][[direction]])$factor <- factor
    factor_genes_names[[factor]][[direction]] <- V(factor_graphs[[factor]][[direction]])$Gene_name

    # Ensure the graph is directed
    if (!is.directed(factor_graphs[[factor]][[direction]])) {
      factor_graphs[[factor]][[direction]] <- as.directed(factor_graphs[[factor]][[direction]], mode = "mutual")
    }
  }
}

factor_genes_names_df<-lapply(factor_genes_names, stack)
factor_genes_names_df <- bind_rows(factor_genes_names_df, .id = "factor")

# factor_graphs

# ├── Factor1
# │   ├── up   -> igraph object for Factor1 increased condition
# │   └── down -> igraph object for Factor1 decreased condition
# ├── Factor2
# │   ├── up   -> igraph object for Factor2 increased condition
# │   └── down -> igraph object for Factor2 decreased condition
# └── Factor3
# ├── up   -> igraph object for Factor3 increased condition
# └── down -> igraph object for Factor3 decreased condition

#get mofa weights
MOFAobject.trained<-load_model(file = "./results/mofa/mofa_object.hdf5")

weights <- get_weights(MOFAobject.trained,
                       views = "all",
                       as.data.frame = TRUE
)
#for clarity swap Factor1 so it means going down with treatment
weights[weights$factor=="Factor1",]$value<-weights[weights$factor=="Factor1",]$value * -1
new_up<-factor_graphs$Factor1$down
new_down<-factor_graphs$Factor1$up
factor_graphs$Factor1$down<-new_down
factor_graphs$Factor1$up<-new_up


weights$node<-unlist(map(str_split(weights$feature, pattern = "_"),1))
weights$node<-unlist(map(str_split(weights$node, pattern = ";"),1))
factor_weights<-weights[weights$factor == factor_to_vis,]
# Collapse the data frame by summing the value for each node
collapsed_factor_weights <- factor_weights %>%
  group_by(node) %>%
  summarise(value = sum(value))

# Define the centrality_plot function with an additional title parameter
centrality_plot <- function(graph_in, title, factor_weights) {
  # Get centrality
  centrality_out <- page.rank(graph_in)
  centrality_df <- stack(centrality_out$vector)

  # Order and rank centrality
  centrality_df <- centrality_df[order(centrality_df$values, decreasing = TRUE), ]
  centrality_df$rank <- 1:nrow(centrality_df)

  # Retrieve gene names
  genename_df <- AnnotationDbi::select(EnsDb.Hsapiens.v86, keys = as.character(centrality_df$ind),
                                       keytype = "UNIPROTID",
                                       columns = "GENENAME")
  centrality_df$names <- genename_df$GENENAME[match(centrality_df$ind, genename_df$UNIPROTID)]

  centrality_df$factor_weight <- factor_weights$value[match(centrality_df$names, factor_weights$node)]

  centrality_df$names[centrality_df$rank > 10] <- ""

  # Create the plot with the given title
  pal <- wes_palette("Zissou1", 100, type = "continuous")
  # Define your plot with ggplot
  out <- ggplot(centrality_df, aes(x = rank, y = values, label = names, colour = factor_weight)) +
    geom_point() +
    cowplot::theme_cowplot() +
    geom_text_repel(colour = "black") +
    scale_colour_gradientn(colours = pal) +
    ggtitle(title) +  # Add the plot title
    labs(x = "Rank in network", y = "Centrality (PageRank)") +  # Adding axis titles
    theme(legend.position = "bottom")  # Placing legend at the bottom


  return(centrality_df)
}

# Create the centrality plots with titles
factor1_up_centrality<-centrality_plot(factor_graphs[[factor_to_vis]]$up, paste0(factor_to_vis, " Up - Centrality Plot"), collapsed_factor_weights)
factor1_down_centrality<-centrality_plot(factor_graphs[[factor_to_vis]]$down, paste0(factor_to_vis, " Down - Centrality Plot"), collapsed_factor_weights)
factor_centrality<-rbind(factor1_up_centrality, factor1_down_centrality)


# Order and rank centrality
factor_centrality <- factor_centrality[order(factor_centrality$values, decreasing = TRUE), ]
factor_centrality$rank <- 1:nrow(factor_centrality)

# Retrieve gene names
genename_df <- AnnotationDbi::select(EnsDb.Hsapiens.v86, keys = as.character(factor_centrality$ind),
                                     keytype = "UNIPROTID",
                                     columns = "GENENAME")
factor_centrality$names <- genename_df$GENENAME[match(factor_centrality$ind, genename_df$UNIPROTID)]

factor_centrality$names[factor_centrality$rank > 30] <- ""

centrality_plots<-ggplot(factor_centrality, aes(x = rank, y = values, label = names)) +
  geom_point(color = "darkred") +
  geom_text_repel(colour = "black", force = 15) +
  scale_colour_gradientn(colours = pal) +
  labs(x = "Rank in network", y = "Centrality (PageRank)",
       title = "Centrality of Network describing combined drug-agnostic changes") +
  cowplot::theme_cowplot() +
  theme(
    plot.title = element_text(size = 12, face = "bold"),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 1)
  ) +
  grids(linetype = "dashed")

pdf(file = "./paper/Figures/drug_agnostic_centrality.pdf", width = 6)
centrality_plots
dev.off()




union2<-function(g1, g2){

  #Internal function that cleans the names of a given attribute
  CleanNames <- function(g, target){
    #get target names
    gNames <- parse(text = (paste0(target,"_attr_names(g)"))) %>% eval
    #find names that have a "_1" or "_2" at the end
    AttrNeedsCleaning <- grepl("(_\\d)$", gNames )
    #remove the _x ending
    StemName <- gsub("(_\\d)$", "", gNames)

    NewnNames <- unique(StemName[AttrNeedsCleaning])
    #replace attribute name for all attributes
    for( i in NewnNames){

      attr1 <- parse(text = (paste0(target,"_attr(g,'", paste0(i, "_1"),"')"))) %>% eval
      attr2 <- parse(text = (paste0(target,"_attr(g,'", paste0(i, "_2"),"')"))) %>% eval

      g <- parse(text = (paste0("set_",target,"_attr(g, i, value = ifelse(is.na(attr1), attr2, attr1))"))) %>%
        eval

      g <- parse(text = (paste0("delete_",target,"_attr(g,'", paste0(i, "_1"),"')"))) %>% eval
      g <- parse(text = (paste0("delete_",target,"_attr(g,'", paste0(i, "_2"),"')"))) %>% eval

    }

    return(g)
  }


  g <- igraph::union(g1, g2)
  V(g)$consensus_direction <- paste0(V(g)$direction_1, "__", V(g)$direction_2)
  V(g)$consensus_factor <- paste0(V(g)$factor_1, "__", V(g)$factor_2, "__", V(g)$factor_3)

  #loop through each attribute type in the graph and clean
  for(i in c("graph", "edge", "vertex")){
    g <- CleanNames(g, i)
  }

  return(g)
}

# Merge graphs into a super graph
union_factor_graphs <- list()

for (factor in factorS) {
  up_graph <- largest_component(factor_graphs[[factor]][["up"]])
  down_graph <- largest_component(factor_graphs[[factor]][["down"]])
  V(up_graph)$direction<-"up"
  V(down_graph)$direction<-"down"

  # Merge up and down graphs for each factor
  union_factor_graphs[[factor]] <- union2(up_graph, down_graph)
}


g<-union_factor_graphs[[factor_to_vis]]

g1<-union_factor_graphs[["Factor1"]]
g2<-union_factor_graphs[["Factor2"]]
g3<-union_factor_graphs[["Factor3"]]

# Open PDF to save plots
pdf(file = "paper/Supplementary_plots/big_networks.pdf", width = 10, height = 10)

# Plot each factor network with gene names
ggraph(g1, layout = "fr") +
  geom_edge_link(aes(edge_alpha = 0.3), color = "grey") +
  geom_node_point(aes(color = consensus_direction), size = 3) +
  geom_node_text(aes(label = Gene_name), repel = TRUE, size = 5) +
  labs(title = "Network Visualization of Factor1",
       subtitle = "Key genes and interactions associated with Factor1") +
  cowplot::theme_cowplot()+
  theme(
    axis.text = element_blank(),  # Removes axis text
    axis.title = element_blank(),  # Removes axis titles
    axis.ticks = element_blank(),  # Removes axis ticks
    axis.line = element_blank(),
    legend.position = "none"
  )

ggraph(g2, layout = "fr") +
  geom_edge_link(aes(edge_alpha = 0.3), color = "grey") +
  geom_node_point(aes(color = consensus_direction), size = 3) +
  geom_node_text(aes(label = Gene_name), repel = TRUE, size = 4) +
  labs(title = "Network Visualization of Factor2",
       subtitle = "Key genes and interactions associated with Factor2") +
  cowplot::theme_cowplot()+
  theme(
    axis.text = element_blank(),  # Removes axis text
    axis.title = element_blank(),  # Removes axis titles
    axis.ticks = element_blank(),  # Removes axis ticks
    axis.line = element_blank(),
    legend.position = "none"
  )

ggraph(g3, layout = "fr") +
  geom_edge_link(aes(edge_alpha = 0.3), color = "grey") +
  geom_node_point(size = 3, aes(color = consensus_direction)) +
  geom_node_text(aes(label = Gene_name), repel = TRUE, size = 5) +
  labs(title = "Network Visualization of Factor3",
       subtitle = "Key genes and interactions associated with Factor3") +
  cowplot::theme_cowplot()+
  theme(
    axis.text = element_blank(),  # Removes axis text
    axis.title = element_blank(),  # Removes axis titles
    axis.ticks = element_blank(),  # Removes axis ticks
    axis.line = element_blank(),
    legend.position = "none"
  )

# Close PDF device
dev.off()


pdf(file = "paper/Supplementary_plots/big_networks_centrality.pdf", width = 6, height = 6)

# Load required libraries
library(igraph)
library(ggplot2)
library(dplyr)
library(ggrepel) # for non-overlapping text labels

# Define function to compute PageRank centrality and plot it
plot_pagerank_centrality <- function(graph, graph_name) {
  # Compute PageRank centrality
  pagerank_values <- page.rank(graph)$vector

  # Create a data frame with centrality values and ranks
  df <- data.frame(Node = V(graph)$Gene_name, Centrality = pagerank_values, direction = V(graph)$consensus_direction) %>%
    arrange(desc(Centrality)) %>%
    mutate(Rank = row_number())
  # Replace values in the 'direction' column
  df <- df %>%
    mutate(direction = case_when(
      direction == "up__NA" ~ "upregulated network",
      direction == "NA__down" ~ "downregulated network",
      TRUE ~ direction  # Keep other values as they are
    ))

  # Plot using ggplot2
  ggplot(df, aes(x = Rank, y = Centrality)) +
    geom_point(color = "darkred") +
    geom_text_repel(data = df %>% dplyr::filter(Rank <= 50), # Filter for top 50 nodes
                    aes(label = Node),
                    size = 3,
                    box.padding = 0.3,
                    max.overlaps = 50) + # Limits label overlaps to 50 labels
    labs(title = paste("PageRank Centrality vs. Rank for", graph_name),
         x = "Rank",
         y = "PageRank Centrality") +
    facet_wrap(~direction) +
    cowplot::theme_cowplot() +
    theme(
      plot.title = element_text(size = 15, face = "bold"),
      panel.border = element_rect(colour = "black", fill = NA, linewidth = 1),
      axis.text.x = element_text(angle = 90, hjust = 1) # Rotate x-axis labels
    ) +
    grids(linetype = "dashed")
}
# Now call the function for each graph separately
plot_pagerank_centrality(g1, "Factor1")
plot_pagerank_centrality(g2, "Factor2")
plot_pagerank_centrality(g3, "Factor3")

dev.off()

library(ggvenn)
ggvenn(
  list(Factor1=V(factor_graphs[["Factor1"]][["up"]])$name,
       Factor2=V(factor_graphs[["Factor2"]][["up"]])$name,
       Factor3=V(factor_graphs[["Factor3"]][["down"]])$name),
  fill_color = c("#0073C2FF", "#EFC000FF", "#868686FF", "#CD534CFF"),
  stroke_size = 0.5, set_name_size = 4
)
