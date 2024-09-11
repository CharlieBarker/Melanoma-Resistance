
# Load the necessary libraries
library(igraph)
library(ggraph)
library(ggplot2)
library(dplyr)

# Set environment and working directory
packLib <- "/usr/lib/R"
setwd(dir = "~/Desktop/Melanoma_Resistance/")
if (file.exists(packLib)) {
  reticulate::use_condaenv("py37", required = TRUE)
}

source("./src/functions/pathwayLayout.R")

# Define the function to plot the graph with a flow cutoff
plot_graph_with_flow_cutoff <- function(graph, flow_cutoff, print_graph=T) {

  # Filter edges based on flow value
  filtered_graph <- delete_edges(graph, E(graph)[E(graph)$flow < flow_cutoff])

  # Delete vertices that are not connected (no edges)
  filtered_graph <- delete_vertices(filtered_graph, V(filtered_graph)[degree(filtered_graph) == 0])

  # Extract node and edge data for plotting
  nodes <- igraph::as_data_frame(filtered_graph, what = "vertices")
  edges <- igraph::as_data_frame(filtered_graph, what = "edges")

  V(filtered_graph)$name <- V(filtered_graph)$Gene_name

  path_layout <- pathwayLayout(filtered_graph, receptor_y_position = .5)$layout_matrix

  if(print_graph){
    # Plot the graph using ggraph
    ggraph(filtered_graph, layout = path_layout) +
      geom_node_point() +
      geom_edge_hive(aes(alpha = flow),
                     start_cap = circle(3, 'mm'),
                     end_cap = circle(3, 'mm')) +
      coord_fixed() +
      theme_void() +
      geom_node_label(aes(label = name), size = 4, repel = T) +
      theme(legend.position = "bottom") +
      scale_colour_gradientn(colours = pal, na.value = "lightgrey")
  }else{
    return(filtered_graph)
  }
}

# Load necessary packages
library(igraph)
library(ggraph)
library(dplyr)

# Define your edges for communication with JUN
edges_to_keep <- c("ROS1--PRKD1",
                   "PRKD1--JUN",
                   "EGF--JUN",
                   "EGF--EGFR",
                   "EGFR--PRKD1",
                   "EGFR--FGFR1",
                   "ROS1--PTPN6",
                   "EGFR--PTPN6",
                   "PTPN6--PRKD1")

# Split the edges into two columns: from and to
edges_split <- do.call(rbind, strsplit(edges_to_keep, "--"))
colnames(edges_split) <- c("from", "to")

pdf(file = paste0("./results/heatdiffusion/paths_cartoon.pdf"),
    width = 9, height = 6)

maxflowSubgraph<-make_shortest_paths_plot(
  sources = c("EGFR", "ROS1", "FGFR1"),
  sinks =  c("JUN"),
  g = subnet,
  conv_nodes = conv_nodes,
  factor_genes_names_df = factor_genes_names_df,
  factor_weights_wide = factor_weights_wide,
  pal = pal, return_df = T
)
# Example usage with flow cutoff
plot_graph_with_flow_cutoff(graph = maxflowSubgraph, flow_cutoff = 700, print_graph = T)+
  ggtitle("EGFR, ROS1 and FGFR1 to JUN")  # Example flow cutoff

maxflowSubgraph<-make_shortest_paths_plot(
  sources = c("EGFR", "ROS1", "FGFR1"),
  sinks =  c("MAPK1"),
  g = subnet,
  conv_nodes = conv_nodes,
  factor_genes_names_df = factor_genes_names_df,
  factor_weights_wide = factor_weights_wide,
  pal = pal,return_df = T
)
# Example usage with flow cutoff
plot_graph_with_flow_cutoff(graph = maxflowSubgraph, flow_cutoff = 1000)+
  ggtitle("EGFR, ROS1 and FGFR1 to MAPK1")  # Example flow cutoff

maxflowSubgraph<-make_shortest_paths_plot(
  sources = c("ITGA4"),
  sinks =  c("MAPK1"),
  g = subnet,
  conv_nodes = conv_nodes,
  factor_genes_names_df = factor_genes_names_df,
  factor_weights_wide = factor_weights_wide,
  pal = pal,return_df = T
)
# Example usage with flow cutoff
plot_graph_with_flow_cutoff(graph = maxflowSubgraph, flow_cutoff = 100)+
  ggtitle("ITGA4, ROS1 and FGFR1 to MAPK1")  # Example flow cutoff


# Create an igraph object from the edge list
subgraph <- graph_from_data_frame(edges_split, directed = F)
path_layout <- pathwayLayout(subgraph, receptor_y_position = .5)$layout_matrix
# Plot the graph with custom layout
londonUnderground_plot(subgraph, layout = path_layout, color = "#0098D4")

# Define your edges for communication with MAPK1. MAPK3
edges_to_keep <- c("ITGA4--PRKACA",
                   "ITGA4--RPS6KA3",
                   "RPS6KA3--FOS",
                   "RPS6KA3--MAPK1",
                   "RPS6KA3--MAPK3",
                   "FOS--MAPK1",
                   "FOS--PRKACA",
                   "FOS--MAPK3",
                   "ETS2--MAPK1",
                   "ETS2--MAPK3",
                   "ETS2--PRKACA")
# Step 1: Extract unique nodes
nodes <- unique(unlist(strsplit(edges_to_keep, "--")))
# Split the edges into two columns: from and to
edges_split <- do.call(rbind, strsplit(edges_to_keep, "--"))
colnames(edges_split) <- c("from", "to")
# Create an igraph object from the edge list
subgraph <- graph_from_data_frame(edges_split, directed = F)
# Create an igraph object from the edge list
subgraph <- graph_from_data_frame(edges_split, directed = F)
path_layout <- pathwayLayout(subgraph, receptor_y_position = .5)$layout_matrix
# Plot the graph with custom layout
londonUnderground_plot(subgraph, layout = path_layout, color = "#E32017")
dev.off()
