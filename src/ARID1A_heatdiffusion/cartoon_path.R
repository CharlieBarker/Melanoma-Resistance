

maxflowSubgraph<-make_shortest_paths_plot(
  sources = c("EGFR", "ROS1", "FGFR1"),
  sinks =  c("JUN"),
  g = subnet,
  conv_nodes = conv_nodes,
  factor_genes_names_df = factor_genes_names_df,
  factor_weights_wide = factor_weights_wide,
  pal = pal,return_df = T
)

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
# Define node types (manual assignment for demonstration)
node_types <- data.frame(
  name = c("ROS1", "PRKD1", "JUN", "EGF", "EGFR", "FGFR1", "PTPN6"),
  type = c("receptor", "kinase", "transcription_factor", "receptor", "receptor", "receptor", "kinase")
)

# Split the edges into two columns: from and to
edges_split <- do.call(rbind, strsplit(edges_to_keep, "--"))
colnames(edges_split) <- c("from", "to")

# Create an igraph object from the edge list
subgraph <- graph_from_data_frame(edges_split, directed = F)



# Add node types as a vertex attribute
V(subgraph)$type <- node_types$type[match(V(subgraph)$name, node_types$name)]

# Function to create a layout for nodes based on their type
create_custom_layout <- function(graph, type_labels, y_position) {
  type_nodes <- V(graph)$name[V(graph)$type %in% type_labels]
  subgraph_type <- induced_subgraph(graph, type_nodes)
  layout_type <- layout_with_sugiyama(subgraph_type)$layout
  
  # Adjust layout: y_position determines where the layout will be placed vertically
  layout_type <- as.data.frame(layout_type)
  colnames(layout_type) <- c("x", "y")
  layout_type$y <- layout_type$y + y_position
  rownames(layout_type) <- type_nodes
  
  return(layout_type)
}

# Create layouts for different types of nodes
layout_receptors <- create_custom_layout(subgraph, "receptor", y_position = 1)
layout_tfs <- create_custom_layout(subgraph, "transcription_factor", y_position = -1)
layout_others <- create_custom_layout(subgraph, "kinase", y_position = 0)

# Combine layouts
all_layouts <- rbind(layout_receptors, layout_tfs, layout_others)
all_layouts<-all_layouts[V(subgraph)$name,]

pdf(file = paste0("./results/heatdiffusion/paths_cartoon.pdf"), 
    width = 9, height = 6)

# Plot the graph with custom layout
ggraph(subgraph, layout = all_layouts) + 
  geom_edge_link(aes(color = "line"), width = 2) + # Customize the edge color and thickness
  geom_node_point(color = "black", size = 8, shape = 21, fill = "white", stroke = 2) + # Customize the node color, outline, and size
  geom_node_text(aes(label = name), nudge_y = 0.5, size = 4, fontface = "bold", color = "black") + # Avoid label overlap
  scale_edge_color_manual(values = "#0098D4") + # Set edge color
  theme_void() +
  theme(legend.position = "none") + # Remove legend
  coord_fixed() +   
  scale_x_continuous(expand = expansion(c(1, 1))) +
  scale_y_continuous(expand = expansion(c(1, 1)))


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

# Step 2: Manually assign node types
# This is an example assignment based on the known biological roles of these proteins.
node_types <- data.frame(
  name = nodes,
  type = c("receptor",            # ITGA4
           "kinase",              # PRKACA
           "kinase",              # RPS6KA3
           "transcription_factor",# FOS
           "kinase",              # MAPK1
           "kinase",              # MAPK3
           "transcription_factor" # ETS2
  )
)
# Split the edges into two columns: from and to
edges_split <- do.call(rbind, strsplit(edges_to_keep, "--"))
colnames(edges_split) <- c("from", "to")
# Create an igraph object from the edge list
subgraph <- graph_from_data_frame(edges_split, directed = F)
# Add node types as a vertex attribute
V(subgraph)$type <- node_types$type[match(V(subgraph)$name, node_types$name)]
# Create layouts for different types of nodes
layout_receptors <- create_custom_layout(subgraph, "receptor", y_position = 2)
layout_tfs <- create_custom_layout(subgraph, "transcription_factor", y_position = -1)
layout_others <- create_custom_layout(subgraph, "kinase", y_position = 0)
# Combine layouts
all_layouts <- rbind(layout_receptors, layout_tfs, layout_others)
all_layouts<-all_layouts[V(subgraph)$name,]
# Plot the graph with custom layout
ggraph(subgraph, layout = all_layouts) + 
  geom_edge_link(aes(color = "line"), width = 2) + # Customize the edge color and thickness
  geom_node_point(color = "black", size = 8, shape = 21, fill = "white", stroke = 2) + # Customize the node color, outline, and size
  geom_node_text(aes(label = name), nudge_y = 0.5, size = 4, fontface = "bold", color = "black") + # Avoid label overlap
  scale_edge_color_manual(values = "red") + # Set edge color
  theme_void() +
  theme(legend.position = "none") + # Remove legend
  coord_fixed() +   
  scale_x_continuous(expand = expansion(c(1, 1))) +
  scale_y_continuous(expand = expansion(c(1, 1)))

dev.off()