load('./results/heatdiffusion/data_for_heat_diffusion.Rdata')
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
library(diffusr)
library(OmnipathR)

source("./src/functions/pathwayLayout.R")


# for (variable in V(subnet)$Gene_name) {
#   cat(variable, "\n")
# }

#######DIFFUSION#######

#get list of receptors
receptors <-
  OmnipathR::import_omnipath_intercell(
    parent = 'receptor',
    topology = 'pmtm',
    consensus_percentile = 50,
    loc_consensus_percentile = 30,
    entity_types = 'protein'
  ) %>%
  pull(genesymbol) %>% unique()




process_receptors_lfc<-function(file, receptor_list, g){
  mRNA_lfc<-read.csv(file = file)
  gene_affected <- mRNA_lfc[mRNA_lfc$X %in% V(g)$Gene_name,]
  receptors_affected <- gene_affected[gene_affected$X %in% receptors,]
  receptors_affected <- receptors_affected[receptors_affected$adj.P.Val<0.01,]
  return(receptors_affected)
}

combined_lfc<-process_receptors_lfc(file = "./results/lfc/mRNA/combination_lfc.csv",
                                    receptor_list = receptors,
                                    g= g)

arid1a_lfc<-process_receptors_lfc(file = "./results/lfc/mRNA/arid1a_lfc.csv",
                                  receptor_list = receptors,
                                  g= g)

normalize_vector <- function(vec) {
  # Compute the sum of the vector
  total_sum <- sum(vec)

  # Check if the total_sum is zero to avoid division by zero
  if (total_sum == 0) {
    stop("The sum of the vector elements is zero. Cannot normalize.")
  }

  # Normalize the vector
  normalized_vec <- vec / total_sum

  return(normalized_vec)
}

# Perform a random walk and get stationary distribution
perform_random_walk <- function(graph, start_vector) {
  adj_matrix <- as_adjacency_matrix(graph)
  pt <- random.walk(start_vector, data.matrix(adj_matrix),
                    r = 0.02,
                    correct.for.hubs = T,do.analytical = T)
  return(pt$p.inf)
}


####ARID1A####
all(L_df$Gene_name == V(subnet)$Gene_name)

# Extract and normalize up-regulated receptors
up_receptors <- arid1a_lfc$logFC[match(V(subnet)$Gene_name, arid1a_lfc$Gene)]
up_receptors[is.na(up_receptors) | up_receptors < 0] <- 0
arid1a_up_receptors <- normalize_vector(up_receptors)

# Extract and normalize down-regulated receptors
down_receptors <- arid1a_lfc$logFC[match(V(subnet)$Gene_name, arid1a_lfc$Gene)]
down_receptors[is.na(down_receptors) | down_receptors > 0] <- 0
arid1a_down_receptors <- normalize_vector(down_receptors)

# Perform random walks and get stationary distributions
L_df$arid1a_up_probabilities <- perform_random_walk(subnet, arid1a_up_receptors)
L_df$arid1a_down_probabilities <- perform_random_walk(subnet, arid1a_down_receptors)


####COMBINED####

# Extract and normalize up-regulated receptors
up_receptors <- combined_lfc$logFC[match(V(subnet)$Gene_name, combined_lfc$Gene)]
up_receptors[is.na(up_receptors) | up_receptors < 0] <- 0
combined_up_receptors <- normalize_vector(up_receptors)

# Extract and normalize down-regulated receptors
down_receptors <- combined_lfc$logFC[match(V(subnet)$Gene_name, combined_lfc$Gene)]
down_receptors[is.na(down_receptors) | down_receptors > 0] <- 0
combined_down_receptors <- normalize_vector(down_receptors)

# Perform random walks and get stationary distributions
L_df$combined_up_probabilities <- perform_random_walk(subnet, combined_up_receptors)
L_df$combined_down_probabilities <- perform_random_walk(subnet, combined_down_receptors)



# Example function to position source nodes above
adjust_layout <- function(g, sources_uniprot, original_layout) {
  # Check if the graph and layout have the same number of nodes
  if (nrow(original_layout) != vcount(g)) {
    stop("Layout does not match the number of nodes in the graph.")
  }

  # Get the indices of the source nodes
  source_indices <- which(V(g)$name %in% sources_uniprot)

  # Set a y-coordinate above the main graph for source nodes
  # Assuming original layout is a matrix with two columns (x, y)
  adjusted_layout <- original_layout
  adjusted_layout[source_indices, 2] <- max(original_layout[, 2]) + 1.5  # Set y to be above the max y of the original layout

  # Evenly space the source nodes along the x-axis
  num_sources <- length(source_indices)
  spacing <- (max(original_layout[, 1]) - min(original_layout[, 1])) / (num_sources + 1)
  adjusted_layout[source_indices, 1] <- seq(min(original_layout[, 1]) + spacing, length.out = num_sources, by = spacing)

  return(adjusted_layout)
}

max_flow_subgraph <- function(graph, source, sink, prob = 0.95) {
  # Compute the maximum flow
  mf_result <- max_flow(graph, source = source, target = sink, capacity = E(graph)$weight)

  # Get the flow values on the edges
  edge_flows <- abs(mf_result$flow)

  # Determine the specified percentile of the flow values
  threshold <- quantile(edge_flows, probs = prob)

  # Identify edges with flow values above the threshold
  high_flow_edges <- E(graph)[edge_flows > threshold]

  # Create a subgraph containing only these edges
  subgraph <- induced_subgraph(graph, unique(c(ends(graph, high_flow_edges))))

  # Attach flow values to the edges in the subgraph
  E(subgraph)$flow <- edge_flows[match(as_ids(E(subgraph)), as_ids(E(graph)))]

  # Find connected components
  components <- components(subgraph)

  # Get the component that includes the source and sink
  component_with_source_sink <- which(components$membership[V(subgraph)$name == source] == components$membership[V(subgraph)$name == sink])

  # Extract the nodes in this component
  component_nodes <- V(subgraph)[components$membership == component_with_source_sink]

  # Create the final subgraph based on these nodes
  final_subgraph <- induced_subgraph(subgraph, component_nodes)

  # Attach flow values to the edges in the final subgraph
  E(final_subgraph)$flow <- E(subgraph)$flow[match(E(final_subgraph), E(subgraph))]

  # Return the final subgraph with flow values as edge attributes
  return(final_subgraph)
}


# Define the make_shortest_paths_plot function
make_shortest_paths_plot <- function(sources, sinks, g,
                                     conv_nodes, factor_genes_names_df, factor_weights_wide, pal,
                                     return_df=F) {
  # Match gene names to UniProt IDs
  sources_uniprot <- conv_nodes$uniprt[match(sources, conv_nodes$gene_name)]
  sinks_uniprot <- conv_nodes$uniprt[match(sinks, conv_nodes$gene_name)]
  # Ensure all nodes are in the graph
  if (!all(c(sources_uniprot, sinks_uniprot) %in% V(g)$name)) {
    stop("Some nodes are not present in the graph.")
  }
  # Initialize a list to hold subgraphs of shortest paths
  shortest_paths_list <- list()
  # Function to get edges from path nodes
  edge_from_path <- function(path_nodes) {
    edges <- c()
    for (i in 1:(length(path_nodes) - 1)) {
      edges <- c(edges, which(V(g)$name == path_nodes[i]), which(V(g)$name == path_nodes[i + 1]))
    }
    return(matrix(edges, ncol = 2, byrow = TRUE))
  }
  super_graph <- make_empty_graph(directed = FALSE)
  # Compute shortest paths and create subgraphs
  for (source in sources_uniprot) {
    for (sink in sinks_uniprot) {
      # Compute the shortest path
      sp <- max_flow_subgraph(graph = g, source = source, sink = sink, prob = 0.9)
      # Combine subgraphs
      if (is.null(super_graph) || length(V(super_graph)) == 0) {
        super_graph <- sp
      } else {
        super_graph <- union2(super_graph, sp)
      }
    }
  }

  # Compute the flow passing through each node by summing the flow of incident edges
  V(super_graph)$flow <- sapply(V(super_graph), function(v) {
    # Sum the flow values for all edges incident to this node
    incident_edges <- incident(super_graph, v)
    sum(E(super_graph)[incident_edges]$flow, na.rm = TRUE)
  })
  # Compute the flow passing through each node by summing the flow of incident edges
  V(super_graph)$flow <- sapply(V(super_graph), function(v) {
    # Sum the flow values for all edges incident to this node
    incident_edges <- incident(super_graph, v)
    sum(E(super_graph)[incident_edges]$flow, na.rm = TRUE)
  })

  # Display the list of subgraphs
  test <- super_graph
  subnet_l <- tryCatch({
    igraph::layout_with_sugiyama(test)$layout
  }, error = function(e) {
    stop("Layout calculation failed: ", e$message)
  })

  subnet_l <- adjust_layout(test, sources_uniprot, subnet_l)
  V(test)$Gene_name <- conv_nodes$gene_name[match(V(test)$name, conv_nodes$uniprt)]
  L_df <- as.data.frame(subnet_l)
  colnames(L_df) <- c("x", "y")
  L_df$Gene_name <- V(test)$Gene_name
  L_df$Factor1 <- L_df$Gene_name %in% factor_genes_names_df$values[factor_genes_names_df$factor == "Factor1"]
  L_df$Factor2 <- L_df$Gene_name %in% factor_genes_names_df$values[factor_genes_names_df$factor == "Factor2"]
  L_df$Factor3 <- L_df$Gene_name %in% factor_genes_names_df$values[factor_genes_names_df$factor == "Factor3"]

  # Add columns for each factor and view
  L_df <- L_df %>%
    add_factor_view_columns(factor_weights_wide, "Factor1") %>%
    add_factor_view_columns(factor_weights_wide, "Factor2") %>%
    add_factor_view_columns(factor_weights_wide, "Factor3")

  out_ggplot <- ggraph(test, layout = subnet_l) +
    geom_point(data = L_df[L_df$Gene_name != "Centroid_All_Factors",],
               aes(x = x, y = y, colour = protein_Factor3),
               alpha = 1, size = 10,
               stroke = 5) +
    geom_point(data = L_df[L_df$Gene_name != "Centroid_All_Factors",],
               aes(x = x, y = y, colour = mRNA_Factor3),
               alpha = 1, size = 6,
               stroke = 5) +
    geom_edge_hive(start_cap = circle(3, 'mm'),
                   end_cap = circle(3, 'mm'),
                   aes(alpha = flow)) +
    coord_fixed() + theme_void() +
    geom_node_label(aes(label = Gene_name), size = 4, repel = FALSE) +
    theme(legend.position = "bottom") + # Placing legend at the bottom
    scale_colour_gradientn(colours = pal, na.value = "lightgrey")
  if (return_df) {
    return(test)
  }else{
    return(out_ggplot)
  }
}

# Function to dynamically generate titles
create_title <- function(sources, sinks) {
  sources_str <- paste(sources, collapse = ", ")
  sinks_str <- paste(sinks, collapse = ", ")
  return(paste("Shortest Signaling Paths from", sources_str, "to", sinks_str))
}

# Function to dynamically generate the flow plot title
create_flow_title <- function(sources, sinks) {
  sources_str <- paste(sources, collapse = ", ")
  sinks_str <- paste(sinks, collapse = ", ")
  return(paste("Flow Through Nodes in Signaling Paths from", sources_str, "to", sinks_str))
}

# Function to generate the shortest paths plot
plot_shortest_paths <- function(sources, sinks, g, conv_nodes, factor_genes_names_df, factor_weights_wide, pal) {
  title <- create_title(sources, sinks)
  p1 <- make_shortest_paths_plot(
    sources = sources,
    sinks = sinks,
    g = g,
    conv_nodes = conv_nodes,
    factor_genes_names_df = factor_genes_names_df,
    factor_weights_wide = factor_weights_wide,
    pal = pal
  ) +
    ggtitle(title)
  return(p1)
}

# Function to generate the flow plot
plot_flow <- function(sources, sinks, short_path_graph) {
  flow_df <- data.frame(name = V(short_path_graph)$Gene_name, flow = V(short_path_graph)$flow)
  title <- create_flow_title(sources, sinks)

  p2 <- ggplot(data = flow_df, aes(x = reorder(name, -flow), y = flow)) +  # Reorder by flow descending
    geom_bar(stat = "identity", fill = "steelblue", color = "black") +   # Custom color for the bars
    labs(title = title,
         x = "Gene Name",
         y = "Flow") +    # Add labels
    geom_text(aes(label = round(flow, 2)), vjust = -0.5, size = 3.5)  +
    cowplot::theme_cowplot()  +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))          # Rotate x-axis labels
  return(p2)
}


create_combined_plots <- function(sources, sinks, g, conv_nodes, factor_genes_names_df, factor_weights_wide, pal) {

  # Generate the shortest path graph data
  short_path_graph <- make_shortest_paths_plot(
    sources = sources,
    sinks = sinks,
    g = g,
    conv_nodes = conv_nodes,
    factor_genes_names_df = factor_genes_names_df,
    factor_weights_wide = factor_weights_wide,
    pal = pal,
    return_df = TRUE
  )

  # Generate both plots
  p1 <- plot_shortest_paths(sources, sinks, g, conv_nodes, factor_genes_names_df, factor_weights_wide, pal)
  p2 <- plot_flow(sources, sinks, short_path_graph)

  # Combine the two plots, each taking half of the height
  combined_plot <- p1 / p2 + plot_layout(heights = c(2, 1))

  return(combined_plot)
}

plot_scaled_graph <- function(graph) {
  # Scaling function to normalize flows within a reasonable range
  scale_to_range <- function(x, min_size, max_size) {
    scaled <- ((x - min(x)) / (max(x) - min(x))) * (max_size - min_size) + min_size
    return(scaled)
  }

  # Get the flow values for nodes and edges
  node_flow <- V(graph)$flow
  edge_flow <- E(graph)$flow

  # Check if flow attributes exist, and handle potential errors
  if (is.null(node_flow)) stop("Node attribute 'flow' is missing")
  if (is.null(edge_flow)) stop("Edge attribute 'flow' is missing")

  # Scale node sizes (e.g., between 10 and 30) and edge widths (e.g., between 4 and 10)
  scaled_node_size <- scale_to_range(node_flow, 10, 30)
  scaled_edge_width <- scale_to_range(edge_flow, 4, 10)

  # Plot the graph with dynamically scaled sizes and widths
  plot(
    graph,
    vertex.size = scaled_node_size,            # Scaled node sizes
    edge.width = scaled_edge_width,            # Scaled edge thickness
    vertex.label = V(graph)$name,              # Label nodes by their 'name' attribute
    edge.label = E(graph)$flow,                # Label edges by their 'flow' values
    vertex.color = "lightblue",                # Set node color
    edge.color = "gray",                       # Set edge color
    vertex.frame.color = "transparent",        # Remove node outlines
    vertex.label.family = "sans",              # Set text to Arial,
    vertex.label.color="black"
  )
}


# Load necessary libraries
library(ggplot2)
library(patchwork)

pdf(file = paste0("./results/heatdiffusion/shortestpath_results.pdf"),
    width = 15, height = 20)

######FOR EGFR etc to JUN######

# Plot shortest paths from specific receptors to transcription factors
sources <- c("EGFR", "ROS1", "FGFR1")
sinks =  c("JUN")


# Call the function to generate and display the combined plot grid
combined_plot <- create_combined_plots(
  sources = sources,
  sinks = sinks,
  g = subnet,
  conv_nodes = conv_nodes,
  factor_genes_names_df = factor_genes_names_df,
  factor_weights_wide = factor_weights_wide,
  pal = pal
)

# Display the combined plot
print(combined_plot)

# Generate the shortest path graph data
short_path_graph <- make_shortest_paths_plot(
  sources = sources,
  sinks = sinks,
  g = subnet,
  conv_nodes = conv_nodes,
  factor_genes_names_df = factor_genes_names_df,
  factor_weights_wide = factor_weights_wide,
  pal = pal,
  return_df = TRUE
)
# List of gene names with high flow
selected_genes <- c("PRKD1", "PTPN6", "EGFR", "JUN",
                    "FYN", "ROS1")

# Find the vertex ids corresponding to the selected gene names
vids <- which(V(short_path_graph)$Gene_name %in% selected_genes)

# Extract the subgraph using the vertex ids
sub_short_path_graph <- subgraph(short_path_graph, vids = vids)
V(sub_short_path_graph)$name <- V(sub_short_path_graph)$Gene_name
path_layout <- pathwayLayout(sub_short_path_graph, receptor_y_position = .3)$layout_matrix
londonUnderground_plot(sub_short_path_graph, layout = path_layout, color = "#E32017")
plot_scaled_graph(sub_short_path_graph)

######FOR EGFR etc to MAPK1######

# Plot shortest paths from specific receptors to transcription factors
sources <- c("EGFR", "ROS1")
sinks =  c("MAPK1")


# Call the function to generate and display the combined plot grid
combined_plot <- create_combined_plots(
  sources = sources,
  sinks = sinks,
  g = subnet,
  conv_nodes = conv_nodes,
  factor_genes_names_df = factor_genes_names_df,
  factor_weights_wide = factor_weights_wide,
  pal = pal
)

# Display the combined plot
print(combined_plot)


# Generate the shortest path graph data
short_path_graph <- make_shortest_paths_plot(
  sources = sources,
  sinks = sinks,
  g = subnet,
  conv_nodes = conv_nodes,
  factor_genes_names_df = factor_genes_names_df,
  factor_weights_wide = factor_weights_wide,
  pal = pal,
  return_df = TRUE
)
# List of gene names with high flow
selected_genes <- c("EGFR",
"PTPN6",
"CRK",
"GSK3B",
"MAPK1",
"FYN",
"FGFR2",
"MAPK3",
"PTK2",
"ROS1", "PRKD1", "DAPK1",
"ERBB3",
"FGFR1",
"HCK",
"SYK",
"PRKCD"
)

# Find the vertex ids corresponding to the selected gene names
vids <- which(V(short_path_graph)$Gene_name %in% selected_genes)

# Extract the subgraph using the vertex ids
sub_short_path_graph <- subgraph(short_path_graph, vids = vids)
V(sub_short_path_graph)$name <- V(sub_short_path_graph)$Gene_name

path_layout<-layout_with_sugiyama(sub_short_path_graph
                                  #,layers = c(1, 2, 1, 4, 4, 2, 5, 2, 5, 3, 4, 1)
                                  )$layout #manually define layout layers
londonUnderground_plot(sub_short_path_graph, layout = path_layout, color = "#FFD300")
plot_scaled_graph(sub_short_path_graph)


######FOR ITGA4 etc to MAPK1######


# Plot shortest paths from ITGA4 and NTRK3 to MAPK1
sources <- c("ITGA4", "NTRK3")
sinks =  c("MAPK1")

# Call the function to generate and display the combined plot grid
combined_plot <- create_combined_plots(
  sources = sources,
  sinks = sinks,
  g = subnet,
  conv_nodes = conv_nodes,
  factor_genes_names_df = factor_genes_names_df,
  factor_weights_wide = factor_weights_wide,
  pal = pal
)

# Display the combined plot
print(combined_plot)

# Generate the shortest path graph data
short_path_graph <- make_shortest_paths_plot(
  sources = sources,
  sinks = sinks,
  g = subnet,
  conv_nodes = conv_nodes,
  factor_genes_names_df = factor_genes_names_df,
  factor_weights_wide = factor_weights_wide,
  pal = pal,
  return_df = TRUE
)
# List of gene names with high flow
selected_genes <- c("MAPK1","ITGA4", "PRKACA","RPS6KA3","ETS2", "MAPK3")

# Find the vertex ids corresponding to the selected gene names
vids <- which(V(short_path_graph)$Gene_name %in% selected_genes)

# Extract the subgraph using the vertex ids
sub_short_path_graph <- subgraph(short_path_graph, vids = vids)
V(sub_short_path_graph)$name <- V(sub_short_path_graph)$Gene_name
df <- classify_nodes(V(sub_short_path_graph)$name)
layers <- ifelse(df$type == "receptor", 1, ifelse(df$type == "other", 2, 3))
path_layout<-layout_with_sugiyama(sub_short_path_graph, layers = layers)$layout
londonUnderground_plot(sub_short_path_graph, layout = path_layout, color = "#0098D4")
plot_scaled_graph(sub_short_path_graph)
dev.off()

