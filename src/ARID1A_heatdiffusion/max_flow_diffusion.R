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


pdf(file = paste0("./results/heatdiffusion/heatdiffusion_networks.pdf"),
    width = 20, height = 20)

levels <- c("protein_Factor3", "mRNA_Factor3")

# Define a function that returns a list of geom_point layers for different levels

add_levels <- function(df, x, y, size = 6, alpha = 1, stroke = 1, size_increment = 4, ...) {
  # Initialize a size multiplier
  current_size <- size

  # Create a list to store the layers
  layers <- list()

  # Loop through each level and create a geom_point layer
  for (level in levels) {
    layer <- geom_point(
      data = df,
      aes_string(x = x, y = y, colour = level),
      alpha = alpha,
      size = current_size,
      stroke = stroke,
      ... # Pass additional arguments to geom_point
    )

    # Add the layer to the list
    layers <- append(layers, list(layer))

    # Increase the size for the next level
    current_size <- current_size - size_increment
  }

  return(layers)
}

ggraph(subnet, layout = subnet_l) +
  add_levels(df = L_df,
             x = "x",
             y = "y",
             alpha = 1,
             stroke = 5,
             levels = levels,
             size = 10) + # Now this works in the pipeline
  geom_edge_link(start_cap = circle(3, 'mm'),
                 end_cap = circle(3, 'mm'),
                 aes(alpha = weight)) +
  geom_node_point(size = 5) +
  coord_fixed() + theme_void() +
  geom_node_label(aes(label = Gene_name), size = 5, repel = FALSE) +
  theme(legend.position = "bottom") + # Placing legend at the bottom
  scale_colour_gradientn(colours = pal, na.value = "lightgrey")


library(patchwork) # for plot_spacer()

create_custom_legend <- function(df, levels,
                                 size = 6, size_increment = 3, alpha = 1, stroke = 1, title = "Custom Legend", ...) {
  # Create a simple graph with a single node
  graph <- make_empty_graph(1)  # Create a graph with one node

  # Create a layout for the single node (centered)
  layout <- matrix(c(0, 0), ncol = 2)

  # Create a data frame for the single node with an empty placeholder
  node_df <- data.frame(x = 0, y = 0)

  # Add columns to node_df for each level, using the first non-NA value from df for each level
  for (level in levels) {
    # Get the first non-NA value from df for the corresponding level
    node_df[[level]] <- df[[level]][!is.na(df[[level]])][1]
  }

  # Initialize the ggraph object with the single node
  g <- ggraph(graph, layout = layout) +
    theme_void() +  # Remove all background, gridlines, etc.
    coord_fixed()   # Ensure equal aspect ratio

  # Initialize a size multiplier
  current_size <- size

  # Loop through each level to add corresponding geom_point layers
  for (level in levels) {
    g <- g +
      geom_point(
        data = node_df,
        aes_string(x = "x", y = "y", colour = level),
        size = current_size,
        stroke = stroke,
        alpha = alpha,
        ...
      )
    # Decrease the size for the next level (or adjust as needed)
    current_size <- current_size - size_increment
  }

  # Add a title and a square around the node
  g <- g +
    ggtitle(title) +  # Add the title above the legend
    theme(plot.title = element_text(hjust = 0.5),
          legend.position = "none") +  # Center the title
    geom_rect(aes(xmin = -1, xmax = 1, ymin = -1, ymax = 1),
              fill = NA, colour = "black", size = 0.5)  # Draw a square around the node

  return(g)
}

# Function to create a combined legend grid
create_combined_legend <- function(L_df, levels, colours, nrow = 1, ncol = 2) {
  # Initialize an empty list to store legends
  legend_list <- list()

  # Loop through each level to create individual plots and extract legends
  for (i in seq_along(levels)) {
    level <- levels[i]

    # Filter out rows with missing values in the current level
    level_data <- L_df[!is.na(L_df[[level]]), ]

    # Create the ggplot for the current level
    p <- ggplot(data = level_data) +
      geom_point(aes_string(x = level, y = level, fill = level), shape = 21) +
      scale_fill_gradientn(
        name = level,
        limits = c(0, max(level_data[[level]], na.rm = TRUE)),
        colours = colours,
        values = c(0, 0.5, 1)
      )

    # Extract the legend from the plot (handling multiple components)
    legend_list[[i]] <- cowplot::get_legend(p)
  }
  # Create a blank plot to align legends
  blank_p <- plot_spacer() + theme_void()

  # Combine all legends into a grid
  combined_legends <- plot_grid(
    plotlist = legend_list,
    blank_p,
    nrow = nrow,
    ncol = ncol
  )

  # Return the combined legends
  return(combined_legends)
}

# Example usage
levels <- c("mRNA_Factor3", "protein_Factor3")
colours <- rev(c("#000000", "#FFFFFF", "#BA0000"))

combined_legends <- create_combined_legend(L_df = L_df, levels = levels, colours = pal, nrow = 3, ncol = 1)

plot_list<-list(create_custom_legend(L_df, levels = levels, size = 50, size_increment = 20, title = "Graph Legend")+ # Placing legend at the bottom
                  scale_colour_gradientn(colours = pal, na.value = "lightgrey"),
                combined_legends)

plot_grid(
  plotlist = plot_list,
  nrow = 1
)

ggraph(subnet, layout = subnet_l) +
  geom_point(data = L_df[L_df$Gene_name != "Centroid_All_Factors",],
             aes(x = x, y = y, colour = arid1a_up_probabilities),
             alpha = 1, size = 10,
             stroke = 5) +
  geom_edge_link(start_cap = circle(3, 'mm'),
                 end_cap = circle(3, 'mm'),
                 aes(alpha = weight)) +
  geom_node_point(size = 5) +
  coord_fixed() + theme_void() +
  geom_node_label(aes(label = Gene_name),size=4, repel = FALSE) +
  ggtitle("Heat diffusion of receptors up in ARID1A ko") +
  theme(legend.position = "bottom")  # Placing legend at the bottom

ggraph(subnet, layout = subnet_l) +
  geom_point(data = L_df[L_df$Gene_name != "Centroid_All_Factors",],
             aes(x = x, y = y, colour = combined_up_probabilities),
             alpha = 1, size = 10,
             stroke = 5) +
  geom_edge_link(start_cap = circle(3, 'mm'),
                 end_cap = circle(3, 'mm'),
                 aes(alpha = weight)) +
  geom_node_point(size = 5) +
  coord_fixed() + theme_void() +
  geom_node_label(aes(label = Gene_name),size=4, repel = FALSE) +
  ggtitle("Heat diffusion of receptors up in combination drugs") +
  theme(legend.position = "bottom")  # Placing legend at the bottom

ggraph(subnet, layout = subnet_l) +
  geom_point(data = L_df[L_df$Gene_name != "Centroid_All_Factors",],
             aes(x = x, y = y, colour = arid1a_down_probabilities),
             alpha = 1, size = 10,
             stroke = 5) +
  geom_edge_link(start_cap = circle(3, 'mm'),
                 end_cap = circle(3, 'mm'),
                 aes(alpha = weight)) +
  geom_node_point(size = 5) +
  coord_fixed() + theme_void() +
  geom_node_label(aes(label = Gene_name),size=4, repel = FALSE) +
  ggtitle("Heat diffusion of receptors down in ARID1A ko") +
  theme(legend.position = "bottom")  # Placing legend at the bottom

ggraph(subnet, layout = subnet_l) +
  geom_point(data = L_df[L_df$Gene_name != "Centroid_All_Factors",],
             aes(x = x, y = y, colour = combined_down_probabilities),
             alpha = 1, size = 10,
             stroke = 5) +
  geom_edge_link(start_cap = circle(3, 'mm'),
                 end_cap = circle(3, 'mm'),
                 aes(alpha = weight)) +
  geom_node_point(size = 5) +
  coord_fixed() + theme_void() +
  geom_node_label(aes(label = Gene_name),size=4, repel = FALSE) +
  ggtitle("Heat diffusion of receptors down in combination drugs") +
  theme(legend.position = "bottom")  # Placing legend at the bottom
dev.off()

# basic scatterplot
L_df$is_receptor = L_df$Gene_name %in% receptors

# Define the function
process_data_for_labeling <- function(data, x_col, y_col, label_col,
                                      threshold_x = NULL, threshold_y = NULL,
                                      quantile_threshold = 0.9) {

  # Determine thresholds if not provided
  if (is.null(threshold_x)) {
    threshold_x <- quantile(data[[x_col]], quantile_threshold)
  }
  if (is.null(threshold_y)) {
    threshold_y <- quantile(data[[y_col]], quantile_threshold)
  }

  # Create a new column to flag high values
  data$label_flag <- with(data, ifelse(
    data[[x_col]] > threshold_x |
      data[[y_col]] > threshold_y,
    data[[label_col]],
    NA
  ))

  return(data)
}

pdf(file = paste0("./results/heatdiffusion/heatdiffusion_scatters.pdf"),
    width = 10, height = 10)

# Usage example
L_df_processed <- process_data_for_labeling(L_df, "arid1a_up_probabilities", "combined_up_probabilities", "Gene_name")
ggplot(L_df_processed, aes(x=arid1a_up_probabilities, y=combined_up_probabilities, label=label_flag, color = is_receptor)) +
  geom_point() +
  ggtitle("Comparison of 'heat' of nodes in network : up receptors vs up") +
  ggrepel::geom_text_repel(size = 8, colour = "black") + cowplot::theme_cowplot()+ scale_color_manual(values=wes_palette("Chevalier1"))

# basic scatterplot
L_df_processed <- process_data_for_labeling(L_df, "arid1a_down_probabilities", "combined_down_probabilities", "Gene_name")
ggplot(L_df_processed, aes(x=arid1a_down_probabilities, y=combined_down_probabilities, label=label_flag, color = is_receptor)) +
  geom_point() +
  ggtitle("Comparison of 'heat' of nodes in network : down receptors vs down") +
  ggrepel::geom_text_repel(size = 8, colour = "black") + cowplot::theme_cowplot()+ scale_color_manual(values=wes_palette("Chevalier1"))

# basic scatterplot
L_df_processed <- process_data_for_labeling(L_df, "arid1a_up_probabilities", "combined_down_probabilities", "Gene_name")
ggplot(L_df_processed, aes(x=arid1a_up_probabilities, y=combined_down_probabilities, label=label_flag, color = is_receptor)) +
  geom_point() +
  ggtitle("Comparison of 'heat' of nodes in network : up receptors vs down") +
  ggrepel::geom_text_repel(size = 8, colour = "black") + cowplot::theme_cowplot()+ scale_color_manual(values=wes_palette("Chevalier1"))

# basic scatterplot
L_df_processed <- process_data_for_labeling(L_df, "arid1a_down_probabilities", "combined_up_probabilities", "Gene_name")
ggplot(L_df_processed, aes(x=arid1a_down_probabilities, y=combined_up_probabilities, label=label_flag, color = is_receptor)) +
  geom_point() +
  ggtitle("Comparison of 'heat' of nodes in network : down receptors vs up") +
  ggrepel::geom_text_repel(size = 8, colour = "black") + cowplot::theme_cowplot() + scale_color_manual(values=wes_palette("Chevalier1"))
dev.off()


vector_to_binary_matrix <- function(vec) {
  # Find indices of non-zero elements
  nonzero_indices <- which(vec != 0)

  # Number of non-zero elements
  num_nonzero <- length(nonzero_indices)

  # Create an empty matrix
  result_matrix <- matrix(0, nrow = length(vec), ncol = num_nonzero)

  # Populate the matrix with the actual values from vec
  for (i in 1:num_nonzero) {
    result_matrix[nonzero_indices[i], i] <- vec[nonzero_indices[i]]
  }

  # Return the resulting matrix
  return(result_matrix)
}


# Define the function to create a heatmap
create_heatmap <- function(data_frame, to_test, nodes_to_see, subnet, rtks_list, plot_title="") {
  # Ensure that 'L_df.Gene_name' column is present in 'data_frame'
  if (!"L_df.Gene_name" %in% colnames(data_frame)) {
    stop("The data frame must contain the column 'L_df.Gene_name'")
  }
  # Prepare the data frame for heatmap
  heatmap_data <- data_frame
  heatmap_data <- heatmap_data[match(nodes_to_see, heatmap_data$L_df.Gene_name),]
  # Extract the matrix for heatmap
  data_matrix <- as.matrix(heatmap_data[, -1])
  colnames(data_matrix) <- rtks_list
  rownames(data_matrix) <- heatmap_data$L_df.Gene_name
  # Create and return the heatmap
  Heatmap(
    data_matrix,
    name = "HEAT from receptors", na_col = "darkgrey",
    rect_gp = gpar(col = "white", lwd = 5), cluster_columns = F, cluster_rows = F,
    column_gap=unit(.05, "npc"),
    row_title = "Nodes in network", row_title_gp = gpar(fontsize = 10),
    column_title = plot_title, column_title_gp = gpar(fontsize = 10),
    col = viridis::magma(256)
  )
}


# pdf(file = paste0("./results/heatdiffusion/individual_receptor_heat.pdf"),
#     width = 9, height = 3)
#
# library(ComplexHeatmap)
# library(viridis)
#
# # Example usage (you need to provide the actual objects for `arid1a_up_receptors`, `subnet`, and `vector_to_binary_matrix`):
# to_test<-arid1a_up_receptors
# df_arid1a_up<-data.frame(L_df$Gene_name, perform_random_walk(subnet, vector_to_binary_matrix(to_test)))
# rtks_list<-L_df$Gene_name[as.logical(rowSums(vector_to_binary_matrix(to_test)))]
# arid1a_up_heatmap <- create_heatmap(data_frame = df_arid1a_up,
#                           to_test = arid1a_up_receptors,
#                           nodes_to_see = c("JUN", "PRKD1", "MAPK1", "MAPK3"),
#                           subnet = subnet,
#                           rtks_list = rtks_list,
#                           plot_title="Receptors upregulated after ARID1A KO")
# to_test<-combined_up_receptors
# df_combined_up<-data.frame(L_df$Gene_name, perform_random_walk(subnet, vector_to_binary_matrix(to_test)))
# rtks_list<-L_df$Gene_name[as.logical(rowSums(vector_to_binary_matrix(to_test)))]
# combined_heatmap <- create_heatmap(data_frame = df_combined_up,
#                           to_test = combined_up_receptors,
#                           nodes_to_see = c("JUN", "PRKD1","MAPK1", "MAPK3"),
#                           subnet = subnet,
#                           rtks_list = rtks_list,
#                           plot_title="Receptors upregulated after combined therapy")
#
# ht_list = arid1a_up_heatmap + combined_heatmap
# draw(ht_list)
#
# to_test<-arid1a_down_receptors
# heatmap_arid1a_up<-data.frame(L_df$Gene_name, perform_random_walk(subnet, vector_to_binary_matrix(to_test)))
# rtks_list<-L_df$Gene_name[as.logical(rowSums(vector_to_binary_matrix(to_test)))]
# heatmap <- create_heatmap(data_frame = heatmap_arid1a_up,
#                           to_test = arid1a_up_receptors,
#                           nodes_to_see = c("SPRED1", "NGFR","MAPK1", "MAPK3"),
#                           subnet = subnet,
#                           rtks_list = rtks_list)
#
#
# dev.off()

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

# Now you can simply call this function with your graph




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

