
# Define the function to classify nodes
classify_nodes <- function(nodes) {
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
  
  #get list of transcription factors
  transcriptionFactors <-
    import_omnipath_annotations(
      resources = 'TFcensus',
      entity_types = 'protein'
    ) %>%
    pull(genesymbol) %>%
    unique
  
  # Create a data frame with nodes and their types
  node_types <- data.frame(
    name = nodes,
    type = ifelse(nodes %in% receptors, "receptor", 
                  ifelse(nodes %in% transcriptionFactors, "transcription_factor", "other"))
  )
  
  return(node_types)
}

# Helper function to create a layout for nodes based on their type
create_custom_layout <- function(graph, type_labels, y_position) {
  # Get the nodes of the specified type
  type_nodes <- V(graph)$name[V(graph)$type %in% type_labels]
  
  # Create a subgraph with only these nodes
  subgraph_type <- induced_subgraph(graph, type_nodes)
  
  # Use the Sugiyama layout for this subgraph
  layout_type <- layout_with_sugiyama(subgraph_type)$layout
  
  # Convert the layout to a data frame and adjust y positions
  layout_type <- as.data.frame(layout_type)
  colnames(layout_type) <- c("x", "y")
  layout_type$y <- layout_type$y + y_position
  rownames(layout_type) <- type_nodes
  
  return(layout_type)
}

# Define the custom pathway layout function
pathwayLayout <- function(graph,
                          receptor_y_position = 1,
                          transcription_factor_y_position = -1,
                          other_y_position = 0) {
  # Extract node names
  nodes <- V(graph)$name
  
  # Classify nodes
  node_types <- classify_nodes(nodes)
  
  # Add 'type' attribute to vertices
  V(graph)$type <- node_types$type[match(V(graph)$name, node_types$name)]
  
  # Create layouts for different types of nodes
  layout_receptors <- create_custom_layout(graph, "receptor", y_position = receptor_y_position)
  layout_tfs <- create_custom_layout(graph, "transcription_factor", y_position = transcription_factor_y_position)
  layout_others <- create_custom_layout(graph, "other", y_position = other_y_position)
  
  # Combine all layouts into one data frame
  all_layouts <- rbind(layout_receptors, layout_tfs, layout_others)
  
  # Reorder the layout to match the original graph's node order
  all_layouts <- all_layouts[V(graph)$name,]
  
  # Convert the data frame to a matrix for igraph layout
  layout_matrix <- as.matrix(all_layouts)
  
  return(layout_matrix)
}
# 
# # Example usage
# # Assuming 'graph' is your igraph object with a 'type' attribute for each vertex
# flow_graph<-plot_graph_with_flow_cutoff(graph = maxflowSubgraph, flow_cutoff = 1000, print_graph = F)
# path_layout <- pathwayLayout(flow_graph)
# 
# 
# # Plot the graph with custom layout
# ggraph(flow_graph, layout = path_layout) + 
#   geom_edge_link(aes(color = "line"), width = 2) + # Customize the edge color and thickness
#   geom_node_point(color = "black", size = 8, shape = 21, fill = "white", stroke = 2) + # Customize the node color, outline, and size
#   geom_node_text(aes(label = name), nudge_y = 0.5, size = 4, fontface = "bold", color = "black") + # Avoid label overlap
#   scale_edge_color_manual(values = "#0098D4") + # Set edge color
#   theme_void() +
#   theme(legend.position = "none") + # Remove legend
#   coord_fixed() +   
#   scale_x_continuous(expand = expansion(c(1, 1))) +
#   scale_y_continuous(expand = expansion(c(1, 1)))

