

library(ggraph)
library(igraph)

# Assuming you have already created 'graph' and 'prior' as in your code

to_plot_df <- full_max_flow_out_df %>%
  filter(normalized_flow != 0) %>%
  filter(factor == "Factor1")
to_plot_df$weight<-NULL
to_plot_df<-to_plot_df[to_plot_df$regulation == "ARID1A KO, Increased signalling",]
graph <- as_tbl_graph(to_plot_df) |> 
  mutate(Centrality = centrality_degree(mode = 'in'))



# Find the vertex IDs for the source and sink nodes
source_node <- V(graph)[name == "source"]$name
sink_node <- V(graph)[name == "sink"]$name

# Calculate the shortest paths from source to sink
shortest_paths <- all_shortest_paths(graph, from = source_node, to = sink_node)

# Extract the shortest paths
paths <- shortest_paths$res

library(ggraph)
library(igraph)

pdf(# The directory you want to save the file in
  width = 10, # The width of the plot in inches
  height = 10,
  file = "./paper/plots/max_flow_shortest_path.pdf")
# Plot each shortest path as a separate graph
graph_list<-list()
for (i in 1:length(paths)) {
  # Extract vertex names from the path
  vertex_names <- paths[[i]]$name
  g <- graph_from_data_frame(data.frame(from = vertex_names[-length(vertex_names)], to = vertex_names[-1]))
  
  # Construct a subgraph containing only the vertices in the path
  subgraph <- induced_subgraph(g, vertex_names)
  graph_list[[i]]<-subgraph
}
path1<-ggraph(graph_list[[1]], layout = "linear") + 
  geom_edge_link(aes(color = "line"), width=5) + # Customize the edge color and thickness
  geom_node_point(color = "black", size = 8, shape = 21, fill = "white", stroke = 3) + # Customize the node color, outline, and size
  geom_node_text(aes(label = name),nudge_y = -.3, size = 4, fontface = "bold", color = "black") + # Avoid label overlap
  scale_edge_color_manual(values = "#0098D4") + # Set edge color
  theme_void() +
  theme(legend.position = "none") + # Remove legend
  coord_fixed() +   
  scale_x_continuous(expand = expansion(c(1, 1))) +
  scale_y_continuous(expand = expansion(c(1, 1)))
path2<-ggraph(graph_list[[2]], layout = "linear") + 
  geom_edge_link(aes(color = "line"), width=5) + # Customize the edge color and thickness
  geom_node_point(color = "black", size = 8, shape = 21, fill = "white", stroke = 3) + # Customize the node color, outline, and size
  geom_node_text(aes(label = name),nudge_y = -.3, size = 4, fontface = "bold", color = "black") + # Avoid label overlap
  scale_edge_color_manual(values = "#000000") + # Set edge color
  theme_void() +
  theme(legend.position = "none") + # Remove legend
  coord_fixed() +   
  scale_x_continuous(expand = expansion(c(1, 1))) +
  scale_y_continuous(expand = expansion(c(1, 1)))
path3<-ggraph(graph_list[[3]], layout = "linear") + 
  geom_edge_link(aes(color = "line"), width=5) + # Customize the edge color and thickness
  geom_node_point(color = "black", size = 8, shape = 21, fill = "white", stroke = 3) + # Customize the node color, outline, and size
  geom_node_text(aes(label = name),nudge_y = -.3, size = 4, fontface = "bold", color = "black") + # Avoid label overlap
  scale_edge_color_manual(values = "#E32017") + # Set edge color
  theme_void() +
  theme(legend.position = "none") + # Remove legend
  coord_fixed() +   
  scale_x_continuous(expand = expansion(c(1, 1))) +
  scale_y_continuous(expand = expansion(c(1, 1)))
path4<-ggraph(graph_list[[4]], layout = "linear") + 
  geom_edge_link(aes(color = "line"), width=5) + # Customize the edge color and thickness
  geom_node_point(color = "black", size = 8, shape = 21, fill = "white", stroke = 3) + # Customize the node color, outline, and size
  geom_node_text(aes(label = name),nudge_y = -.3, size = 4, fontface = "bold", color = "black") + # Avoid label overlap
  scale_edge_color_manual(values = "#FFD300") + # Set edge color
  theme_void() +
  theme(legend.position = "none") + # Remove legend
  coord_fixed() +   
  scale_x_continuous(expand = expansion(c(1, 1))) +
  scale_y_continuous(expand = expansion(c(1, 1)))
path5<-ggraph(graph_list[[5]], layout = "linear") + 
  geom_edge_link(aes(color = "line"), width=5) + # Customize the edge color and thickness
  geom_node_point(color = "black", size = 8, shape = 21, fill = "white", stroke = 3) + # Customize the node color, outline, and size
  geom_node_text(aes(label = name),nudge_y = -.3, size = 4, fontface = "bold", color = "black") + # Avoid label overlap
  scale_edge_color_manual(values = "#00782A") + # Set edge color
  theme_void() +
  theme(legend.position = "none") + # Remove legend
  coord_fixed() +   
  scale_x_continuous(expand = expansion(c(1, 1))) +
  scale_y_continuous(expand = expansion(c(1, 1)))
path6<-ggraph(graph_list[[6]], layout = "linear") + 
  geom_edge_link(aes(color = "line"), width=5) + # Customize the edge color and thickness
  geom_node_point(color = "black", size = 8, shape = 21, fill = "white", stroke = 3) + # Customize the node color, outline, and size
  geom_node_text(aes(label = name),nudge_y = -.3, size = 4, fontface = "bold", color = "black") + # Avoid label overlap
  scale_edge_color_manual(values = "#9B0056") + # Set edge color
  theme_void() +
  theme(legend.position = "none") + # Remove legend
  coord_fixed() +   
  scale_x_continuous(expand = expansion(c(1, 1))) +
  scale_y_continuous(expand = expansion(c(1, 1)))
path7<-ggraph(graph_list[[7]], layout = "linear") + 
  geom_edge_link(aes(color = "line"), width=5) + # Customize the edge color and thickness
  geom_node_point(color = "black", size = 8, shape = 21, fill = "white", stroke = 3) + # Customize the node color, outline, and size
  geom_node_text(aes(label = name),nudge_y = -.3, size = 4, fontface = "bold", color = "black") + # Avoid label overlap
  scale_edge_color_manual(values = "#95CDBA") + # Set edge color
  theme_void() +
  theme(legend.position = "none") + # Remove legend
  coord_fixed() +   
  scale_x_continuous(expand = expansion(c(1, 1))) +
  scale_y_continuous(expand = expansion(c(1, 1)))
ggarrange(path1, 
          path2,
          path3,
          path4,
          path5,
          path6,
          path7,
          ncol = 1,  common.legend = T)

dev.off()

# Make the union of the graphs in the list
union_graph <- Reduce(union2, graph_list)
test<-weights[weights$node %in% V(union_graph)$name,]
test[order(test$node),]