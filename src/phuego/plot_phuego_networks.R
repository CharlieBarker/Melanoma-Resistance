

# Load required libraries
library(readxl)
library(ggplot2)
library(dplyr)
library(ggpubr)
library(ggrepel)
library(wesanderson)
library(igraph)
library(cowplot)
library(purrr)
library(stringr)
library(ghibli)
knitr::opts_chunk$set(dev = "ragg_png")

# Set environment and working directory
packLib <- "/usr/lib/R"
setwd(dir = "~/Desktop/Melanoma_Resistance/")
if (file.exists(packLib)) {
  reticulate::use_condaenv("py37", required = TRUE)
}

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
  #loop through each attribute type in the graph and clean
  for(i in c("graph", "edge", "vertex")){
    g <- CleanNames(g, i)
  }

  return(g)
}


# Merge graphs into a super graph
super_graph_down <- make_empty_graph(directed = TRUE)
super_graph_up <- make_empty_graph(directed = TRUE)
super_graph <- make_empty_graph(directed = TRUE)

conv_nodes<-data.frame(uniprot = V(super_graph)$name,
                       genename= V(super_graph)$Gene_name)
# write.csv(conv_nodes, file = "./data/name_converter_for_network.csv")

for (factor in factorS) {
  up_graph<-factor_graphs[[factor]][["up"]]
  down_graph <- factor_graphs[[factor]][["down"]]

  up_graph <- largest_component(factor_graphs[[factor]][["up"]])
  down_graph <- largest_component(factor_graphs[[factor]][["down"]])

  # Merge up and down graphs for each factor
  combined_graph <- union2(up_graph, down_graph)

  # Add the combined graph to the super graph
  super_graph <- union2(super_graph, combined_graph)
  super_graph_down <- union2(super_graph_down, down_graph)
  super_graph_up <- union2(super_graph_up, up_graph)

}

largest_component_graph <- super_graph_down

library(ggraph)
library(igraph)
library(tidygraph)
random_number <- sample(1:1000, 1)
print(random_number)
# Calculate the degree for each vertex
vertex_degrees <- degree(largest_component_graph)
# Add the degree as a vertex attribute
V(largest_component_graph)$degree <- vertex_degrees
# Assuming `V(largest_component_graph)$factor` is stored in a variable
factors <- V(largest_component_graph)$factor
# Replace the factor levels
factors <- ifelse(factors == "Factor1", "Factor 1 (Drug Agnostic)",
                  ifelse(factors == "Factor2", "Factor 2 (Combination Specific)",
                         ifelse(factors == "Factor3", "Factor 3 (ARID1A rewired)", factors)))
# Assign the updated factors back to the graph
V(largest_component_graph)$factor <- factors

#902
set.seed(997)
l <- igraph::layout.graphopt(largest_component_graph)


big_graph<-ggraph(largest_component_graph, layout = l ) +
  geom_edge_link(alpha=0.3) +
  geom_node_point(aes(colour = factor), size = 2) +
  coord_fixed() +
  theme_graph(base_family = "sans",base_size = 12)+
  theme(legend.position = "none",
        plot.margin = unit(c(0,0,0,0), 'lines'))+
  # ghibli stuff
  scale_colour_ghibli_d("YesterdayMedium", direction = -1) +
  ggtitle("graph showing all the factors")


# Facetted graph with reduced node size
facetted_graph <- ggraph(largest_component_graph, layout = l) +
  geom_edge_link(alpha = 0.3) +
  geom_node_point(aes(colour = factor), size = 2) + # Adjust size here
  coord_fixed() +
  facet_graph(~ factor, switch = "x") +
  theme_graph(base_family = "sans") +
  theme(legend.position = "none",
        plot.margin = unit(c(0, 0, 0, 0), 'lines')) +
  scale_colour_ghibli_d("YesterdayMedium", direction = -1) +
  ggtitle("graph showing all the factors facetted")



pdf("./paper/networks/test.pdf", width = 8, height = 6)
big_graph
facetted_graph
# Get network greyed out
if (exists("subnet")) {
  zoom_in_network<-largest_component_graph
  V(zoom_in_network)$zoom_in <- V(zoom_in_network)$name %in% V(subnet)$name
  # Facetted graph with reduced node size
  zoom_in_network <- ggraph(zoom_in_network, layout = l) +
    geom_edge_link(alpha = 0.3) +
    geom_node_point(aes(colour = zoom_in), size = 2) + # Adjust size here
    coord_fixed() +
    theme_graph(base_family = "sans") +
    theme(legend.position = "none",
          plot.margin = unit(c(0, 0, 0, 0), 'lines')) +
    scale_color_manual(values = c("lightgrey", "red")) +
    ggtitle("graph showing the pcsf from figure 2 highlighted ")
  zoom_in_network

} else {
  message("The variable 'subnet' is not defined.")
}

V(largest_component_graph)$is_subgraph = V(largest_component_graph)$factor == "Factor 1 (Drug Agnostic)"
ggraph(largest_component_graph, layout = l ) +
  geom_edge_link(alpha=0.3) +
  geom_node_point(aes(colour=is_subgraph), size = 2) +
  coord_fixed() +
  theme_graph(base_family = "sans",base_size = 12)+
  theme(legend.position = "none",
        plot.margin = unit(c(0,0,0,0), 'lines')) +
  scale_colour_manual(values=c("lightgrey", "#92bbd9ff")) +
  ggtitle("graph showing factor 1 highlighted")

V(largest_component_graph)$is_subgraph = V(largest_component_graph)$factor == "Factor 2 (Combination Specific)"
ggraph(largest_component_graph, layout = l ) +
  geom_edge_link(alpha=0.3) +
  geom_node_point(aes(colour=is_subgraph), size = 2) +
  coord_fixed() +
  theme_graph(base_family = "sans",base_size = 12)+
  theme(legend.position = "none",
        plot.margin = unit(c(0,0,0,0), 'lines')) +
  scale_colour_manual(values=c("lightgrey", "#92bbd9ff")) +
  ggtitle("graph showing factor 2 highlighted")

V(largest_component_graph)$is_subgraph = V(largest_component_graph)$factor == "Factor 3 (ARID1A rewired)"
ggraph(largest_component_graph, layout = l ) +
  geom_edge_link(alpha=0.3) +
  geom_node_point(aes(colour=is_subgraph), size = 2) +
  coord_fixed() +
  theme_graph(base_family = "sans",base_size = 12)+
  theme(legend.position = "none",
        plot.margin = unit(c(0,0,0,0), 'lines')) +
  scale_colour_manual(values=c("lightgrey", "#92bbd9ff")) +
  ggtitle("graph showing factor 3 highlighted")

dev.off()
