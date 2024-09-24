
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
seed<-runif(n = 1, min = 1, max = 1000)
set.seed(680.2898)
print(seed) #680.2898
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

factor_to_vis<-"Factor3"

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


# Define the union2 function
union2 <- function(g1, g2) {
  # Internal function to clean the names of a given attribute
  CleanNames <- function(g, target) {
    # Get target names
    gNames <- parse(text = paste0(target, "_attr_names(g)")) %>% eval
    # Find names that have a "_1" or "_2" at the end
    AttrNeedsCleaning <- grepl("(_\\d)$", gNames)
    # Remove the _x ending
    StemName <- gsub("(_\\d)$", "", gNames)

    NewnNames <- unique(StemName[AttrNeedsCleaning])
    # Replace attribute name for all attributes
    for (i in NewnNames) {
      attr1 <- parse(text = paste0(target, "_attr(g,'", paste0(i, "_1"), "')")) %>% eval
      attr2 <- parse(text = paste0(target, "_attr(g,'", paste0(i, "_2"), "')")) %>% eval

      g <- parse(text = paste0("set_", target, "_attr(g, i, value = ifelse(is.na(attr1), attr2, attr1))")) %>% eval
      g <- parse(text = paste0("delete_", target, "_attr(g,'", paste0(i, "_1"), "')")) %>% eval
      g <- parse(text = paste0("delete_", target, "_attr(g,'", paste0(i, "_2"), "')")) %>% eval
    }

    return(g)
  }

  g <- igraph::union(g1, g2)
  V(g)$consensus_direction <- paste0(V(g)$direction_1, "__", V(g)$direction_2)
  V(g)$consensus_factor <- paste0(V(g)$factor_1, "__", V(g)$factor_2, "__", V(g)$factor_3)

  # Loop through each attribute type in the graph and clean
  for (i in c("graph", "edge", "vertex")) {
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

  return(centrality_df)
}

# Create the centrality plots with titles
factor1_up_centrality<-centrality_plot(factor_graphs[["Factor1"]]$up, paste0("Factor1", " Up - Centrality Plot"), collapsed_factor_weights)
factor1_down_centrality<-centrality_plot(factor_graphs[["Factor1"]]$down, paste0("Factor1", " Down - Centrality Plot"), collapsed_factor_weights)

factor3_up_centrality<-centrality_plot(factor_graphs[["Factor3"]]$up, paste0("Factor3", " Up - Centrality Plot"), collapsed_factor_weights)
factor3_down_centrality<-centrality_plot(factor_graphs[["Factor3"]]$down, paste0("Factor3", " Down - Centrality Plot"), collapsed_factor_weights)


factor1_centrality<-rbind(factor1_up_centrality, factor1_down_centrality)
factor3_centrality<-rbind(factor3_up_centrality, factor3_down_centrality)


# Order and rank centrality
factor1_centrality <- factor1_centrality[order(factor1_centrality$values, decreasing = TRUE), ]
factor1_centrality$rank <- 1:nrow(factor1_centrality)

factor3_centrality <- factor3_centrality[order(factor3_centrality$values, decreasing = TRUE), ]
factor3_centrality$rank <- 1:nrow(factor3_centrality)



# Retrieve gene names
genename_df <- AnnotationDbi::select(EnsDb.Hsapiens.v86, keys = c(as.character(factor1_centrality$ind), as.character(factor3_centrality$ind)),
                                     keytype = "UNIPROTID",
                                     columns = "GENENAME")
factor1_centrality$names <- genename_df$GENENAME[match(factor1_centrality$ind, genename_df$UNIPROTID)]
factor3_centrality$names <- genename_df$GENENAME[match(factor3_centrality$ind, genename_df$UNIPROTID)]


factor1_centrality$names[factor1_centrality$rank > 30] <- ""
factor3_centrality$names[factor3_centrality$rank > 30] <- ""

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


g<-union2(union_factor_graphs[["Factor1"]],
          union_factor_graphs[["Factor3"]])

conv_nodes<-data.frame(uniprt=V(g)$name,
                       gene_name=V(g)$Gene_name)

#######PCSF#######
###PCSF on the most central nodes (pagerank)

library(OmnipathR)
library(PCSF)
library(readr)

cut_off<-50
names_for_subnet <- unique(c(factor1_centrality$ind[factor1_centrality$rank < cut_off],
                             factor3_centrality$ind[factor3_centrality$rank < cut_off]))

####prep data for pcsf#####
terminal<-rep(1, length(names_for_subnet))
names(terminal)<-names_for_subnet
####parameters ####
n<-4000 #no. of runs30
r<-5 #adding random noise to edge costs  5
w<-40 #number of trees in output 40
b<-8 #tuning node prizes 1
mu<-0.005 #hub penalisation #0.005

######perfor pcsf n times adding some noise to produce more robust network. #####
subnet <- PCSF_rand(g,
                    terminal,
                    n = n,
                    r = r,
                    w = w,
                    b = b,
                    mu = mu)

subnet_l <- igraph::layout_with_graphopt(subnet)
V(subnet)$Gene_name <- conv_nodes$gene_name[match(V(subnet)$name, conv_nodes$uniprt)]
L_df<-as.data.frame(subnet_l)
colnames(L_df)<-c("x", "y")
L_df$Gene_name <- V(subnet)$Gene_name
L_df$Factor1 <- L_df$Gene_name %in% factor_genes_names_df$values[factor_genes_names_df$factor == "Factor1"]
L_df$Factor2 <- L_df$Gene_name %in% factor_genes_names_df$values[factor_genes_names_df$factor == "Factor2"]
L_df$Factor3 <- L_df$Gene_name %in% factor_genes_names_df$values[factor_genes_names_df$factor == "Factor3"]

# Remove the 'feature' column and filter out 'phospho' views
factor_weights_wide <- weights %>%
  dplyr::select(-feature) %>%
  dplyr::filter(view != "phospho") %>%
  pivot_wider(names_from = view, values_from = value)

# Function to create columns for each factor and view
add_factor_view_columns <- function(df, weights_in, factor) {
  factor_weights_subset <- weights_in %>%
    dplyr::filter(factor == !!factor) %>%
    dplyr::select(node, mRNA, protein)
  df <- df %>%
    left_join(factor_weights_subset, by = c("Gene_name" = "node")) %>%
    rename_with(~ paste0(., "_", factor), c(mRNA, protein))
  return(df)
}

# Add columns for each factor and view
L_df <- L_df %>%
  add_factor_view_columns(factor_weights_wide, "Factor1") %>%
  add_factor_view_columns(factor_weights_wide, "Factor2") %>%
  add_factor_view_columns(factor_weights_wide, "Factor3")

save.image(file='./results/heatdiffusion/data_for_heat_diffusion.Rdata')

library(diffusr)

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
    width = 8, height = 8)

# Usage example
L_df_processed <- process_data_for_labeling(L_df, "arid1a_up_probabilities", "combined_up_probabilities", "Gene_name")
L_df_processed<-L_df_processed[arid1a_up_receptors == 0 | combined_up_receptors == 0,]

ggplot(L_df_processed, aes(x=arid1a_up_probabilities, y=combined_up_probabilities, label=label_flag, color = is_receptor)) +
  geom_point() +
  ggtitle("Comparison of 'heat' of nodes in network : up receptors vs up") +
  ggrepel::geom_text_repel(size = 6, colour = "black") +
  cowplot::theme_cowplot() +
  theme(
    plot.title = element_text(size = 15, face = "bold"),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 1),
  ) +
  grids(linetype = "dashed")+ scale_color_manual(values=wes_palette("Chevalier1"))

# basic scatterplot
L_df_processed <- process_data_for_labeling(L_df, "arid1a_down_probabilities", "combined_down_probabilities", "Gene_name")
L_df_processed<-L_df_processed[arid1a_down_receptors == 0 | combined_down_receptors == 0,]

ggplot(L_df_processed, aes(x=arid1a_down_probabilities, y=combined_down_probabilities, label=label_flag, color = is_receptor)) +
  geom_point() +
  ggtitle("Comparison of 'heat' of nodes in network : down receptors vs down") +
  ggrepel::geom_text_repel(size = 8, colour = "black") +
  cowplot::theme_cowplot() +
  theme(
    plot.title = element_text(size = 15, face = "bold"),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 1),
  ) +
  grids(linetype = "dashed")+ scale_color_manual(values=wes_palette("Chevalier1"))

# basic scatterplot
L_df_processed <- process_data_for_labeling(L_df, "arid1a_up_probabilities", "combined_down_probabilities", "Gene_name")
L_df_processed<-L_df_processed[arid1a_up_receptors == 0 | combined_down_receptors == 0,]

ggplot(L_df_processed, aes(x=arid1a_up_probabilities, y=combined_down_probabilities, label=label_flag, color = is_receptor)) +
  geom_point() +
  ggtitle("Comparison of 'heat' of nodes in network : up receptors vs down") +
  ggrepel::geom_text_repel(size = 8, colour = "black") +
  cowplot::theme_cowplot() +
  theme(
    plot.title = element_text(size = 15, face = "bold"),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 1),
  ) +
  grids(linetype = "dashed")+ scale_color_manual(values=wes_palette("Chevalier1"))

# basic scatterplot
L_df_processed <- process_data_for_labeling(L_df, "arid1a_down_probabilities", "combined_up_probabilities", "Gene_name")
L_df_processed<-L_df_processed[arid1a_down_receptors == 0 | combined_up_receptors == 0,]

ggplot(L_df_processed, aes(x=arid1a_down_probabilities, y=combined_up_probabilities, label=label_flag, color = is_receptor)) +
  geom_point() +
  ggtitle("Comparison of 'heat' of nodes in network : down receptors vs up") +
  ggrepel::geom_text_repel(size = 8, colour = "black") +
  cowplot::theme_cowplot() +
  theme(
    plot.title = element_text(size = 15, face = "bold"),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 1),
  ) +
  grids(linetype = "dashed")+ scale_color_manual(values=wes_palette("Chevalier1"))
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

  # Ensure super_graph is an igraph object
  if (!is_igraph(super_graph)) {
    stop("super_graph is not a valid igraph object.")
  }

  # Debugging output
  cat("Type of super_graph:", class(super_graph), "\n")
  cat("Number of vertices in super_graph:", length(V(super_graph)), "\n")
  cat("Number of edges in super_graph:", length(E(super_graph)), "\n")

  # Display the list of subgraphs
  test <- super_graph

  subnet_l <- tryCatch({
    igraph::layout_with_sugiyama(test)$layout
  }, error = function(e) {
    stop("Layout calculation failed: ", e$message)
  })

  if (is.null(subnet_l)) {
    stop("Layout calculation failed.")
  }
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
    geom_node_point(size = 5) +
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

pdf(file = paste0("./results/heatdiffusion/shortestpath_results.pdf"),
    width = 30, height = 40)

# Plot shortest paths from specific receptors to transcription factors JUN and PRKD1
make_shortest_paths_plot(
  sources = c("EGFR", "ROS1", "FGFR1"),
  sinks =  c("JUN"),
  g = subnet,
  conv_nodes = conv_nodes,
  factor_genes_names_df = factor_genes_names_df,
  factor_weights_wide = factor_weights_wide,
  pal = pal
) +
  ggtitle("Shortest Paths from EGFR, ROS1, and FGFR1 to JUN and PRKD1")

# Plot shortest paths from specific receptors to MAPK1
make_shortest_paths_plot(
  sources = c("EGFR", "ROS1", "FGFR1"),
  sinks =  c("MAPK1"),
  g = subnet,
  conv_nodes = conv_nodes,
  factor_genes_names_df = factor_genes_names_df,
  factor_weights_wide = factor_weights_wide,
  pal = pal
) +
  ggtitle("Shortest Paths from EGFR, ROS1, and FGFR1 to MAPK1")

# Plot shortest paths from ITGA4 and NTRK3 to MAPK1
make_shortest_paths_plot(
  sources = c("ITGA4", "NTRK3"),
  sinks =  c("MAPK1", "MAPK3"),
  g = subnet,
  conv_nodes = conv_nodes,
  factor_genes_names_df = factor_genes_names_df,
  factor_weights_wide = factor_weights_wide,
  pal = pal
) +
  ggtitle("Shortest Paths from ITGA4, NTRK3, SIRPA and ERBB3 to MAPK1")

# Plot shortest paths from NGFR to SPRED1
make_shortest_paths_plot(
  sources = c("NGFR"),
  sinks =  c("SPRED1"),
  g = subnet,
  conv_nodes = conv_nodes,
  factor_genes_names_df = factor_genes_names_df,
  factor_weights_wide = factor_weights_wide,
  pal = pal
) +
  ggtitle("Shortest Path from NGFR to SPRED1")

# Plot shortest paths from NGFR to MAPK1
make_shortest_paths_plot(
  sources = c("NGFR"),
  sinks =  c("MAPK1"),
  g = subnet,
  conv_nodes = conv_nodes,
  factor_genes_names_df = factor_genes_names_df,
  factor_weights_wide = factor_weights_wide,
  pal = pal
) +
  ggtitle("Shortest Path from NGFR to MAPK1")

dev.off()
