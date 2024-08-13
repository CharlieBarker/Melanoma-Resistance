

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
  geom_point() + 
  cowplot::theme_cowplot() + 
  geom_text_repel(colour = "black", force = 15) +
  scale_colour_gradientn(colours = pal) + 
  labs(x = "Rank in network", y = "Centrality (PageRank)", 
       title = "Centrality of Network describing combination-specific changes") +  # Adding axis titles
  theme(legend.position = "bottom", plot.title = element_text(size=10))  # Placing legend at the bottom





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


#get mofa weights 
MOFAobject.trained<-load_model(file = "./results/mofa/mofa_object.hdf5")

weights <- get_weights(MOFAobject.trained, 
                       views = "all", 
                       as.data.frame = TRUE 
)


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

####PCSF on the most central nodes (pagerank)####

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
n<-30 #no. of runs30
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


library(diffusr)

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
  pt <- random.walk(start_vector, data.matrix(adj_matrix), r = 0.02, correct.for.hubs = T)
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


pdf(file = paste0("/Users/charliebarker/Desktop/test.pdf"), 
    width = 18, height = 15)

ggraph(subnet, layout = subnet_l) + 
  geom_point(data = L_df[L_df$Gene_name != "Centroid_All_Factors",], 
             aes(x = x, y = y, colour = protein_Factor3), 
             alpha = 1, size = 10,
             stroke = 5) +
  geom_point(data = L_df[L_df$Gene_name != "Centroid_All_Factors",], 
             aes(x = x, y = y, colour = mRNA_Factor3), 
             alpha = 1, size = 6,
             stroke = 5) + 
  geom_edge_link(start_cap = circle(3, 'mm'),
                 end_cap = circle(3, 'mm'), 
                 aes(alpha = weight)) + 
  geom_node_point(size = 5) + 
  coord_fixed() + theme_void() +
  geom_node_label(aes(label = Gene_name),size=4, repel = FALSE) +
  theme(legend.position = "bottom") + # Placing legend at the bottom
  scale_colour_gradientn(colours = pal, na.value = "lightgrey")

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

# basic scatterplot
L_df$is_receptor = L_df$Gene_name %in% receptors
ggplot(L_df, aes(x=arid1a_up_probabilities, y=combined_up_probabilities, label=Gene_name, color = is_receptor)) + 
  geom_point() +
  ggtitle("Comparison of 'heat' of nodes in network : up receptors vs up") +
  ggrepel::geom_text_repel(size = 8) + cowplot::theme_cowplot()

# basic scatterplot
ggplot(L_df, aes(x=arid1a_down_probabilities, y=combined_down_probabilities, label=Gene_name, color = is_receptor)) + 
  geom_point() +
  ggtitle("Comparison of 'heat' of nodes in network : down receptors vs down") +
  ggrepel::geom_text_repel(size = 8) + cowplot::theme_cowplot()

# basic scatterplot
ggplot(L_df, aes(x=arid1a_up_probabilities, y=combined_down_probabilities, label=Gene_name, color = is_receptor)) + 
  geom_point() +
  ggtitle("Comparison of 'heat' of nodes in network : up receptors vs down") +
  ggrepel::geom_text_repel(size = 8) + cowplot::theme_cowplot()

# basic scatterplot
ggplot(L_df, aes(x=arid1a_down_probabilities, y=combined_up_probabilities, label=Gene_name, color = is_receptor)) + 
  geom_point() +
  ggtitle("Comparison of 'heat' of nodes in network : down receptors vs up") +
  ggrepel::geom_text_repel(size = 8) + cowplot::theme_cowplot()

dev.off()


# Function to convert vector to binary matrix
vector_to_binary_matrix <- function(vec) {
  # Find indices of non-zero elements
  nonzero_indices <- which(vec != 0)
  
  # Number of non-zero elements
  num_nonzero <- length(nonzero_indices)
  
  # Create an empty matrix
  result_matrix <- matrix(0, nrow = length(vec), ncol = num_nonzero)
  
  # Populate the matrix
  for (i in 1:num_nonzero) {
    result_matrix[nonzero_indices[i], i] <- 1
  }
  
  # Return the resulting matrix
  return(result_matrix)
}

to_test<-combined_down_receptors
test<-data.frame(L_df$Gene_name, perform_random_walk(subnet, vector_to_binary_matrix(to_test)))

test[test$L_df.Gene_name %in% c("JUN", "PRKD1"),]

L_df$Gene_name[as.logical(rowSums(vector_to_binary_matrix(to_test)))]