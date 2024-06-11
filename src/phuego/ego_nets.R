

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

source("./src/functions/default_variables.R")

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
  out <- ggplot(centrality_df, aes(x = rank, y = values, label = names, colour=factor_weight)) + 
    geom_point() + 
    cowplot::theme_cowplot() + 
    geom_text_repel(colour="black") +
    scale_colour_gradientn(colours = pal) + 
    ggtitle(title) # Add the title here 
  
  return(out)
}

# Create the centrality plots with titles
centrality_plots <- plot_grid(
  centrality_plot(factor_graphs[[factor_to_vis]]$up, "Factor 1 Up - Centrality Plot", collapsed_factor_weights), 
  centrality_plot(factor_graphs[[factor_to_vis]]$down, "Factor 1 Down - Centrality Plot", collapsed_factor_weights),
  ncol = 2, align = "h", rel_widths = c(3, 3)
)



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

conv_nodes<-data.frame(uniprt=V(g)$name,
                       gene_name=V(g)$Gene_name)


# Define the function to create an ego network plot
create_ego_net_plot <- function(gene_name, g, conv_nodes, factor_genes_names_df, weights) {
  # Find the node
  node <- conv_nodes$uniprt[conv_nodes$gene_name == gene_name]
  
  # Create the ego network graph
  ego_net <- make_ego_graph(g, order = 1, nodes = node)[[1]]
  
  # Generate the plot
  l <- igraph::layout.circle(ego_net)
  L_df<-as.data.frame(l)
  colnames(L_df)<-c("x", "y")
  L_df$Gene_name <- V(ego_net)$Gene_name
  L_df$Factor1 <- L_df$Gene_name %in% factor_genes_names_df$values[factor_genes_names_df$factor == "Factor1"]
  L_df$Factor2 <- L_df$Gene_name %in% factor_genes_names_df$values[factor_genes_names_df$factor == "Factor2"]
  L_df$Factor3 <- L_df$Gene_name %in% factor_genes_names_df$values[factor_genes_names_df$factor == "Factor3"]
  
  # Calculate centroid for all factors combined
  centroids <- L_df %>%
    filter(Factor1 == TRUE | Factor2 == TRUE | Factor3 == TRUE) %>%
    summarize(
      x = mean(x),
      y = mean(y),
      Gene_name = "Centroid_All_Factors",
      Factor1 = TRUE,
      Factor2 = TRUE,
      Factor3 = TRUE
    )
  # Bind the centroid rows to the original data frame
  L_df <- bind_rows(L_df, unique(centroids))
  
  # Remove the 'feature' column and filter out 'phospho' views
  factor_weights_wide <- weights %>%
    select(-feature) %>%
    filter(view != "phospho") %>%
    pivot_wider(names_from = view, values_from = value)

  # Function to create columns for each factor and view
  add_factor_view_columns <- function(df, weights_in, factor) {
    factor_weights_subset <- weights_in %>%
      filter(factor == !!factor) %>%
      select(node, mRNA, protein)
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

  # Generate the plot
  ego_plot <- ggraph(ego_net, layout = l, circular = TRUE) + 
    geom_convexhull(data = L_df[L_df$Factor1==1,], aes(x = x, y = y, fill = "Factor 1"), alpha = 0.3) + 
    geom_convexhull(data = L_df[L_df$Factor2==1,], aes(x = x, y = y, fill = "Factor 2"), alpha = 0.3) + 
    geom_convexhull(data = L_df[L_df$Factor3==1,], aes(x = x, y = y, fill = "Factor 3"), alpha = 0.3) + 
    geom_point(data = L_df[L_df$Gene_name != "Centroid_All_Factors",], 
               aes(x = x, y = y, colour = protein_Factor1), 
               alpha = 1, size = 30,
               stroke = 5) +
    geom_point(data = L_df[L_df$Gene_name != "Centroid_All_Factors",], 
               aes(x = x, y = y, colour = mRNA_Factor1), 
               alpha = 1, size = 20,
               stroke = 5) + 
    geom_edge_arc(start_cap = circle(3, 'mm'),
                  end_cap = circle(3, 'mm'), 
                  aes(alpha = weight^2)) + 
    geom_node_point(size = 5) + 
    coord_fixed() + theme_void() +
    geom_node_label(aes(label = Gene_name), repel = FALSE) +
    scale_fill_manual(
      name = "Factors",
      values = c("Factor 1" = "#92bbd9ff", "Factor 2" = "#dcca2cff", "Factor 3" = "#6fb382ff")
    ) +
    labs(title = paste0(gene_name, " ego net")) +
    scale_colour_gradientn(colours = pal) 
  
  
  return(ego_plot)
}

conv_nodes<-data.frame(uniprt=V(g)$name,
                       gene_name=V(g)$Gene_name)

list_of_inputs<-list(
  phospho=data.frame(read.csv("./data/input_data/phosphosites.csv"))
  ,protein=data.frame(read.csv("./data/input_data/proteins.csv"))
  ,mRNA=data.frame(read.csv("./data/input_data/rna_expression.csv"))
)

# Define the drugs and their corresponding colors
drug_colors <- c(
  "Untreated" = "#A2AEBB",
  "Vemurafenib" = "#FFBA08",
  "Trametinib" = "#D00000",
  "Combination" = "#3F88C5"
)


plot_raw_data_ego<-function(node_gene_name, #node which the ego is based on
                            feature_set #mRNA or protein
                            )
{
  node= conv_nodes$uniprt[conv_nodes$gene_name == node_gene_name]
  ego_net <- make_ego_graph(g, order = 1, nodes = node)[[1]]
  
  subset_proteins<-V(ego_net)$Gene_name
  
  if (!feature_set=="protein") {
    rtk_mrna<-reshape2::melt(list_of_inputs[[feature_set]][list_of_inputs[[feature_set]]$X %in% subset_proteins,])
    rtk_mrna$name <- factor(rtk_mrna$X, levels= subset_proteins)
    y_axis_label = "RNA expression"
    main_title = paste(y_axis_label, " of ", node_gene_name, " interacting partners")
    
  }else{
    uniprot_codes <- conv_nodes$uniprt[match(subset_proteins, conv_nodes$gene_name)]
    rtk_mrna<-reshape2::melt(list_of_inputs[[feature_set]][list_of_inputs[[feature_set]]$X %in% uniprot_codes,])
    rtk_mrna$name <- conv_nodes$gene_name[match(rtk_mrna$X, conv_nodes$uniprt)]
    y_axis_label = "Protein abundance"
    main_title = paste(y_axis_label, " of ", node_gene_name, " interacting partners")
    
  }
  rtk_mrna$variable<-sub("*\\.[0-9]", "", rtk_mrna$variable)
  rtk_mrna <- rtk_mrna |>
    separate_wider_delim(variable, delim = "__", names = c("drug", "ko"))
  
  rtk_mrna$drug<-names(replacement_Vec)[match(rtk_mrna$drug, unname(replacement_Vec))]
  # Define the desired order of levels
  desired_order <- c("Untreated", "Vemurafenib", "Trametinib", "Combinations")
  # Convert my_column to a factor with the specified order
  rtk_mrna$drug <- factor(rtk_mrna$drug, levels = desired_order)
  plot_out<-rtk_mrna %>%
    ggplot( aes(x=drug, y=value, fill=ko)) +
    geom_boxplot() +
    geom_jitter(color="black", size=0.4, alpha=0.9) +
    scale_fill_manual(values = gene_colours) +
    cowplot::theme_cowplot() +
    theme(
      plot.title = element_text(size=18),
      axis.text.x = element_text(angle = 70, hjust = 1, size = rel(1)),
    ) +
    xlab("") +
    facet_wrap(~name) + 
    grids(linetype = "dashed")+
    labs(
      x = "Drug treatment",
      y = y_axis_label,
      title = main_title
    )
  return(plot_out)
}


pdf(file = paste0("~/Desktop/Melanoma_Resistance/paper/for_sumana/", 
                  factor_to_vis, 
                  "/centrality.pdf"), 
    width = 15, height = 15)

centrality_plots

# Example gene names
gene_names <- c("MAPK1", "EGFR", "IGF1R", "PRKD1")

# Generate plots for each gene name
plots <- lapply(gene_names, create_ego_net_plot, g = g, conv_nodes = conv_nodes, factor_genes_names_df = factor_genes_names_df, weights=weights)

# Display the plots
for (plot in plots) {
  print(plot)
}

plot_raw_data_ego(node_gene_name = "MAPK1", "protein")
plot_raw_data_ego(node_gene_name = "MAPK1", "mRNA")

plot_raw_data_ego(node_gene_name = "EGFR", "protein")
plot_raw_data_ego(node_gene_name = "EGFR", "mRNA")


plot_raw_data_ego(node_gene_name = "IGF1R", "protein")
plot_raw_data_ego(node_gene_name = "IGF1R", "mRNA")


plot_raw_data_ego(node_gene_name = "PRKD1", "protein")
plot_raw_data_ego(node_gene_name = "PRKD1", "mRNA")

dev.off()
