
#find correct environment 
packLib="/usr/lib/R"
if (file.exists(packLib)) {
  reticulate::use_condaenv("py37", required = T)
  setwd(dir = "~/Desktop/Melanoma_Resistance/")
}else {
  setwd(dir = "~/Desktop/Melanoma_Resistance/")
}

library(igraph)
library(EnsDb.Hsapiens.v86)
library("ggVennDiagram")
library(MOFA2)
library(stringr)
library(purrr)
library(ggrepel)
library(OmnipathR)
library(magrittr)
library(dplyr)
library(tidyr)

library(clusterProfiler)
library(enrichplot)
# we use ggplot2 to add x axis labels (ex: ridgeplot)
library(ggplot2)
library(msigdbr)
library(RCy3)

readr::local_edition(1) 

#get uninverse (background) for enrichment
universe<-read.delim(file = "~/Desktop/phuego_support/uniprot_to_gene.tab", sep = "\t", header = F)
universe<-unique(universe$V2)
#get mofa weights 
MOFAobject.trained<-load_model(file = "./results/mofa/mofa_object.hdf5")

weights <- get_weights(MOFAobject.trained, 
                       views = "all", 
                       as.data.frame = TRUE 
)
weights$node<-unlist(map(str_split(weights$feature, pattern = "_"),1))
weights$node<-unlist(map(str_split(weights$node, pattern = ";"),1))


#get phuego graphs
KDE <- "0.5"
factorS <- c("Factor1", "Factor2", "Factor3")
results_dir <- "./results/phuego/results/"
factor_graphs<-list()

#provisionally we need to find a way to make the phuego networks directed. 

for (factor in factorS) {
  file_graphml_up<-paste0(results_dir, factor, "/increased/KDE_", KDE, "/networks/KDE.graphml")
  file_graphml_down<-paste0(results_dir, factor, "/decreased/KDE_", KDE, "/networks/KDE.graphml")
  factor_graphs[[factor]][["up"]]<-read_graph(file=file_graphml_up,format = "graphml")
  factor_graphs[[factor]][["down"]]<-read_graph(file=file_graphml_down,format = "graphml")
  E(factor_graphs[[factor]][["up"]])$source <- "phuego"
  E(factor_graphs[[factor]][["down"]])$source <- "phuego"
}

get_nodes<-function(igraph_object){
  nodes<-V(igraph_object)$name
  genename_df <- AnnotationDbi::select(EnsDb.Hsapiens.v86, keys = nodes, 
                                       keytype = "UNIPROTID", 
                                       columns = "GENENAME")
  return(sort(genename_df$GENENAME))
}

#get nodes from all graphs
factor_nodes <- lapply(factor_graphs, function(sublist) {
  sapply(sublist, get_nodes)
})



enr_nodes<-function(nodes_vec, m_t2g, universe=NULL){
  enr_out<-enricher(nodes_vec, TERM2GENE = m_t2g)
  return(enr_out)
}

#perform enrichment

#see msig db collections
#msigdbr_collections()

# C2 (curated gene sets, 7233 gene sets)
# 
# CGP (chemical and genetic pertubations, 3438 gene sets)
# CP (canonical pathways, 3795 gene sets)
# CP:BIOCARTA (BioCarta gene sets, 292 gene sets)
# CP:PID (PID gene sets, 196 gene sets)
# CP:REACTOME (Reactome gene sets, 1692 gene sets)
# CP:WIKIPATHWAYS (WikiPathways gene sets, 791 gene sets)
# CP:KEGG (KEGG Legacy gene sets, 186 gene sets)

msig_cat<-"C2"
sub_cat<-"CP:REACTOME"
m_t2g <- msigdbr(species = "Homo sapiens", category = msig_cat, subcategory = sub_cat) %>% 
  dplyr::select(gs_name, human_gene_symbol)
colnames(m_t2g) <- c("term", "name")
#get nodes from all graphs
factor_enr <- lapply(factor_nodes, function(sublist) {
  sapply(sublist, enr_nodes, m_t2g, universe)
})

#graph construction

# Function to remove nodes from a graph while preserving paths
transform_graph <- function(graph) {
  # Extract nodes of type A and B
  nodes_A <- ends(graph, E(graph))[, 1]
  nodes_B <- ends(graph, E(graph))[, 2]

  # Identify edges between A-B-A
  # Function to get edges for each node in nodes_B
  get_edges_ABA <- function(node) {
    edges_ABA <- neighbors(graph, node, mode = "all")
    expand.grid(names(edges_ABA), names(edges_ABA))
  }

  # Apply the function to each node in nodes_B
  dfs_list <- map(nodes_B, get_edges_ABA)

  # Concatenate the resulting data frames
  result_df <- bind_rows(dfs_list)
  # Count occurrences of A-B-A edges
  condensed_data <- result_df %>%
  group_by_all() %>%
  summarise(weight = n())

  # Create a new weighted graph with nodes of type A
  new_graph <- graph_from_data_frame(condensed_data, directed = F)
  return(new_graph)
}

# Function to filter edges based on a quantile of edge weights
filter_edges_by_quantile <- function(graph, quantile_value) {
  # Check if the "weight" attribute is present in the graph
  if ("weight" %in% names(edge.attributes(graph))) {
    # Get the edge weights
    edge_weights <- E(graph)$weight
    
    # Determine the quantile threshold based on the provided quantile value
    threshold <- quantile(edge_weights, quantile_value)
    
    # Remove edges with weights under the quantile threshold
    graph_filtered <- delete_edges(graph, E(graph)[edge_weights < threshold])
    
    return(graph_filtered)
  } else {
    stop("The graph does not have a 'weight' attribute for edges.")
  }
}


# Function to process a direction ("up" or "down")
print_graph_to_cytoscape <- function(graph_data, direction, factor, sub_cat, msig_cat, n_pathways = 8, quantile_value=.25) {
  
  df.g <- graph.data.frame(d = graph_data[graph_data$direction == direction & graph_data$rank <= n_pathways, ], directed = FALSE)
  
  enr_subnet <- transform_graph(df.g)
  enr_subnet <- filter_edges_by_quantile(graph = enr_subnet, quantile_value = quantile_value)
  
  V(enr_subnet)$p.adj <- graph_data$p.adjust[match(V(enr_subnet)$name, graph_data$target)]
  V(enr_subnet)$genes <- graph_data$source[match(V(enr_subnet)$name, graph_data$target)]
  V(enr_subnet)$log10padj <- -log10(V(enr_subnet)$p.adj)
  
  RCy3::createNetworkFromIgraph(enr_subnet, title = direction, collection = paste0(factor, "_", sub_cat, "_", msig_cat))
  setVisualStyle("Curved")
  
  return(df.g)
}

# Loop over factors
enr_graphs<-list()
for (factor in factorS) {
  df <- factor_enr[[factor]]
  df<-factor_enr[[factor]]
  enr_network<-rbind(data.frame(target=df$up$ID,
                                source=df$up$geneID,
                                p.adjust=df$up$p.adjust,
                                GeneRatio=df$up$GeneRatio,
                                direction="up",
                                rank=c(1:length(df$up$ID))),
                     data.frame(target=df$down$ID,
                                source=df$down$geneID,
                                p.adjust=df$down$p.adjust,
                                GeneRatio=df$down$GeneRatio,
                                direction="down",
                                rank=c(1:length(df$down$ID))))
  enr_net_expanded<-tidyr::separate_rows(enr_network, source, sep="/")
  enr_net_expanded$target<-str_remove_all(enr_net_expanded$target, pattern = "REACTOME_")
  # Process "up" direction
  enr_graphs[[factor]]$up <- print_graph_to_cytoscape(enr_net_expanded, "up", factor, sub_cat, msig_cat)
  
  # Process "down" direction
  enr_graphs[[factor]]$down <- print_graph_to_cytoscape(enr_net_expanded, "down", factor, sub_cat, msig_cat)
}

