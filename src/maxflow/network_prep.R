
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


readr::local_edition(1) 

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
for (factor in factorS) {
  file_graphml_up<-paste0(results_dir, factor, "/increased/KDE_", KDE, "/networks/KDE.graphml")
  file_graphml_down<-paste0(results_dir, factor, "/decreased/KDE_", KDE, "/networks/KDE.graphml")
  factor_graphs[[factor]][["up"]]<-read_graph(file=file_graphml_up,format = "graphml")
  factor_graphs[[factor]][["down"]]<-read_graph(file=file_graphml_down,format = "graphml")
}

#get tf activities 
arid1a_tf_acts<-read.csv(file = "./results/transcriptomics/tf_activity/arid1a_tf_acts.csv")
arid1a_tf_pval<-read.csv(file = "./results/transcriptomics/tf_activity/arid1a_tf_pval.csv")
arid1a_tf<-merge(arid1a_tf_acts,arid1a_tf_pval,by="X")
colnames(arid1a_tf)<-c("TF", "Activity", "P_val")
arid1a_tf$diffexpressed <- arid1a_tf$Activity < 0
arid1a_tf$diffexpressed[arid1a_tf$P_val>0.001] <- "Not significant"
arid1a_tf$diffexpressed [arid1a_tf$diffexpressed == T] <- "Downregulated"
arid1a_tf$diffexpressed [arid1a_tf$diffexpressed == F] <- "Upregulated"
arid1a_tf$label<-arid1a_tf$TF
arid1a_tf$TF[arid1a_tf$diffexpressed=="Not significant"]<- ""
sig_tfs<-arid1a_tf[arid1a_tf$diffexpressed!="Not significant",]$TF
ppi <- read.csv(file = "./data/omnipath.csv")
#get incoming edges to sig tfs
tf_uniprt <- AnnotationDbi::select(EnsDb.Hsapiens.v86, keys = sig_tfs, keytype = "GENENAME", columns = "UNIPROTID")
ppi_tf<-ppi[ppi$target %in% unique(tf_uniprt$UNIPROTID),]

ppi_tf_expanded <- ppi_tf %>%
  mutate(
    source = ifelse(grepl("^COMPLEX:", source), strsplit(source, "_")[[1]], source),
    target = ifelse(grepl("^COMPLEX:", target), strsplit(target, "_")[[1]], target)
  ) %>%
  unique()

# Convert the 'source' and 'target' columns to a directed graph
graph_tf <- graph_from_data_frame(data.frame(source=ppi_tf_expanded$source,
                                          target=ppi_tf_expanded$target),
                               directed = F)

# If you want to add additional attributes to the graph vertices or edges, you can do so like this:
# For example, adding the 'is_directed', 'is_stimulation', 'is_inhibition', and 'consensus_direction' attributes to edges
E(graph_tf)$is_directed <- ppi_tf_expanded$is_directed
E(graph_tf)$is_stimulation <- ppi_tf_expanded$is_stimulation
E(graph_tf)$is_inhibition <- ppi_tf_expanded$is_inhibition
E(graph_tf)$consensus_direction <- ppi_tf_expanded$consensus_direction

tf_name<-AnnotationDbi::select(EnsDb.Hsapiens.v86, keys = unique(ppi_tf_expanded$target), keytype = "UNIPROTID", columns = "GENENAME")
expan_tf_name<-AnnotationDbi::select(EnsDb.Hsapiens.v86, keys = unique(ppi_tf_expanded$source), keytype = "UNIPROTID", columns = "GENENAME")
expan_tf<-unique(c(expan_tf_name$GENENAME, tf_name$GENENAME))

# ggplot(data = arid1a_tf, aes(x = Activity, y = -log10(P_val), col = diffexpressed, label = TF)) +
#   geom_vline(xintercept = c(-0.6, 0.6), col = "gray", linetype = 'dashed') +
#   geom_hline(yintercept = -log10(0.05), col = "gray", linetype = 'dashed') +
#   geom_point(size = 2) +
#   scale_color_manual(values = c("#00AFBB", "grey", "#bb0c00"), # to set the colours of our variable
#                      labels = c("Downregulated", "Not significant", "Upregulated")) + # to set the labels in case we want to overwrite the categories from the dataframe (UP, DOWN, NO)
#   labs(color = 'Severe', #legend_title,
#        x = expression("log"[2]*"FC"), y = expression("-log"[10]*"p-value")) +
#   scale_x_continuous(breaks = seq(-10, 10, 2)) + # to customise the breaks in the x axis
#   ggtitle('TF-activity in ARID1A KO vs WT cells') + # Plot title
#   geom_text_repel(max.overlaps = Inf)+ # To show all labels
#   cowplot::theme_cowplot()

get_nodes<-function(igraph_object){
  nodes<-V(igraph_object)$name
  genename_df <- AnnotationDbi::select(EnsDb.Hsapiens.v86, keys = nodes, 
                                       keytype = "UNIPROTID", 
                                       columns = "GENENAME")
  return(genename_df$GENENAME)
}

join_tfs <- function(graph, tf_graph) {
  reg_graph_union <- igraph::union(graph, tf_graph)
  largest_contiguous_subgraph <- largest_component(reg_graph_union)
  return(largest_contiguous_subgraph)
}

#merge each factor graph with the tf graph 
factor_graphs <- lapply(factor_graphs, function(graph) {
  sapply(graph, join_tfs, graph_tf)
})
#get nodes from all graphs
factor_nodes <- lapply(factor_graphs, function(sublist) {
  sapply(sublist, get_nodes)
})

#get receptors or ligands as the starting points in the maximum flow. 

ligands <-
  OmnipathR::import_omnipath_intercell(
    parent = 'ligand',
    topology = 'sec',
    consensus_percentile = 50,
    loc_consensus_percentile = 30,
    entity_types = 'protein'
  ) %>%
  pull(genesymbol) %>%
  unique %T>%
  length
receptors <-
  OmnipathR::import_omnipath_intercell(
    parent = 'receptor',
    topology = 'pmtm',
    consensus_percentile = 50,
    loc_consensus_percentile = 30,
    entity_types = 'protein'
  ) %>%
  pull(genesymbol) %>%
  unique %T>%
  length

vennList<-list()
for (factor in factorS) {
  nodes<-factor_nodes[[factor]]
  nodes$starting_points<-unique(c(receptors, ligands))
  nodes$end_points<-expan_tf
  # Default plot
  vennPlot<-ggVennDiagram(nodes)
  vennList[[factor]]<-vennPlot
}

# Convert nested list to data frame with the factor weight from mofa
add_mofa_weight <- function(nodes,weights, factor) {
  df_out<-data.frame(
    Node = nodes,
    stringsAsFactors = FALSE,
    Factor=factor,
    receptor_ligands = nodes %in% c(receptors, ligands),
    tf= nodes %in% unique(tf_uniprt$GENENAME)
  )
  df_out<-merge(x=df_out,y=weights,
            by.x=c("Node", "Factor"), by.y=c("node", "factor"))
  return(df_out)
}

#add the sink (transcription factors) and source (receptor and ligands) nodes 

#get nodes from all graphs
sink_source_egdeList<-list()
for (factor in factorS) {
  out<-list()
  out$up<-add_mofa_weight(factor_nodes[[factor]]$up, weights, factor)
  out$down<-add_mofa_weight(factor_nodes[[factor]]$down, weights, factor)
  sink_source_egdeList[[factor]]<-out
}



