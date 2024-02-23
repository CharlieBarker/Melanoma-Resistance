
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
library(ggpubr)
library(wesanderson)
#globally improtant variables

Abs<-F #are we looking at up and down regulation seperately?



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

#provisionally we need to find a way to make the phuego networks directed. 

for (factor in factorS) {
  file_graphml_up<-paste0(results_dir, factor, "/increased/KDE_", KDE, "/networks/KDE.graphml")
  file_graphml_down<-paste0(results_dir, factor, "/decreased/KDE_", KDE, "/networks/KDE.graphml")
  factor_graphs[[factor]][["up"]]<-read_graph(file=file_graphml_up,format = "graphml")
  factor_graphs[[factor]][["down"]]<-read_graph(file=file_graphml_down,format = "graphml")
  E(factor_graphs[[factor]][["up"]])$source <- "phuego"
  E(factor_graphs[[factor]][["down"]])$source <- "phuego"
  E(factor_graphs[[factor]][["up"]])$capacity <- E(factor_graphs[[factor]][["up"]])$weight
  E(factor_graphs[[factor]][["down"]])$capacity <- E(factor_graphs[[factor]][["down"]])$weight
}

#change Factor1 weights to be the opposite because up is down and down is up. 
factor1_down<-factor_graphs[["Factor1"]][["up"]]
factor1_up<-factor_graphs[["Factor1"]][["down"]]

factor_graphs[["Factor1"]][["up"]]<-factor1_up
factor_graphs[["Factor1"]][["down"]]<-factor1_down

#get tf activities 
arid1a_tf_acts<-read.csv(file = "./results/transcriptomics/tf_activity/arid1a_tf_acts.csv")
arid1a_tf_pval<-read.csv(file = "./results/transcriptomics/tf_activity/arid1a_tf_pval.csv")
arid1a_tf<-merge(arid1a_tf_acts,arid1a_tf_pval,by="X")
colnames(arid1a_tf)<-c("TF", "Activity", "P_val")
arid1a_tf$diffexpressed <- arid1a_tf$Activity < 0
arid1a_tf$diffexpressed[arid1a_tf$P_val>0.1] <- "Not significant"
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
                               directed = T)

# If you want to add additional attributes to the graph vertices or edges, you can do so like this:
# For example, adding the 'is_directed', 'is_stimulation', 'is_inhibition', and 'consensus_direction' attributes to edges
E(graph_tf)$is_directed <- ppi_tf_expanded$is_directed
E(graph_tf)$is_stimulation <- ppi_tf_expanded$is_stimulation
E(graph_tf)$is_inhibition <- ppi_tf_expanded$is_inhibition
E(graph_tf)$n_ref <- ppi_tf_expanded$n_references

E(graph_tf)$source <- "omnipath_tf_network"

get_nodes<-function(igraph_object){
  nodes<-V(igraph_object)$name
  genename_df <- AnnotationDbi::select(EnsDb.Hsapiens.v86, keys = nodes, 
                                       keytype = "UNIPROTID", 
                                       columns = "GENENAME")
  return(genename_df$GENENAME)
}

join_tfs <- function(graph, tf_graph, list_of_tfs) {
  #find the intersection of the tf network and the pheugo network 
  reg_graph_intersect <- igraph::intersection(as.directed(graph, mode = "mutual"), tf_graph)
  list_to_keep<-unique(V(reg_graph_intersect)$name, list_of_tfs)
  #only keep the parts of the tf graph that are connecting to the phuego graph 
  #find the union of these two graphs
  reg_graph_union <- igraph::union(as.directed(graph, mode = "mutual"), tf_graph)
  #remove these graphs
  #largest_contiguous_subgraph <- largest_component(reg_graph_union)
  #E(largest_contiguous_subgraph)$capacity<-E(largest_contiguous_subgraph)$weight
  return(reg_graph_union)
}
#get nodes from all graphs
factor_nodes <- lapply(factor_graphs, function(sublist) {
  sapply(sublist, get_nodes)
})
#merge each factor graph with the tf graph 
factor_graphs <- lapply(factor_graphs, function(graph) {
  sapply(graph, join_tfs, graph_tf, sig_tfs)
})
#get nodes from all graphs
factor_nodes <- lapply(factor_graphs, function(sublist) {
  sapply(sublist, get_nodes)
})
#turn all graphs into genenames not uniprot ids
turn_to_genenames<-function(igraph_object){
  # Extract the names of the vertices (nodes)
  nodes <- V(igraph_object)$name
  
  # Query the annotation database to retrieve gene names corresponding to UNIPROT IDs
  genename_df <- AnnotationDbi::select(EnsDb.Hsapiens.v86, keys = nodes, 
                                       keytype = "UNIPROTID", 
                                       columns = "GENENAME")
  
  # Create a mapping between UNIPROT IDs and gene names
  name_mapping <- setNames(genename_df$GENENAME, genename_df$UNIPROTID)
  
  # Manually correct gene names for specific UNIPROT IDs
  name_mapping["P42771"] <- "CDKN2A"
  name_mapping["Q8N726"] <- "ARF"
  
  # Update vertex names using the name_mapping
  V(igraph_object)$name <- ifelse(V(igraph_object)$name %in% names(name_mapping),
                                  name_mapping[V(igraph_object)$name],
                                  V(igraph_object)$name)
  return(igraph_object)
}
factor_graphs<-lapply(factor_graphs, function(graph) {
  sapply(graph, turn_to_genenames)
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

# Convert nested list to data frame with the factor weight from mofa
add_mofa_weight <- function(nodes,
                            weights, 
                            factor, 
                            Abs=T) {
  df_out<-data.frame(
    Node = nodes,
    stringsAsFactors = FALSE,
    Factor=factor,
    receptor_ligands = nodes %in% c(receptors, ligands),
    tf= nodes %in% unique(tf_uniprt$GENENAME)
  )
  df_out<-merge(x=df_out,y=weights,
            by.x=c("Node", "Factor"), by.y=c("node", "factor"))
  df_out$meta_node <- ifelse(df_out$tf == TRUE, "sink", ifelse(df_out$receptor_ligands == TRUE, "source", NA))
  if (Abs) {
    df_out$value <- abs(df_out$value)
  }
  return(df_out[!is.na(df_out$meta_node),])
}

#get nodes from all graphs
sink_source_egdeList<-list()
for (factor in factorS) {
  out<-list()
  out$up<-add_mofa_weight(factor_nodes[[factor]]$up, weights, factor, Abs = Abs)
  out$down<-add_mofa_weight(factor_nodes[[factor]]$down, weights, factor, Abs = Abs)
  sink_source_egdeList[[factor]]<-out
}


#function that makes a source and sink node and weights them based on their factor in arid1a, or TF activity. 
make_source_sink_graph<-function(df,
                                 tf_activity_df,
                                 Abs=T){

  # Create an empty graph
  g_source_sink <- make_empty_graph()
  # Add source and sink nodes
  g_source_sink <- add_vertices(g_source_sink, name = c("source", "sink"), nv = 2,
                                attr=list(nodeType=c("meta","meta")))
  #add tfs 
  tf_list<-unique(df$Node[df$tf])
  g_source_sink <- add_vertices(g_source_sink, name = tf_list, nv = length(tf_list),
                                attr=list(nodeType=rep("tf", length(tf_list))))
  #..and receptors
  receptors_list<-unique(df$Node[df$receptor_ligands])
  receptors_df<-df[df$Node %in% receptors_list,]
  receptors_df$capacity <- receptors_df$value
  g_source_sink <- add_vertices(g_source_sink, name = receptors_list, nv = length(receptors_list),
                                attr=list(nodeType=rep("receptors", length(receptors_list))))
  
  tf_df<-tf_activity_df[tf_activity_df$label %in% tf_list,]
  #absolute values.
  if(Abs){
    tf_df$capacity <- abs(tf_df$Activity)
  }else{
    tf_df$capacity <- tf_df$Activity
  }
  # Add edges from TFs to "sink"
  for (tfs in tf_list) {
    g_source_sink <- add_edges(g_source_sink, edges = c(tfs, "sink"),
                               attr = list(capacity=sum(tf_df$capacity[tf_df$TF==tfs])*1,
                                           source="dummy"))
  }
  # Add edges from "source" to receptors
  for (receptor in receptors_list) {
    g_source_sink <- add_edges(g_source_sink, edges = c("source", receptor),
                               attr = list(capacity=sum(receptors_df$capacity[receptors_df$Node==receptor]),
                                           source="dummy"))
  }
  # Plot the graph
  return(g_source_sink)
}

#make source and sink graphs and then merge them with the graphs derived from phuego. 
# add the sink (transcription factors) and source (receptor and ligands) nodes 


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

max_flow_graphs<-list()
for (factor in factorS) {
  df_in<-sink_source_egdeList[["Factor3"]]
  phuego_graphs<-factor_graphs[[factor]]
  
  #make graphs that define the inputs and outputs of the maxflow
  source_sink_up  <- make_source_sink_graph(df_in$up, tf_activity_df = arid1a_tf, Abs = Abs)
  source_sink_down  <- make_source_sink_graph(df_in$down, tf_activity_df = arid1a_tf, Abs = Abs)
  
  #combine these with the phuego graphs. 
  union_graph_up <- union2(phuego_graphs$up, source_sink_up)
  union_graph_down <- union2(phuego_graphs$down, source_sink_down)
  
  #add them to various lists. 
  max_flow_graphs[[factor]]$up <- union_graph_up
  max_flow_graphs[[factor]]$down <- union_graph_down
}


