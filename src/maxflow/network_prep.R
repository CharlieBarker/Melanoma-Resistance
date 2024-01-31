
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
                               directed = T)

# If you want to add additional attributes to the graph vertices or edges, you can do so like this:
# For example, adding the 'is_directed', 'is_stimulation', 'is_inhibition', and 'consensus_direction' attributes to edges
E(graph_tf)$is_directed <- ppi_tf_expanded$is_directed
E(graph_tf)$is_stimulation <- ppi_tf_expanded$is_stimulation
E(graph_tf)$is_inhibition <- ppi_tf_expanded$is_inhibition
E(graph_tf)$n_ref <- ppi_tf_expanded$n_references

E(graph_tf)$source <- "omnipath_tf_network"

tf_name<-AnnotationDbi::select(EnsDb.Hsapiens.v86, keys = unique(ppi_tf_expanded$target), keytype = "UNIPROTID", columns = "GENENAME")
expan_tf_name<-AnnotationDbi::select(EnsDb.Hsapiens.v86, keys = unique(ppi_tf_expanded$source), keytype = "UNIPROTID", columns = "GENENAME")
expan_tf<-unique(c(expan_tf_name$GENENAME, tf_name$GENENAME))


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
  #code code 
  #find the union of these two graphs
  reg_graph_union <- igraph::union(as.directed(graph, mode = "mutual"), tf_graph)
  #remove these graphs
  largest_contiguous_subgraph <- largest_component(reg_graph_union)
  E(largest_contiguous_subgraph)$capacity<-E(largest_contiguous_subgraph)$weight
  return(largest_contiguous_subgraph)
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
  nodes <- V(igraph_object)$name
  genename_df <- AnnotationDbi::select(EnsDb.Hsapiens.v86, keys = nodes, 
                                       keytype = "UNIPROTID", 
                                       columns = "GENENAME")
  genename_df <- genename_df[!duplicated(genename_df[, "UNIPROTID"]), ]
  
  # Match the names based on UNIPROTID
  name_mapping <- setNames(genename_df$GENENAME, genename_df$UNIPROTID)
  
  # Update vertex names using the name_mapping
  V(igraph_object)$name <- name_mapping[V(igraph_object)$name]
  
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
  receptors_df$capacity <- receptors_df$value / sum(receptors_df$value)
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
                               attr = list(capacity=sum(tf_df$capacity[tf_df$TF==tfs])*10,
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
  df_in<-sink_source_egdeList[[factor]]
  phuego_graphs<-factor_graphs[[factor]]
  
  source_sink_up  <- make_source_sink_graph(df_in$up, tf_activity_df = arid1a_tf, Abs = Abs)
  source_sink_down  <- make_source_sink_graph(df_in$down, tf_activity_df = arid1a_tf, Abs = Abs)
  
  union_graph_up <- union2(phuego_graphs$up, source_sink_up)
  union_graph_down <- union2(phuego_graphs$down, source_sink_down)
  max_flow_graphs[[factor]]$up <- union_graph_up
  max_flow_graphs[[factor]]$down <- union_graph_down
}

# Function to calculate max flow
calculate_max_flow <- function(graph, source_name, sink_name,
                               mode="absolute") {
  source_node <- which(V(graph)$name == source_name)
  sink_node <- which(V(graph)$name == sink_name)
  if(mode=="absolute"){
    max_flow_result <- max_flow(graph, source = source_node, target = sink_node)
  }
  if(mode=="directed"){
    graph_up_reg<-graph
    graph_down_reg<-graph
    sink_source_capacities_up<-E(graph)$capacity[is.na(E(graph)$source)]
    sink_source_capacities_down<-E(graph)$capacity[is.na(E(graph)$source)]
    
    sink_source_capacities_up[sink_source_capacities_up < 0] <- 0
    sink_source_capacities_down[sink_source_capacities_down > 0] <- 0
    sink_source_capacities_down<-abs(sink_source_capacities_down)
    
    E(graph_up_reg)$capacity[is.na(E(graph_up_reg)$source)]<-sink_source_capacities_up
    E(graph_down_reg)$capacity[is.na(E(graph_down_reg)$source)]<-sink_source_capacities_down
    
    max_flow_result_up <- max_flow(graph_up_reg, source = source_node, target = sink_node)
    max_flow_result_down <- max_flow(graph_down_reg, source = source_node, target = sink_node)
    max_flow_result<-list(up_reg=max_flow_result_up, down_reg=max_flow_result_down)
  }
  return(max_flow_result)
}

# Function to create DataFrame from graph and max flow result
create_df_from_graph <- function(graph, max_flow_result, factor, direction,
                                 mode="absolute") {
  df_graph <- igraph::as_data_frame(graph)
  if(mode=="absolute"){
    df_graph$flow <- max_flow_result$flow
  }
  if(mode=="directed"){
    df_graph_increased<-df_graph
    df_graph_increased$flow <- max_flow_result$up_reg$flow
    df_graph_increased$regulation <- "ARID1A KO, Increased signalling"
    
    df_graph_decreased<-df_graph
    df_graph_decreased$flow <- max_flow_result$down_reg$flow
    df_graph_decreased$regulation <- "ARID1A KO, Decreased signalling"
    df_graph <- rbind(df_graph_increased, df_graph_decreased)
  }
  df_graph$factor <- factor
  df_graph$Up_or_Down <- direction
  return(df_graph)
}

# Main loop
max_flow_out <- list()
max_flow_out_df <- list()

for (factor in factorS) {
  flow_graphs <- max_flow_graphs[[factor]]
  
  E(flow_graphs$down)$capacity[!is.na(E(flow_graphs$down)$n_ref)]<-1
  E(flow_graphs$up)$capacity[!is.na(E(flow_graphs$up)$n_ref)]<-1
  
  # Calculate max flow for 'up' graph
  out_up <- calculate_max_flow(flow_graphs$up, "source", "sink", mode = "directed")

  # Create DataFrame for 'up' graph
  df_graph_up <- create_df_from_graph(flow_graphs$up, out_up, factor, "up", mode = "directed")
  
  # Calculate max flow for 'down' graph
  out_down <- calculate_max_flow(flow_graphs$down, "source", "sink", mode = "directed")

  # Create DataFrame for 'down' graph
  df_graph_down <- create_df_from_graph(flow_graphs$down, out_down, factor, "down", mode = "directed")
  
  # Combine DataFrames
  max_flow_out_df[[factor]] <- rbind(df_graph_up, df_graph_down)
}

full_max_flow_out_df<-bind_rows(max_flow_out_df)





#visualisation 

library(ggraph)
#> Loading required package: ggplot2
library(tidygraph)
#> 
#> Attaching package: 'tidygraph'
#> The following object is masked from 'package:stats':
#> 
#>     filter

# Create graph of highschool friendships


# Assuming 'full_max_flow_out_df' is your data frame



to_plot_df <- full_max_flow_out_df %>%
  filter(flow != 0) %>%
  filter(factor == "Factor1")
to_plot_df$weight<-NULL
graph <- as_tbl_graph(to_plot_df) |> 
  mutate(Centrality = centrality_degree(mode = 'in'))
# plot using ggraph
graph_Factor1<-ggraph(graph, layout = 'kk') + 
  geom_edge_link(aes(colour = factor(regulation))) + 
  geom_edge_fan(aes(alpha = flow), show.legend = FALSE) + 
  geom_node_point(aes(size = Centrality)) + 
  facet_edges(~regulation) + 
  theme_graph(foreground = 'darkgrey', fg_text_colour = 'white') + 
  geom_node_text(aes(label = name), color = 'grey', 
                 size = 3)

to_plot_df <- full_max_flow_out_df %>%
  filter(flow != 0) %>%
  filter(factor == "Factor3")
to_plot_df$weight<-NULL
graph <- as_tbl_graph(to_plot_df) |> 
  mutate(Centrality = centrality_degree(mode = 'in'))
# plot using ggraph
graph_Factor2<-ggraph(graph, layout = 'kk') + 
  geom_edge_link(aes(colour = factor(regulation))) + 
  geom_edge_fan(aes(alpha = flow), show.legend = FALSE) + 
  geom_node_point(aes(size = Centrality)) + 
  facet_edges(~regulation) + 
  theme_graph(foreground = 'darkgrey', fg_text_colour = 'white') + 
  geom_node_text(aes(label = name), color = 'grey', 
                 size = 3)
ggarrange(graph_Factor1, graph_Factor2, ncol = 1, nrow = 2)
