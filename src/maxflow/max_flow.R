
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
library(stringr)
library(purrr)
library(ggrepel)
library(OmnipathR)
library(magrittr)
library(dplyr)
library(tidyr)
library(ggpubr)
library(wesanderson)

load("./results/maxflow/max_flow_graphs.Rdata")
factorS <- c("Factor1", "Factor2", "Factor3")

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
    #there is negative flow here and i dont know whats going on.  
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

capcity_for_phuego<-20
for (factor in factorS) {
  flow_graphs <- max_flow_graphs[[factor]]
  
  E(flow_graphs$down)$capacity[!is.na(E(flow_graphs$down)$n_ref)]<-capcity_for_phuego
  E(flow_graphs$up)$capacity[!is.na(E(flow_graphs$up)$n_ref)]<-capcity_for_phuego
  
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
full_max_flow_out_df<-full_max_flow_out_df[full_max_flow_out_df$Up_or_Down=="up",]

#normalise flow
full_max_flow_out_df <- full_max_flow_out_df %>%
  group_by(regulation, factor) %>%
  mutate(normalized_flow = flow)

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


layout<-"sugiyama"


# Function to create DataFrame from graph and max flow result
plot_max_flow <- function(df_in, factor_to_study, layout) {
  to_plot_df <- df_in %>%
    filter(normalized_flow != 0) %>%
    filter(factor == factor_to_study)
  to_plot_df$weight<-NULL
  graph <- as_tbl_graph(to_plot_df) |> 
    mutate(Centrality = centrality_degree(mode = 'in'))
  prior <- create_layout(graph, layout)
  prior$x[prior$name %in% c("sink", "source")]<-round(mean(prior$x[!prior$name %in% c("sink", "source")]))
  
  prior$y[prior$name=="source"]<-prior$y[prior$name=="source"] + 4
  prior$y[prior$name=="sink"]<-prior$y[prior$name=="sink"] - 4
  
  # plot using ggraph
  graph_Factor<-ggraph(prior) + 
    geom_edge_hive(aes(colour = factor(regulation),
                       width = (normalized_flow),
                       alpha = (normalized_flow)),
                   show.legend = FALSE) + 
    geom_node_point(aes(size = Centrality)) + 
    facet_edges(~regulation) + 
    theme_graph(base_family="sans", foreground = 'darkgrey', fg_text_colour = 'white') + 
    geom_node_text(aes(label = name), color = 'grey', 
                   size = 3, repel = T) + 
    scale_edge_color_manual(values = wes_palette("Darjeeling1"))
  return(graph_Factor)
}
graph_Factor1<-plot_max_flow(df_in = full_max_flow_out_df, factor_to_study = "Factor1", layout = layout)
graph_Factor2<-plot_max_flow(df_in = full_max_flow_out_df, factor_to_study = "Factor2", layout = layout)



pdf(# The directory you want to save the file in
  width = 18, # The width of the plot in inches
  height = 12,
  file = "./paper/plots/max_flow.pdf")
# Assuming pal is defined somewhere in your code

ggarrange(graph_Factor1, graph_Factor2, ncol = 1,  common.legend = T)


dev.off()