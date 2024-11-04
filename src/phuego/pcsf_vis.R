library(OmnipathR)
library(PCSF)
library(readr)
library(igraph)
library(dplyr)
library(tidyr)

names_for_subnet <-factor_centrality$ind[factor_centrality$rank < 50]

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



pdf(file = paste0("/Users/charliebarker/Desktop/Melanoma_Resistance/paper/networks/test.pdf"),
    width = 18, height = 15)
# Generate the plot
ggraph(subnet, layout = subnet_l) +
  geom_point(data = L_df[L_df$Gene_name != "Centroid_All_Factors",],
             aes(x = x, y = y, colour = protein_Factor1),
             alpha = 1, size = 10,
             stroke = 5) +
  geom_point(data = L_df[L_df$Gene_name != "Centroid_All_Factors",],
             aes(x = x, y = y, colour = mRNA_Factor1),
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
dev.off()
