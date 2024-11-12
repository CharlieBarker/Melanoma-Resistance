load('./results/heatdiffusion/data_for_heat_diffusion.Rdata')
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
library(diffusr)
library(OmnipathR)

source("./src/functions/pathwayLayout.R")


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

# basic scatterplot
L_df$is_receptor = L_df$Gene_name %in% receptors


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
pdf(file = paste0("./results/heatdiffusion/individual_receptor_heat.pdf"),
    width = 9, height = 3)

library(ComplexHeatmap)
library(viridis)

# Example usage (you need to provide the actual objects for `arid1a_up_receptors`, `subnet`, and `vector_to_binary_matrix`):
to_test<-arid1a_up_receptors
df_arid1a_up<-data.frame(L_df$Gene_name, perform_random_walk(subnet, vector_to_binary_matrix(to_test)))
rtks_list<-L_df$Gene_name[as.logical(rowSums(vector_to_binary_matrix(to_test)))]
arid1a_up_heatmap <- create_heatmap(data_frame = df_arid1a_up,
                          to_test = arid1a_up_receptors,
                          nodes_to_see = c("JUN", "PRKD1", "PTPN6", "FYN"),
                          subnet = subnet,
                          rtks_list = rtks_list,
                          plot_title="Receptors upregulated after ARID1A KO")
to_test<-combined_up_receptors
df_combined_up<-data.frame(L_df$Gene_name, perform_random_walk(subnet, vector_to_binary_matrix(to_test)))
rtks_list<-L_df$Gene_name[as.logical(rowSums(vector_to_binary_matrix(to_test)))]
combined_heatmap <- create_heatmap(data_frame = df_combined_up,
                          to_test = combined_up_receptors,
                          nodes_to_see = c("JUN", "PRKD1","PTPN6", "FYN"),
                          subnet = subnet,
                          rtks_list = rtks_list,
                          plot_title="Receptors upregulated after combined therapy")

ht_list = arid1a_up_heatmap + combined_heatmap
draw(ht_list)

dev.off()
