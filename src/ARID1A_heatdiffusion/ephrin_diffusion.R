
#find correct environment
packLib="/usr/lib/R"
if (file.exists(packLib)) {
  reticulate::use_condaenv("py37", required = T)
  setwd(dir = "~/Desktop/Melanoma_Resistance/")
}else {
  setwd(dir = "~/Desktop/Melanoma_Resistance/")
}
load('./results/heatdiffusion/data_for_heat_diffusion.Rdata')

library(readxl)
library(ggplot2)
library(dplyr)
library(ggpubr)
library(diffusr)
library(ggrepel)
library(wesanderson)
library(igraph)
library(purrr)
library(stringr)
library(cowplot)
library(dplyr)
library(tidyr)
library(EnsDb.Hsapiens.v86)


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

perform_random_walk <- function(graph, start_vector, num_permutations = 100, correct.for.hubs = T) {
  # Calculate the stationary distribution for the original graph, suppressing messages
  adj_matrix <- as_adjacency_matrix(graph)
  pt_original <- suppressMessages(random.walk(start_vector, data.matrix(adj_matrix),
                                              r = 0.02,
                                              correct.for.hubs = correct.for.hubs,
                                              do.analytical = TRUE))
  stationary_distribution_original <- pt_original$p.inf  # Stationary distribution for each node

  # Get the degree distribution of the original graph
  degree_distribution <- degree(graph)

  # Initialize a matrix to hold stationary distributions from random graphs (each row is a node, each column is a permutation)
  num_nodes <- length(stationary_distribution_original)
  stationary_distributions_permuted <- matrix(0, nrow = num_nodes, ncol = num_permutations)

  # Create a text progress bar
  pb <- txtProgressBar(min = 0, max = num_permutations, style = 3)

  for (i in 1:num_permutations) {
    repeat {
      # Create a random graph that preserves the degree distribution
      random_graph <- sample_degseq(degree_distribution, method = "simple")

      # Check if the random graph is connected (ergodic)
      if (is_connected(random_graph)) {
        # Calculate the stationary distribution for the random graph, suppressing messages
        adj_matrix_random <- as_adjacency_matrix(random_graph)
        pt_random <- suppressMessages(random.walk(start_vector, data.matrix(adj_matrix_random),
                                                  r = 0.02,
                                                  correct.for.hubs = correct.for.hubs,
                                                  do.analytical = TRUE))

        # Store the stationary distribution of the random graph for all nodes in the matrix
        stationary_distributions_permuted[, i] <- pt_random$p.inf

        # Update the progress bar
        setTxtProgressBar(pb, i)

        # Break the repeat loop since we have an ergodic graph
        break
      }
      # If the graph is not connected, it will repeat the loop to create a new random graph
    }
  }

  # Close the progress bar
  close(pb)

  # Calculate p-values for each node: the proportion of times the permuted stationary distribution for a node is >= the original one
  p_values <- sapply(1:num_nodes, function(node_idx) {
    mean(stationary_distributions_permuted[node_idx, ] >= stationary_distribution_original[node_idx])
  })

  # Return the original stationary distribution for each node and the node-specific p-values
  return(list(probability = stationary_distribution_original, p_values = p_values))
}



kinomics_files<-list.files(path = "./data/kinomics//", full.names = T, recursive = T)

results_list<-list()
for (UKA_file in kinomics_files[grep(kinomics_files, pattern = "vs")]) {
  kinomics_results<-read_xlsx(path =UKA_file)
  results_list[[UKA_file]]<-kinomics_results
}
# Apply the conversion to each data frame in the list
results_list <- lapply(results_list, function(x) {
  x$`SD Kinase Statitistic` <- as.numeric(x$`SD Kinase Statitistic`)
  return(x)
})

# Bind the rows
complete_results <- bind_rows(results_list, .id = "file")
complete_results$experiment<-unlist(map(str_split(complete_results$file, pattern = "/"), last))
complete_results$background <- unlist(map(str_split(complete_results$file, pattern = "/"), 6))

# Assuming `complete_results` is your data frame containing the required columns

# Reorder Kinase Name by Mean Kinase Statistic
complete_results$`Kinase Name` <- reorder(complete_results$`Kinase Name`, complete_results$`Mean Kinase Statistic`)

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

get_nodes<-function(igraph_object, conv=F){
  nodes<-V(igraph_object)$name
  weights<-igraph::page_rank(igraph_object)

  if (conv) {
    genename_df <- AnnotationDbi::select(EnsDb.Hsapiens.v86, keys = nodes,
                                         keytype = "UNIPROTID",
                                         columns = "GENENAME")
    return(sort(genename_df$GENENAME))
  }
  else {
    return(sort(weights$vector))
  }
}

#get nodes from all graphs
factor_nodes <- lapply(factor_graphs, function(sublist) {
  sapply(sublist, get_nodes, T)
})

complete_results$arid1a_status<-unlist(map(str_split(complete_results$background, pattern = "__"),2))
complete_results$kinase_type<-unlist(map(str_split(complete_results$background, pattern = "__"),1))

# Transform the dataset
experiment_matrix <- complete_results %>%
  # Spread the data into a matrix format
  dplyr::select(`Kinase Uniprot ID`, file, `Median Kinase Statistic`) %>%
  pivot_wider(names_from = file, values_from = `Median Kinase Statistic`)
genename_df <- AnnotationDbi::select(EnsDb.Hsapiens.v86, keys = experiment_matrix$`Kinase Uniprot ID`, keytype = "UNIPROTID", columns = "GENENAME");
experiment_matrix$`Kinase Name` <- genename_df$GENENAME[match(experiment_matrix$`Kinase Uniprot ID`, genename_df$UNIPROTID)]

# Function to simplify the column names
simplify_column_names <- function(colname) {
  colname %>%
    # Remove the initial path "./data/kinomics///"
    gsub("^.*/kinomics///", "", .) %>%
    # Remove the ".xlsx" extension
    gsub("\\.xlsx", "", .) %>%
    # Replace "__" with " - " for clarity
    gsub("__", " - ", .)
}

# Apply the function to column names
simplified_colnames <- sapply(colnames(experiment_matrix), simplify_column_names)

# Set the new column names
colnames(experiment_matrix) <- simplified_colnames
#get list of receptors
receptors <- OmnipathR::import_omnipath_intercell(
    parent = 'receptor',
    topology = 'pmtm',
    consensus_percentile = 50,
    loc_consensus_percentile = 30,
    entity_types = 'protein'
  ) %>%
  pull(genesymbol) %>% unique()

# Filter the data based on Kinase Name matching Gene_name from graph 'g'
gene_affected <- experiment_matrix %>%
  dplyr::filter(`Kinase Name` %in% unname(unlist(factor_nodes)))

# Further filter data for receptors
receptors_affected <- gene_affected %>%
  dplyr::filter(`Kinase Name` %in% receptors)

# Reshape the data to long format excluding non-numeric/character columns like `Kinase Uniprot ID`
molten_receptors_affected <- receptors_affected %>%
  pivot_longer(cols = -c(`Kinase Name`, `Kinase Uniprot ID`), # Specify the columns to exclude
               names_to = "variable", values_to = "value") %>%
  dplyr::filter(!is.na(value)) %>%
  dplyr::filter(!str_detect(variable, "vs WT"))

# Split the 'variable' column into 'genetic_background' and 'variable',
# and clean up the string patterns
molten_receptors_affected <- molten_receptors_affected %>%
  separate(variable, into = c("genetic_background", "variable"), sep = "/") %>%
  mutate(genetic_background = str_remove_all(genetic_background, "PTK - |STK - "),
         variable = str_remove_all(variable, "ARID1A "))

# Filter rows where the 'variable' column contains 'Combined vs Untreated'
# and the absolute value of 'value' is greater than 1
filtered_receptors_affected <- molten_receptors_affected %>%
  dplyr::filter(str_detect(variable, "Combined vs Untreated"),
                abs(value) > 1)

filtered_receptors_affected <- filtered_receptors_affected[filtered_receptors_affected$genetic_background == "ARID1A_KO",]

# Extract and normalize kinases dysregulated in ARID1A
arid1a_rtks <- filtered_receptors_affected$value[match(V(g)$Gene_name, filtered_receptors_affected$`Kinase Name`)]
arid1a_rtks[is.na(arid1a_rtks) | arid1a_rtks < 0] <- 0
arid1a_kinome_receptors <- normalize_vector(arid1a_rtks)

# Perform random walks and get stationary distributions
rwr_results <- perform_random_walk(g, arid1a_kinome_receptors, num_permutations = 10000)
V(g)$dist <- rwr_results$probability
V(g)$p_values <- rwr_results$p_values + 0.00001

ephrin_heat<-data.frame(name=V(g)$Gene_name,
                        uniprot_id=V(g)$name,
                        dist=as.numeric(V(g)$dist),
                        p_values = V(g)$p_values,
                        in_seed= V(g)$Gene_name %in% filtered_receptors_affected$`Kinase Name`)

# Order the data by 'dist'
ephrin_heat <- ephrin_heat %>%
  arrange(desc(dist))
ephrin_heat<-ephrin_heat[ephrin_heat$p_values < 0.1,]
top_non_seed_nodes <- ephrin_heat[ephrin_heat$in_seed == FALSE, ] %>%
  arrange(desc(dist)) %>%
  head(15)
# Pivot wider to create columns for each genetic background and spread the values
wide_receptors_affected <- molten_receptors_affected %>%
  pivot_wider(names_from = genetic_background, values_from = value)

#look at the levels of proteins of interest

all_abundace<-list(
  protein=read.csv(file = "./data/input_data/proteins.csv"),
  rna=read.csv(file = "./data/input_data/rna_expression.csv")
)

top_non_seed_list<-as.list(top_non_seed_nodes$uniprot_id)
names(top_non_seed_list)<-top_non_seed_nodes$name

# Define the drugs and their corresponding colors
drug_colors <- c(
  "Untreated" = "#A2AEBB",
  "Vemurafenib" = "#FFBA08",
  "Trametinib" = "#D00000",
  "Combination" = "#3F88C5"
)

proteins_df<-reshape2::melt(all_abundace$protein[all_abundace$protein$X %in% unname(top_non_seed_list),])
rna_df<-reshape2::melt(all_abundace$rna[all_abundace$rna$X %in% names(top_non_seed_list),])
replacement_Vec<-c("Untreated","Vermurafenib_1uM","Trametinib_10nM","vemurafenib_and_trametinib")
names(replacement_Vec)<- c("Untreated", "Vemurafenib", "Trametinib", "Combination")

split_exp_conditions<-function(df, error=F){
  df$drug<-unlist(map(str_split(df$variable, pattern = "__"),1))
  df$ko<-unlist(map(str_split(df$variable, pattern = "__"),2))
  df$ko<-gsub("\\..*","",df$ko)
  df$variable<-NULL
  df$drug<-names(replacement_Vec)[match(df$drug, unname(replacement_Vec))]
  # Define the desired order of levels
  desired_order <- c("Untreated", "Vemurafenib", "Trametinib", "Combination")
  # Convert my_column to a factor with the specified order
  df$drug <- factor(df$drug, levels = desired_order)
  if (error) {
    df <- df %>%
      group_by(X, drug, ko) %>%
      summarise(
        n=n(),
        mean=mean(value),
        sd=sd(value)
      ) %>%
      mutate( se=sd/sqrt(n))  %>%
      mutate( ic=se * qt((1-0.05)/2 + .5, n-1))

  }

  return(df)
}


proteins_df<-split_exp_conditions(proteins_df)
rna_df<-split_exp_conditions(rna_df)
rna_df$X <- factor(rna_df$X, levels = names(top_non_seed_list))

proteins_df$X<- names(top_non_seed_list)[match(proteins_df$X, unname(top_non_seed_list))]

complete_df<-rbind(data.frame(proteins_df, data="Protein abundace"),
                   data.frame(rna_df, data="mRNA abundace"))


library(ggraph)

pdf(file = "~/Desktop/Melanoma_Resistance/results/heatdiffusion/ephrin_activity_diffusion.pdf",   # The directory you want to save the file in
    width = 10,  # The width of the plot in inches
    height = 10) # The height of the plot in inches

ephrin_heat$label=""
ephrin_heat$label[ephrin_heat$dist>0.005] = ephrin_heat$name[ephrin_heat$dist>0.005]

ggplot(ephrin_heat[ephrin_heat$in_seed==F,], aes(x = -log10(p_values), y = dist, size = -log10(p_values), color = in_seed)) +
  geom_point() +  # Points with specified fill color
  labs(title = "Random Walk from activated Ephrin receptors",
       x = "-Log10(P-Values)",  # X-axis label
       y = "Stationary Probability") +  # Y-axis label
  scale_color_manual(values = c("TRUE" = "darkorange", "FALSE" = "#446455ff"), name = "In Seed") +  # Custom colors
  geom_text_repel(aes(label = label), colour = "black", size = 8) +  # Add labels with black color
  theme_cowplot() +
  grids(linetype = "dashed") +
  theme(axis.text.y = element_text(size = 12),  # Adjust y-axis text size
        axis.text.x = element_text(size = 12, angle = 45, hjust = 1),  # Tilt x-axis labels
        plot.title = element_text(size = 16, face = "bold"),
        panel.border = element_rect(colour = "black", fill = NA, linewidth = 1)) +
  guides(size = guide_legend(title = "-Log10(P-Values)"), color = guide_legend(title = "In Seed"))

dev.off()

pdf(file = "~/Desktop/Melanoma_Resistance/results/heatdiffusion/ephrin_heat.pdf",   # The directory you want to save the file in
    width = 14,  # The width of the plot in inches
    height = 18) # The height of the plot in inches

ggplot(wide_receptors_affected, aes(x = ARID1A_KO, y = WT)) +
  geom_point(colour="lightblue") +
  cowplot::theme_cowplot() +
  theme(
    plot.title = element_text(size = 15, face = "bold"),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 1),
  ) +
  grids(linetype = "dashed") +
  xlab(paste("Kinase Activity in ARID1A")) +
  ylab("Kinase Activity in WT") +
  facet_wrap(~variable, nrow = 2) +
  # Add axis lines
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black") +

  # Add text labels for points where values exceed thresholds
  geom_text_repel(
    aes(label =
          `Kinase Name`),
    size = 3,  # Adjust text size as needed
    box.padding = 0.35,  # Padding around the text label to avoid overlap
    point.padding = 0.5,  # Space between the label and the point
    max.overlaps = Inf,  # Avoid limit on label overlaps
    min.segment.length = 0  # Draw segment lines for all labels
  )

complete_df %>%
  ggplot(aes(y = value, x = ko, fill = drug)) +
  geom_boxplot() +
  geom_jitter(color = "black", size = 0.4, alpha = 0.9) +
  scale_fill_manual(values = drug_colors) +
  ggtitle("mRNA and protein abundances of top 15 nodes after Ephrin diffusion") +
  xlab("Gene/Protein") +  # Change x-axis label to "Abundance"
  ylab("Abundance") +  # Change y-axis label to "Gene/Protein"
  cowplot::theme_cowplot() +
  facet_wrap(X ~ data, scales="free") +
  grids(linetype = "dashed") +
  theme(
    plot.title = element_text(size = 20, face = "bold"),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 1)  # Use linewidth instead of size
  )


# Step 1: Find the vertex corresponding to "NCK1"
nck1_vertex <- V(g)[V(g)$Gene_name == "FYN"]

# Step 2: Extract the ego network (1-step neighbors)
nck1_ego <- ego(g, order = 1, nodes = nck1_vertex, mode = "all")

# Step 3: Create a subgraph from the ego network
nck1_subgraph <- induced_subgraph(g, unlist(nck1_ego))
ggraph(nck1_subgraph, layout = "fr") +  # Using Fruchterman-Reingold layout
  geom_edge_link() +  # Draw edges
  geom_node_label(aes(label = Gene_name), repel = TRUE) +  # Add labels
  theme_minimal() +  # Minimal theme
  labs(title = "Ego Network of NCK1")  # Title for the plot

# Step 1: Find the vertex corresponding to "NCK1"
nck1_vertex <- V(g)[V(g)$Gene_name == "DUSP18"]

# Step 2: Extract the ego network (1-step neighbors)
nck1_ego <- ego(g, order = 1, nodes = nck1_vertex, mode = "all")

# Step 3: Create a subgraph from the ego network
nck1_subgraph <- induced_subgraph(g, unlist(nck1_ego))
ggraph(nck1_subgraph, layout = "fr") +  # Using Fruchterman-Reingold layout
  geom_edge_link() +  # Draw edges
  geom_node_label(aes(label = Gene_name), repel = TRUE) +  # Add labels
  theme_minimal() +  # Minimal theme
  labs(title = "Ego Network of NCK1")  # Title for the plot


dev.off()




