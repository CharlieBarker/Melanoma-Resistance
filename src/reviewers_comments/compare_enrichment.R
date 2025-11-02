library(MOFA2)
library(dplyr)
library(EnsDb.Hsapiens.v86)
library(purrr)
library(stringr)
library(gprofiler2)
library(clusterProfiler)
library(tidyr)
library(igraph)
library(ggpubr)
library(ggrepel)
setwd("~/Desktop/Melanoma_Resistance//")

MOFAobject.trained<-load_model(file = "./results/mofa/mofa_object.hdf5")

###enrichment of the weights ####

weights <- get_weights(MOFAobject.trained,
                       views = "all",
                       as.data.frame = TRUE
)

weights_clean <- weights %>%
  separate(feature, into = c("gene", "datatype"), sep = "_") %>%
  mutate(gene = sub(";.*", "", gene)) %>%
  select(-datatype) %>%
  group_by(gene, factor) %>%
  summarise(value = sum(value), .groups = 'drop')

# Compute exact quantiles for PRKD1 in each factor
prkd1_exact_quantiles <- weights_clean %>%
  group_by(factor) %>%
  mutate(
    quantile_value = ecdf(value)(value)
  ) %>%
  filter(gene == "JUN") %>%
  select(factor, gene, value, quantile_value)

print(prkd1_exact_quantiles)

# Compute quantile cutoffs for each factor
quantiles <- weights_clean %>%
  group_by(factor) %>%
  summarise(
    q05 = quantile(value, 0.1),
    q95 = quantile(value, 0.9),
    .groups = "drop"
  )

# Join the quantiles back to the main data
weights_with_quantiles <- weights_clean %>%
  left_join(quantiles, by = "factor")

# Filter for bottom 5% and top 5%
top_bottom_genes <- weights_with_quantiles %>%
  filter(value <= q05 | value >= q95) %>%
  mutate(
    quantile = case_when(
      value <= q05 ~ "bottom_5%",
      value >= q95 ~ "top_5%"
    )
  ) %>%
  select(gene, factor, value, quantile)


#get networks 

# Get phuego graphs
KDE <- "0.5"
factorS <- c("Factor1", "Factor2", "Factor3")
results_dir <- "./results/phuego/"
factor_graphs <- list()
factor_genes_names <- list()

factor_to_vis<-"Factor1"

# Process each factor
for (factor in factorS) {
  file_graphml_up <- file.path(results_dir, paste0("KDE_", KDE), factor, "KDE_increased.graphml")
  file_graphml_down <- file.path(results_dir, paste0("KDE_", KDE), factor, "KDE_decreased.graphml")
  
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

library(clusterProfiler)
library(org.Hs.eg.db)
library(dplyr)

run_enrichment <- function(gene_list, gene_universe) {
  enrichGO(
    gene          = gene_list,
    OrgDb         = org.Hs.eg.db,
    keyType       = "SYMBOL",
    universe      = gene_universe,
    ont           = "BP",
    pAdjustMethod = "BH",
    readable      = TRUE,
    pvalueCutoff  = 0.1
  ) %>% as.data.frame()
}

library(ggplot2)
library(ggrepel)

make_comparison_plot <- function(factor_name, universe, factor_genes_df, quantile_genes_df, 
                                 quantile_type = c("bottom_5%","top_5%"), sig.cutOff = 0.01) {
  
  # Genes from factor-specific list
  factor_genes <- factor_genes_df %>%
    filter(factor == factor_name) %>%
    pull(values)
  
  # Genes from quantile-specific list
  quantile_genes <- quantile_genes_df %>%
    filter(factor %in% factor_name, quantile %in% quantile_type) %>%
    pull(gene)
  
  # Run enrichment for both
  enrich_factor <- run_enrichment(factor_genes, universe)
  enrich_quant  <- run_enrichment(quantile_genes, universe)
  
  # Merge by GO Description
  merged <- inner_join(
    enrich_factor %>% select(Description, OR_factor = FoldEnrichment, p_adjust_factor = p.adjust),
    enrich_quant  %>% select(Description, OR_quantile = FoldEnrichment, p_adjust_quantile = p.adjust),
    by = "Description"
  )
  merged$sig.colour <- paste0(merged$p_adjust_factor < sig.cutOff, merged$p_adjust_quantile < sig.cutOff)
  # Plot
  
  # Use a named vector for aesthetic colors
  colour_map <- c(
    "TRUETRUE" = "#E41A1C",    # both significant - red
    "TRUEFALSE" = "#377EB8",   # factor significant only - blue
    "FALSETRUE" = "#4DAF4A",   # quantile significant only - green
    "FALSEFALSE" = "#999999"   # neither - grey
  )
  
  # Plot
  scatterPlot<-ggplot(merged, aes(x = OR_factor, y = OR_quantile, label = Description, colour = sig.colour)) +
    geom_point() +
    geom_abline(slope = 1, linetype = "dashed") +
    scale_colour_manual(values = colour_map, guide = "none") +  # Use custom colours & remove legend
    labs(
      title = paste("GO OR Comparison for", factor_name, "-", quantile_type),
      x = "OR from Network Propagation",
      y = "OR from MOFA output"
    )  +
    cowplot::theme_cowplot() +
    theme(
      plot.title = element_text(size = 12, face = "bold"),
      panel.border = element_rect(colour = "black", fill = NA, linewidth = 1)
    ) +
    grids(linetype = "dashed")

  # Filter and rank enriched terms
  enrich1 <- enrich_factor %>%
    filter(p.adjust < sig.cutOff) %>%
    arrange(p.adjust) %>%
    mutate(rank1 = row_number())
  
  enrich2 <- enrich_quant %>%
    filter(p.adjust < sig.cutOff) %>%
    arrange(p.adjust) %>%
    mutate(rank2 = row_number())
  
  merged_ranks <- inner_join(
    enrich1 %>% select(Description, rank1),
    enrich2 %>% select(Description, rank2),
    by = "Description"
  )
  merged_ranks <- merged_ranks %>%
    mutate(delta_rank = rank1 - rank2) %>%
    arrange(desc(abs(delta_rank)))  # largest movers
  
  # Add labels only for top 10 and bottom 10 delta_rank
  rank_labeled <- merged_ranks %>%
    arrange(delta_rank) %>%
    mutate(label = if_else(row_number() <= 10 | row_number() > (n() - 10), Description, NA_character_))
  
  rankPlot <- ggplot(rank_labeled, aes(x = delta_rank, y = reorder(Description, delta_rank))) +
    geom_segment(aes(xend = 0, yend = Description), color = "gray") +
    geom_point(aes(color = delta_rank > 0), size = 3) +
    geom_text(
      aes(label = label),
      size = 3,
      hjust = -0.1    # nudge text slightly to the right
    ) +
    labs(
      title = "Rank Change of GO Terms",
      x = "Rank Change (Network Propagation - MOFA)",
      y = "Rank"
    ) +
    theme_cowplot() +
    theme(
      plot.title = element_text(size = 12, face = "bold"),
      panel.border = element_rect(colour = "black", fill = NA, linewidth = 1),
      axis.text.y = element_blank(),  # hide y-axis term labels
      axis.ticks.y = element_blank()  # hide y-axis ticks
    ) +
    # dashed grid lines along x-axis only
    theme(panel.grid.major.x = element_line(linetype = "dashed"),
          panel.grid.minor.x = element_blank(),
          panel.grid.major.y = element_blank(),
          panel.grid.minor.y = element_blank())
  
  return(list(rank=rankPlot, scatter = scatterPlot))
}

# Replace with your actual universe (e.g., unique(weights_clean$gene))
geneUniverse <- unique(weights_clean$gene)

# Create plots for three factors
plot1 <- make_comparison_plot("Factor1", geneUniverse, factor_genes_names_df, top_bottom_genes)
plot2 <- make_comparison_plot("Factor2", geneUniverse, factor_genes_names_df, top_bottom_genes, sig.cutOff = 0.1)
plot3 <- make_comparison_plot("Factor3", geneUniverse, factor_genes_names_df, top_bottom_genes)

# Load patchwork
library(patchwork)
library(cowplot)
final_plot <- plot_grid(
  plot1$scatter,
  plot2$scatter,
  plot3$scatter,
  ncol = 1,
  align = "hv",     # Align both horizontal and vertical axes
  axis = "tblr",    # Align top/bottom/left/right
  rel_heights = c(1, 1, 1)
)

# Save to PDF
ggsave("~/Desktop/Melanoma_Resistance/paper/response_to_reviewers/GO_comparison_plots.pdf", 
       final_plot, width = 8, height = 10, units = "in")


final_plot <- plot_grid(
  plot1$rank,
  plot2$rank,
  plot3$rank,
  ncol = 3,
  align = "hv",     # Align both horizontal and vertical axes
  axis = "tblr"    # Align top/bottom/left/right
)



# Save to PDF
ggsave("~/Desktop/Melanoma_Resistance/paper/response_to_reviewers/GO_rank_plots.pdf", 
       final_plot, width = 20, height = 12, units = "in")

