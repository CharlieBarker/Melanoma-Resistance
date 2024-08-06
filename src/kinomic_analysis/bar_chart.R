
#find correct environment 
packLib="/usr/lib/R"
if (file.exists(packLib)) {
  reticulate::use_condaenv("py37", required = T)
  setwd(dir = "~/Desktop/Melanoma_Resistance/")
}else {
  setwd(dir = "~/Desktop/Melanoma_Resistance/")
}

library(readxl)
library(ggplot2)
library(dplyr)
library(ggpubr)
library(ggrepel)
library(wesanderson)
library(igraph)
library(cowplot)

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
  sapply(sublist, get_nodes)
})


# Assuming kinase_list is a vector of kinase names you want to mark with stars
kinase_list_all_factors <- list(Factor1=stack(factor_nodes$Factor1), #factor 1 describes changes in all drugs
                                Factor2=stack(factor_nodes$Factor2), #factor 2 describes changes specific to combination 
                                Factor3=stack(factor_nodes$Factor3) #factor 3 describes ARID1A
                    )
#this one works, we know at least

complete_results$arid1a_status<-unlist(map(str_split(complete_results$background, pattern = "__"),2))
complete_results$kinase_type<-unlist(map(str_split(complete_results$background, pattern = "__"),1))
# Get unique experiment names
experiment_names <- unique(complete_results$experiment)

# Remove "ARID1A " and ".xlsx" from experiment names
clean_experiment_names <- gsub("^ARID1A ", "", experiment_names)
clean_experiment_names <- gsub(".xlsx$", "", clean_experiment_names)

# Update the experiment names in the complete_results dataframe
complete_results$experiment <- gsub("^ARID1A ", "", complete_results$experiment)
complete_results$experiment <- gsub(".xlsx$", "", complete_results$experiment)

complete_results$experiment <- factor(complete_results$experiment, levels = c("Vemurafenib vs Untreated", "Trametinib vs Untreated", "Combined vs Untreated",
                                                                            "Combined vs Vemurafenib", "Combined vs Trametinib"))

drug_targets<-complete_results[complete_results$`Kinase Name` %in% c("ERK1", "BRAF", "PKD1", "JNK1", "JNK2", "JNK3"),]
drug_targets<-drug_targets[!grepl(drug_targets$experiment, pattern="ARID1A"),]
pal <- wes_palette("Zissou1", 100, type = "continuous")


ggplot(drug_targets[!is.na(drug_targets$experiment),], aes(x=`Kinase Name`, y=`Median Kinase Statistic`, 
                         colour=`Median Kinase Statistic`)) +
  scale_colour_gradientn(colours = pal) + 
  geom_hline(yintercept = 0, color = "black") + # Add horizontal line at y=0
  geom_segment(aes(x=`Kinase Name`, xend=`Kinase Name`, y=0, yend=`Median Kinase Statistic`), color="grey") +
  geom_point(aes(size = `Mean Specificity Score`)) +  # Adjust the size of points here
  theme_light() +
  theme(
    panel.grid.major.x = element_blank(),
    panel.border = element_rect(color = "lightgrey", fill = NA, size = 1),
    axis.ticks.x = element_blank()
  ) +
  xlab("") +
  ylab("Median Kinase Statistic") +
  facet_grid(arid1a_status ~ experiment)

ggsave("./results/kinomics_microarray/erk_braf_lolipop.pdf", width = 40, height = 12, units = "cm")

kinase_list<-kinase_list_all_factors$Factor1

# Create the ggplot with ordered Kinase Name and stars for specific kinases
subset_results <- complete_results %>%
  dplyr::filter(arid1a_status == "WT") %>%
  dplyr::filter(kinase_type == "STK") %>%
  mutate(
    in_network = `Kinase Uniprot ID` %in% rownames(kinase_list),
    which_network = case_when(
      `Kinase Uniprot ID` %in% rownames(kinase_list)[kinase_list$ind == "up"] ~ "up regulated network",
      `Kinase Uniprot ID` %in% rownames(kinase_list)[kinase_list$ind == "down"] ~ "down regulated network",
      TRUE ~ "outside network"
    )
  )

subset_results$centrality <- kinase_list$values[match(subset_results$`Kinase Uniprot ID`, rownames(kinase_list))]
#plot kinomics over centrality 

# basic scatterplot with a line of best fit
# Install ggpubr if it's not already installed
# install.packages("ggpubr")

ggplot(subset_results[!is.na(subset_results$centrality),],
       aes(y = abs(`Median Kinase Statistic`), 
           x = centrality, 
           colour=`Median Kinase Statistic`, size = `Mean Specificity Score`)) +
  scale_colour_gradientn(colours = pal) + 
  geom_hline(yintercept = 0, color = "darkgrey", linetype="dashed") + # Add horizontal line at y=0  
  geom_smooth(method = "lm", se = T, size = 1, color="darkgrey") + # Add a line of best fit using linear regression
  geom_point() + 
  facet_wrap( ~ experiment, ncol = 2) + 
  theme_cowplot() +
  theme(legend.position = "bottom") + # Placing legend at the bottom
  geom_text_repel(aes(label = `Kinase Name`), colour = "black", size = 4, max.overlaps = 10) + # Apply label aesthetic here
  coord_cartesian(ylim = c(0, NA)) + # Set the lower limit of y-axis to 0
  stat_cor(method = "pearson", aes(label = paste(..r.label.., ..p.label.., sep = "~`,`~")), label.x = 0, label.y = 1.5) + # Add Pearson correlation coefficient
  labs(x = "Centrality (PageRank)", y = "Absolute Median Kinase Statistic") # Add axis titles
ggsave("./results/kinomics_microarray/factor1_scatter.pdf", width = 25, height = 20, units = "cm")




# Assuming factor_nodes is your list
factor1_ids <- rownames(stack(factor_nodes$Factor1))
factor2_ids <- rownames(stack(factor_nodes$Factor2))

# Classify Kinase Uniprot IDs
to_plot_box <- complete_results %>%
  filter(arid1a_status == "WT") %>%
  mutate(
    factor_network = case_when(
      `Kinase Uniprot ID` %in% factor1_ids & `Kinase Uniprot ID` %in% factor2_ids ~ "both",
      `Kinase Uniprot ID` %in% factor1_ids ~ "Factor1",
      `Kinase Uniprot ID` %in% factor2_ids ~ "Factor2",
      TRUE ~ "neither"
    )
  )

# Plot
p <- to_plot_box[to_plot_box$factor_network != "both",] %>%
  ggplot(aes(x=factor_network, y=abs(`Median Kinase Statistic`), fill=factor_network)) +
  geom_boxplot() +
  geom_jitter(color="black", size=0.4, alpha=0.9) +
  theme(
    legend.position="none",
    plot.title = element_text(size=11)
  ) +
  xlab("") + facet_grid(kinase_type ~ experiment, scales = "free") + theme_cowplot()+
  scale_fill_manual(values = wes_palette("FantasticFox1"))

# Add significance testing
p + stat_compare_means(method = "t.test", comparisons = list(c("Factor1", "Factor2"),c("Factor1", "neither"),c("neither", "Factor2")), label = "p.format")
ggsave("./results/kinomics_microarray/box_plot_network.pdf", width = 40, height = 20, units = "cm")




plot <- ggplot(subset_results) +
  geom_bar(aes(x = reorder(`Kinase Name`, `Mean Kinase Statistic`), y = `Median Kinase Statistic`), 
           stat = "identity", fill = "skyblue", alpha = 0.7) +
  geom_errorbar(aes(x = reorder(`Kinase Name`, `Mean Kinase Statistic`), 
                    ymin = `Median Kinase Statistic` - `SD Kinase Statitistic`, 
                    ymax = `Median Kinase Statistic` + `SD Kinase Statitistic`), 
                width = 0.4, colour = "orange", alpha = 0.9, size = 1.3) + 
  facet_wrap(~ experiment, scales = "free_x", nrow = 5) +
  labs(x = "Kinase Name", y = "Median Kinase Statistic") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

# Add stars to specific kinases in the list kinase_list
plot<-plot + 
  geom_text(data = subset(subset_results, `Kinase Uniprot ID` %in% kinase_list),
            aes(x = reorder(`Kinase Name`, `Mean Kinase Statistic`), 
                y = `Median Kinase Statistic` + (`SD Kinase Statitistic` * 1.2), 
                label = "*"), 
            vjust = 0, size = 3, color = "red", fontface = "bold")



# Save the plot to a PDF file
ggsave("./paper/plots/kinomics_barplot.pdf", plot, width = 30, height = 30, units = "in")



ARID1A_KO<-complete_results[complete_results$arid1a_status == "ARID1A_KO" & complete_results$kinase_type == "STK",]
plot <- ggplot(ARID1A_KO) +
  geom_bar(aes(x = reorder(`Kinase Name`, `Mean Kinase Statistic`), y = `Median Kinase Statistic`), 
           stat = "identity", fill = "skyblue", alpha = 0.7) +
  geom_errorbar(aes(x = reorder(`Kinase Name`, `Mean Kinase Statistic`), 
                    ymin = `Median Kinase Statistic` - `SD Kinase Statitistic`, 
                    ymax = `Median Kinase Statistic` + `SD Kinase Statitistic`), 
                width = 0.4, colour = "orange", alpha = 0.9, size = 1.3) + 
  facet_wrap( ~ experiment, scales = "free") +
  labs(x = "Kinase Name", y = "Median Kinase Statistic") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

# Add stars to specific kinases in the list kinase_list
plot<-plot + 
  geom_text(data = subset(ARID1A_KO, `Kinase Uniprot ID` %in% kinase_list),
            aes(x = reorder(`Kinase Name`, `Mean Kinase Statistic`), 
                y = `Median Kinase Statistic` + (`SD Kinase Statitistic` * 1.2), 
                label = "*"), 
            vjust = 0, size = 3, color = "red", fontface = "bold")
# Save the plot to a PDF file
ggsave("./paper/plots/ARID1A_KO_barplot.pdf", plot, width = 30, height = 30, units = "in")



# Most basic error bar
ggplot(ARID1A_KO) +
  geom_bar( aes(x=name, y=value), stat="identity", fill="skyblue", alpha=0.7) +
  geom_errorbar( aes(x=name, ymin=value-sd, ymax=value+sd), width=0.4, colour="orange", alpha=0.9, size=1.3)



library(tidyr)

# Assuming `complete_results` is your original data frame

# Pivot the data frame to wide format based on experiment
wide_results <- pivot_wider(
  data = complete_results,
  id_cols = c(`Kinase Uniprot ID`, `Kinase Name`),
  names_from = experiment,
  values_from = `Mean Kinase Statistic`,
  names_prefix = "Mean Kinase Statistic "
)
# 
# 
# library(tidyr)
# library(ggplot2)
# library(dplyr)
# library(ggpubr)
# library(GGally)
# library(purrr)
# library(ghibli)
# library(patchwork)
# library(ggExtra)
# 
# ghibli_cs<-"MononokeMedium"
# 
# nCOL<-ncol(wide_results)
# ggpairs(wide_results,                 # Data frame
#         columns = 3:nCOL)
# 
# 
# g1 <- ggplot(df2, aes(Encorafinib, Trametinib, colour = annot)) +
#   geom_point() + cowplot::theme_cowplot() + 
#   geom_hline(yintercept=0, linetype='dotted', col = 'darkred')+
#   geom_vline(xintercept=0, linetype='dotted', col = 'darkred')+
#   scale_color_ghibli_d("YesterdayMedium", direction = -1)+ 
#   theme(legend.position="none")+
#   xlab("Encorafinib (BRAF inhibitor)") + ylab("Trametinib (MAP2K1 inhibitor)")
# new1<-ggMarginal(g1,type = "histogram", groupColour = TRUE, groupFill = TRUE)
# 
# g2 <- ggplot(df2, aes(Neratinib, Trametinib, colour = annot)) +
#   geom_point() + cowplot::theme_cowplot() + 
#   geom_hline(yintercept=0, linetype='dotted', col = 'darkred')+
#   geom_vline(xintercept=0, linetype='dotted', col = 'darkred')+
#   scale_color_ghibli_d("YesterdayMedium", direction = -1)+ 
#   theme(legend.position="none")+
#   xlab("Neratinib (EGFR inhibitor)") + ylab("Trametinib (MAP2K1 inhibitor)")
# new2<-ggMarginal(g2,type = "histogram", groupColour = TRUE, groupFill = TRUE)
# 
# g3 <- ggplot(df2, aes(AZ5363, Trametinib, colour = annot)) +
#   geom_point() + cowplot::theme_cowplot() + 
#   geom_hline(yintercept=0, linetype='dotted', col = 'darkred')+
#   geom_vline(xintercept=0, linetype='dotted', col = 'darkred')+
#   scale_color_ghibli_d("YesterdayMedium", direction = -1)+ 
#   theme(legend.position="none")+
#   xlab("AZ5363 (AKT inhibitor)") + ylab("Trametinib (MAP2K1 inhibitor)")
# new3<-ggMarginal(g3,type = "histogram", groupColour = TRUE, groupFill = TRUE)
# 
