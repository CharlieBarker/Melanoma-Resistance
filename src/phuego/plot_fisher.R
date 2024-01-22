
#find correct environment 
packLib="/usr/lib/R"
if (file.exists(packLib)) {
  reticulate::use_condaenv("py37", required = T)
  setwd(dir = "~/Desktop/Melanoma_Resistance/")
}else {
  setwd(dir = "~/Desktop/Melanoma_Resistance/")
}

library(dplyr)
library(tidyr)
library(viridis)
library(ggplot2)
library(wesanderson)
library(ggpubr)

read.fisher<-function(file){
  
  # Read the lines from the .txt file
  lines <- readLines(file)
  
  # Initialize variables to store data
  data <- list(
    DB_id = character(),
    adjusted_p = numeric(),
    nodes_in_the_network = integer(),
    N_of_seed_nodes = integer(),
    Description = character(),
    gene_names = list()
  )
  
  # Process the header separately
  header <- unlist(strsplit(lines[1], "\t"))
  for (i in 1:(length(header) - 1)) {
    data[[i]] <- header[i]
  }
  
  # Process the rest of the lines
  for (line in lines[-1]) {
    parts <- unlist(strsplit(line, "\t"))
    data$DB_id <- c(data$DB_id, parts[1])
    data$adjusted_p <- c(data$adjusted_p, as.numeric(parts[2]))
    data$nodes_in_the_network <- c(data$nodes_in_the_network, as.integer(parts[3]))
    data$N_of_seed_nodes <- c(data$N_of_seed_nodes, as.integer(parts[4]))
    data$Description <- c(data$Description, parts[5])
    data$gene_names <- append(data$gene_names, paste0(unlist(strsplit(parts[6:length(parts)], "\\s+")), collapse = ","))
  }
  
  # Create a data frame from the extracted data
  df <- data.frame(data)
  df <- data.frame(DB_id=df$DB_id[-1],
                   adjusted_p=as.numeric(df$adjusted_p[-1]),
                   nodes_in_the_network=df$nodes_in_the_network[-1],
                   N_of_seed_nodes=df$N_of_seed_nodes[-1], 
                   Description=df$Description[-1])
  # Display the resulting data frame
  return(df)
  
}


fisher_dir<-"./results/phuego/fisher_summary/"
fisher_summary<-list.files(path = fisher_dir, recursive = T)
all_fisher<-lapply(fisher_summary, function(file_name){
  read.fisher(paste0(fisher_dir, file_name))
})
names(all_fisher)<-fisher_summary

fisher_df<-bind_rows(all_fisher, .id = "column_label")
fisher_df <- fisher_df |>
  separate_wider_delim(column_label, delim = "/", names = c("Factor", "Sign"))
fisher_df$log10p = -log10(fisher_df$adjusted_p)

generate_factor_plot <- function(factor_df) {
  factor_df %>%
    ggplot(aes(x = Sign, y = reorder(Description, -log10p), color = log10p, size = -log10(adjusted_p))) +
    geom_point() +
    scale_color_viridis_c(name = 'Combined.Score') +
    cowplot::theme_cowplot() +
    theme(
      axis.line = element_blank(),
      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = rel(1)),
      axis.text.y = element_text(size = rel(.9)),
      axis.ticks = element_blank(),
      panel.border = element_rect(colour = "black", fill = NA, size = 2)
    ) +
    ylab('') +
    scale_shape_manual(values = c(1, 16)) +
    scale_fill_continuous(guide = guide_legend())+ 
    grids(linetype = "dashed") +
    scale_color_gradientn(colours = pal)+
    labs(
      x = "Sign",
      y = "Description",
      color = "-log10(FDR - P value)",
      size = "-log10(FDR - P value)"
    )
}

pal <- wes_palette("Zissou1", 100, type = "continuous")

Factor1 <- fisher_df %>%
  filter(Factor == "Factor1") %>%
  arrange(adjusted_p) %>%
  slice_head(n = 20) %>%
  generate_factor_plot()

Factor2 <- fisher_df %>%
  filter(Factor == "Factor2") %>%
  arrange(adjusted_p) %>%
  slice_head(n = 20) %>%
  generate_factor_plot()

Factor3 <- fisher_df %>%
  filter(Factor == "Factor3") %>%
  arrange(adjusted_p) %>%
  slice_head(n = 20) %>%
  generate_factor_plot()


output_file<-"./paper/plots/phuego/rwr_enrichment.pdf"
pdf(# The directory you want to save the file in
  width = 8, # The width of the plot in inches
  height = 9,
  file = output_file)
Factor1
Factor2
Factor3
dev.off()