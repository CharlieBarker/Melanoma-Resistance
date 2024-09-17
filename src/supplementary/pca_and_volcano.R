home_dir<-"~/Desktop/Melanoma_Resistance/"
setwd(home_dir)

source("./src/functions/default_variables.R")

library(tidyr)
library(ggplot2)
library(limma)
library(dplyr)
library(ggpubr)
library(scatterplot3d)
library(ghibli)
library(ggrepel)

####Preprocessing####


list_of_inputs<-list(`mRNA RNAseq/transcriptomics`=data.frame(read.csv("./data/input_data/rna_expression.csv")),
                     `Protein abundance`=data.frame(read.csv("./data/input_data/proteins.csv")),
                     `Phosphoproteomic abundance`=data.frame(read.csv("./data/input_data/phosphosites.csv"))
)

#replace UNIPROT for genenames within protein abundances
protein<-list_of_inputs$`Protein abundance`$X
genename_df<-AnnotationDbi::select(EnsDb.Hsapiens.v86, keys = protein, keytype = "UNIPROTID", columns = "GENENAME")
# Collapse gene names
genename_df <- genename_df %>%
  group_by(UNIPROTID) %>%
  summarise(GENENAME = paste(unique(GENENAME), collapse = ", "))
genenames_available<-protein %in% genename_df$UNIPROTID
protein[genenames_available]<-genename_df$GENENAME[match(protein[genenames_available], genename_df$UNIPROTID)]
list_of_inputs$`Protein abundance`$X<-make.unique(protein)

process_columns<-function(counts, replacement_Vec){
  # Use X column to make row names
  rownames(counts) <- counts$X
  counts <- counts[, -1]  # Remove the 'X' column
  # Remove the number after the full stop in the column names
  colnames(counts) <- gsub("\\.\\d+", "", colnames(counts))
  # Extract drug names and genetic elements from column names
  col_names <- colnames(counts)
  drug_names <- sub("^(.*?)__.*", "\\1", col_names)
  gene_names <- sub(".*?__(.*)", "\\1", col_names)
  # Replace drug names with the names in replacement vector
  new_drug_names <- names(replacement_Vec)[match(drug_names, unname(replacement_Vec))]
  # Combine drug names and genetic elements to form new column names
  snames <- paste(new_drug_names, gene_names, sep = ".")
  # Assign new column names
  colnames(counts) <- snames
  out<-list()
  out$counts <- counts
  out$drug_names <- new_drug_names
  out$gene_names <- gene_names
  return(out)
}

####Produce PCA plots####

plot_two_pca<-function(counts,
                       title = NULL){
  to_pca <- counts
  out<-process_columns(to_pca, replacement_Vec)
  to_pca <- out$counts[!rowSums(is.infinite(as.matrix(out$counts))), ]

  pca <- prcomp(t(to_pca), scale.=TRUE )
  pca<-pca$x
  drug_names <- out$drug_names
  gene_names <- out$gene_names
  mds_df<-rbind(data.frame(`Drug treatment`=drug_names, `Gene knockout`=gene_names,
                     x=pca[,1], y=pca[,2], type="PC1 (x axis) vs PC2 (y axis)"),
                data.frame(`Drug treatment`=drug_names, `Gene knockout`=gene_names,
                           x=pca[,1], y=pca[,3], type="PC1 (x axis) vs PC3 (y axis)"))

  plot_mds<-ggplot(mds_df, aes(x=x, y=y, color=Drug.treatment, shape=Gene.knockout)) +
    geom_point(size=2) +
    cowplot::theme_cowplot() +
    scale_colour_manual(values=drug_colours) +
    grids(linetype = "dashed") + facet_wrap(~type)+
    ggtitle(title) +
    theme(plot.title = element_text(size = 12, face = "bold"))

  mds_df$title = title
  return(mds_df)
}


pca_list<-lapply(seq(length(names(list_of_inputs))),
                 function(x){plot_two_pca(list_of_inputs[[x]],
                                          title = paste0("Variation in ",
                                                         names(list_of_inputs)[x]))})

mds_big_df <- do.call("rbind", pca_list)
pdf(# The directory you want to save the file in
  width = 12, # The width of the plot in inches
  height = 12,
  file = "./paper/Supplementary_plots/raw_data.pdf")

ggplot(mds_big_df, aes(x=x, y=y, color=Drug.treatment, shape=Gene.knockout)) +
  geom_point(size=4) +
  cowplot::theme_cowplot() +
  scale_colour_manual(values=drug_colours) +
  grids(linetype = "dashed") +
  facet_grid(title ~ type) +
  theme(plot.title = element_text(size = 8, face = "bold"),
        panel.border = element_rect(colour = "black", fill = NA, linewidth = 1)) # Use linewidth instead of size

dev.off()


####Produce volcano plots####


produce_lfc <- function(counts, contrast_matrix, replacement_Vec) {
  out<-process_columns(counts, replacement_Vec)
  counts<-out$counts
  drug_names <- out$drug_names
  gene_names <- out$gene_names
  group <- interaction(drug_names, gene_names)

  counts <- counts[!rowSums(is.infinite(as.matrix(counts))), ]

  mm <- model.matrix(~0 + group)
  fit <- lmFit(counts, mm)

  contr <- contrast_matrix
  tmp <- contrasts.fit(fit, contr)
  tmp <- eBayes(tmp)
  top.table <- topTable(tmp, sort.by = "P", n = Inf)

  nsig<-length(which(top.table$adj.P.Val < 0.05))

  top.table$Gene <- rownames(top.table)
  top.table <- top.table[,c("Gene", names(top.table)[1:6])]
  return(list(top_table = top.table, nsig = nsig))
}



fit_levels<-c("groupCombinations.ARID1A_KO", "groupTrametinib.ARID1A_KO",   "groupUntreated.ARID1A_KO",
              "groupVemurafenib.ARID1A_KO",  "groupCombinations.WT",        "groupTrametinib.WT",
              "groupUntreated.WT",          "groupVemurafenib.WT")
contr_arid1a <- makeContrasts(groupUntreated.ARID1A_KO - groupUntreated.WT, levels = fit_levels)
contr_combination <- makeContrasts(groupCombinations.WT - groupUntreated.WT, levels = fit_levels)
contr_trametinib <- makeContrasts(groupTrametinib.WT - groupUntreated.WT, levels = fit_levels)
contr_vemurafenib <- makeContrasts(groupVemurafenib.WT - groupUntreated.WT, levels = fit_levels)

#under ARID1A KO
contr_combination_arid1a <- makeContrasts(groupCombinations.ARID1A_KO - groupUntreated.ARID1A_KO, levels = fit_levels)
contr_trametinib_arid1a <- makeContrasts(groupTrametinib.ARID1A_KO - groupUntreated.ARID1A_KO, levels = fit_levels)
contr_vemurafenib_arid1a <- makeContrasts(groupVemurafenib.ARID1A_KO - groupUntreated.ARID1A_KO, levels = fit_levels)


#for ARID1A
#arid1a vs wt
lfc_arid1a<-lapply(list_of_inputs, function(x){
  lfc<-produce_lfc(counts = x,
                        contrast_matrix = contr_arid1a,
                        replacement_Vec = replacement_Vec)
  return(lfc$top_table)
})
#arid1a combination
lfc_combination_arid1a<-lapply(list_of_inputs, function(x){
  lfc<-produce_lfc(counts = x,
                   contrast_matrix = contr_combination_arid1a,
                   replacement_Vec = replacement_Vec)
  return(lfc$top_table)
})
write.csv(x = lfc_combination_arid1a$`mRNA RNAseq/transcriptomics`,
          file = "./results/lfc/mRNA/combination_arid1a_lfc.csv")
write.csv(x = lfc_combination_arid1a$`Protein abundance`,
          file = "./results/lfc/protein/combination_arid1a_lfc.csv")
#arid1a trametinib
lfc_trametinib_arid1a<-lapply(list_of_inputs, function(x){
  lfc<-produce_lfc(counts = x,
                   contrast_matrix = contr_trametinib_arid1a,
                   replacement_Vec = replacement_Vec)
  return(lfc$top_table)
})
write.csv(x = lfc_trametinib_arid1a$`mRNA RNAseq/transcriptomics`,
          file = "./results/lfc/mRNA/trametinib_arid1a_lfc.csv")
write.csv(x = lfc_trametinib_arid1a$`Protein abundance`,
          file = "./results/lfc/protein/trametinib_arid1a_lfc.csv")
#arid1a vemurafenib
lfc_vemurafenib_arid1a<-lapply(list_of_inputs, function(x){
  lfc<-produce_lfc(counts = x,
                   contrast_matrix = contr_vemurafenib_arid1a,
                   replacement_Vec = replacement_Vec)
  return(lfc$top_table)
})
write.csv(x = lfc_vemurafenib_arid1a$`mRNA RNAseq/transcriptomics`,
          file = "./results/lfc/mRNA/vemurafenib_arid1a_lfc.csv")
write.csv(x = lfc_vemurafenib_arid1a$`Protein abundance`,
          file = "./results/lfc/protein/vemurafenib_arid1a_lfc.csv")

#print lfc for RNAseq ARID1A - so that we can generate figure 4
write.csv(x = lfc_arid1a$`mRNA RNAseq/transcriptomics`,
          file = "./results/lfc/mRNA/arid1a_lfc.csv")
write.csv(x = lfc_arid1a$`Protein abundance`,
          file = "./results/lfc/protein/arid1a_lfc.csv")
write.csv(x = lfc_arid1a$`Phosphoproteomic abundance`,
          file = "./results/lfc/phospho//arid1a_lfc.csv")

lfc_arid1a<-bind_rows(lfc_arid1a, .id = "column_label")


#for COMBINATION
lfc_combination<-lapply(list_of_inputs, function(x){
  lfc<-produce_lfc(counts = x,
                   contrast_matrix = contr_combination,
                   replacement_Vec = replacement_Vec)
  return(lfc$top_table)
})
write.csv(x = lfc_combination$`mRNA RNAseq/transcriptomics`,
          file = "./results/lfc/mRNA/combination_lfc.csv")
write.csv(x = lfc_combination$`Protein abundance`,
          file = "./results/lfc/protein/combination_lfc.csv")
write.csv(x = lfc_combination$`Phosphoproteomic abundance`,
          file = "./results/lfc/phospho/combination_lfc.csv")

lfc_combination<-bind_rows(lfc_combination, .id = "column_label")

#for TRAMETINIB
lfc_trametinib<-lapply(list_of_inputs, function(x){
  lfc<-produce_lfc(counts = x,
                   contrast_matrix = contr_trametinib,
                   replacement_Vec = replacement_Vec)
  return(lfc$top_table)
})
write.csv(x = lfc_trametinib$`mRNA RNAseq/transcriptomics`,
          file = "./results/lfc/mRNA/trametinib_lfc.csv")
write.csv(x = lfc_trametinib$`Protein abundance`,
          file = "./results/lfc/protein/trametinib_lfc.csv")
write.csv(x = lfc_trametinib$`Phosphoproteomic abundance`,
          file = "./results/lfc/phospho/trametinib_lfc.csv")

lfc_trametinib<-bind_rows(lfc_trametinib, .id = "column_label")

#for VEMURAFENIB
lfc_vemurafenib<-lapply(list_of_inputs, function(x){
  lfc<-produce_lfc(counts = x,
                   contrast_matrix = contr_vemurafenib,
                   replacement_Vec = replacement_Vec)
  return(lfc$top_table)
})
write.csv(x = lfc_vemurafenib$`mRNA RNAseq/transcriptomics`,
          file = "./results/lfc/mRNA/vemurafenib_lfc.csv")
write.csv(x = lfc_vemurafenib$`Protein abundance`,
          file = "./results/lfc/protein/vemurafenib_lfc.csv")
write.csv(x = lfc_vemurafenib$`Phosphoproteomic abundance`,
          file = "./results/lfc/phospho/vemurafenib_lfc.csv")

lfc_vemurafenib<-bind_rows(lfc_vemurafenib, .id = "column_label")


plot_volcano <- function(to_plot, title,
                         significant_parameters = list("a"=2,   #horizontal asymptote.
                                                       "b"=.5), #vertical asymptote.
                         labelling_parameters = NULL         #seperate parameters for labelling.
                         ) {

  # Define the parameters
  c <- 0   # Y-intercept

  # Internal function to define the mirrored function
  mirrored_asymptotic_function <- function(x,
                                           a=significant_parameters$a,
                                           b=significant_parameters$b
                                           ) {
    c=0
    y <- a / (abs(x) - b) + c
    return(y)
  }

  # Identify points above and below the mirrored function
  which_significant<- mirrored_asymptotic_function(to_plot$logFC)
  to_plot$below <- -log10(to_plot$adj.P.Val) < which_significant
  to_plot$below[abs(to_plot$logFC) < significant_parameters$b] <- TRUE

  # Define alpha values based on 'below'
  alpha_values <- ifelse(to_plot$below, 0.1, 0.5)

  to_plot$label <- ""

  #if labelling parameters are supplied, add labells
  if(!is.null(labelling_parameters)){
    which_labelled<- mirrored_asymptotic_function(to_plot$logFC, a=labelling_parameters$a, b=labelling_parameters$b)
    to_plot$to_label <- -log10(to_plot$adj.P.Val) < which_labelled
    to_plot$below[abs(to_plot$logFC) < labelling_parameters$b] <- TRUE
    to_plot$label[!to_plot$below] <- to_plot$Gene[!to_plot$below]
  }

  volcano_plot <- ggplot(to_plot,
                         aes(logFC, -log10(adj.P.Val), color = below, label = label)) +
    geom_point(alpha = alpha_values) +
    cowplot::theme_cowplot() +
    geom_text_repel(min.segment.length = 0, seed = 42, box.padding = 0.5,
                    color = "black",
                    size=2) +
    geom_vline(xintercept = 0, linetype = 'dotted', col = 'darkred') +
    theme(legend.position = "none") +
    ylab("Log10(Adjusted P value)") + xlab("Log fold change") +
    geom_function(fun = mirrored_asymptotic_function,
                  colour = ghibli_palettes$YesterdayDark[4], alpha = 0.5) +
    scale_color_manual(values = c("TRUE" = "lightgrey", "FALSE" = "#A42820")) +
    ylim(0, max(-log10(to_plot$adj.P.Val))) +
    ggtitle(title) +
    grids(linetype = "dashed") +
    theme(plot.title = element_text(size = 8, face = "bold"),
          panel.border = element_rect(colour = "black", fill = NA, linewidth = 1)) # Use linewidth instead of size


  return(volcano_plot)
}


pdf(
  width = 8.3,  # Width of the plot in inches
  height = 11.7, # Height of the plot in inches
  file = "./paper/Supplementary_plots/volcano_plots_by_data_type.pdf"
)

# Set the parameters for drug labelling
drug_labelling_params <- list("a" = 3, "b" = 0.8)
arid1a_curve_params <- list("a" = 1, "b" = 0.1)

# Page 1: mRNA RNAseq/transcriptomics plots
ggpubr::ggarrange(
  plot_volcano(lfc_trametinib[lfc_trametinib$column_label == "mRNA RNAseq/transcriptomics",],
               title = "Expressed genes, Untreated WT vs Trametinib-treated WT",
               labelling_parameters = drug_labelling_params),
  plot_volcano(lfc_vemurafenib[lfc_vemurafenib$column_label == "mRNA RNAseq/transcriptomics",],
               title = "Expressed genes, Untreated WT vs Vemurafenib-treated WT",
               labelling_parameters = drug_labelling_params),
  plot_volcano(lfc_combination[lfc_combination$column_label == "mRNA RNAseq/transcriptomics",],
               title = "Expressed genes, Untreated WT vs Combination-treated WT",
               labelling_parameters = drug_labelling_params),
  plot_volcano(lfc_arid1a[lfc_arid1a$column_label == "mRNA RNAseq/transcriptomics",],
               title = "Expressed genes, Untreated WT vs Untreated-ARID1A KO",
               significant_parameters = arid1a_curve_params, labelling_parameters = arid1a_curve_params),
  nrow = 4
)

# Page 2: Protein abundance plots
ggpubr::ggarrange(
  plot_volcano(lfc_trametinib[lfc_trametinib$column_label == "Protein abundance",],
               title = "Abundant proteins, Untreated WT vs Trametinib-treated WT",
               labelling_parameters = drug_labelling_params),
  plot_volcano(lfc_vemurafenib[lfc_vemurafenib$column_label == "Protein abundance",],
               title = "Abundant proteins, Untreated WT vs Vemurafenib-treated WT",
               labelling_parameters = drug_labelling_params),
  plot_volcano(lfc_combination[lfc_combination$column_label == "Protein abundance",],
               title = "Abundant proteins, Untreated WT vs Combination-treated WT",
               labelling_parameters = drug_labelling_params),
  plot_volcano(lfc_arid1a[lfc_arid1a$column_label == "Protein abundance",],
               title = "Abundant proteins, Untreated WT vs Untreated-ARID1A KO",
               significant_parameters = arid1a_curve_params, labelling_parameters = arid1a_curve_params),
  nrow = 4
)

# Page 3: Phosphoproteomic abundance plots
ggpubr::ggarrange(
  plot_volcano(lfc_trametinib[lfc_trametinib$column_label == "Phosphoproteomic abundance",],
               title = "Abundant phosphopeptides, Untreated WT vs Trametinib-treated WT",
               labelling_parameters = drug_labelling_params),
  plot_volcano(lfc_vemurafenib[lfc_vemurafenib$column_label == "Phosphoproteomic abundance",],
               title = "Abundant phosphopeptides, Untreated WT vs Vemurafenib-treated WT",
               labelling_parameters = drug_labelling_params),
  plot_volcano(lfc_combination[lfc_combination$column_label == "Phosphoproteomic abundance",],
               title = "Abundant phosphopeptides, Untreated WT vs Combination-treated WT",
               labelling_parameters = drug_labelling_params),
  plot_volcano(lfc_arid1a[lfc_arid1a$column_label == "Phosphoproteomic abundance",],
               title = "Abundant phosphopeptides, Untreated WT vs Untreated-ARID1A KO",
               significant_parameters = arid1a_curve_params, labelling_parameters = arid1a_curve_params),
  nrow = 4
)

# Close the PDF device
dev.off()


####Correlation plots####

protein<-process_columns(list_of_inputs$`Protein abundance`, replacement_Vec)
molten_protein<-reshape2::melt(data.matrix(protein$counts))
colnames(molten_protein)<-c("GENE", "CONDITION", "Protein Abundance")
molten_protein <- molten_protein %>%
  separate_rows(GENE, sep = ", ")
molten_protein <- molten_protein %>%
  group_by(GENE, CONDITION) %>%
  summarise(`Mean Protein Abundance` = mean(`Protein Abundance`))

mRNA<-process_columns(list_of_inputs$`mRNA RNAseq/transcriptomics`, replacement_Vec)
molten_mRNA<-reshape2::melt(data.matrix(mRNA$counts))
colnames(molten_mRNA)<-c("GENE", "CONDITION", "Voom-normalised mRNA abundance")
molten_mRNA <- molten_mRNA %>%
  group_by(GENE, CONDITION) %>%
  summarise(`Mean Voom-normalised mRNA abundance` = mean(`Voom-normalised mRNA abundance`))

phospho_processed<-process_columns(list_of_inputs$`Phosphoproteomic abundance`, replacement_Vec)
molten_phospho_processed<-reshape2::melt(data.matrix(phospho_processed$counts))
colnames(molten_phospho_processed)<-c("GENE", "CONDITION", "peptide abundance")
molten_phospho_processed <- molten_phospho_processed %>%
  group_by(GENE, CONDITION) %>%
  summarise(`Mean peptide abundance` = mean(`peptide abundance`))


#merge protein and mRNA readings
corr_plot_df<-merge(x = molten_mRNA, y = molten_protein,
                    by = c("GENE", "CONDITION"))
pdf(# The directory you want to save the file in
  width = 10, # The width of the plot in inches
  height = 11.7,
  file = "./paper/Supplementary_plots/correlation_plots.pdf")
ggplot(corr_plot_df, aes(x = `Mean Voom-normalised mRNA abundance`,
                         y = `Mean Protein Abundance`)) +
  geom_point(alpha=0.1) +
  geom_density_2d(bins=20) + facet_wrap(~CONDITION, ncol = 2) + cowplot::theme_cowplot() + stat_cor(method = "pearson", label.x = 3, label.y = 20)  +
  grids(linetype = "dashed") +
  theme(plot.title = element_text(size = 8, face = "bold"),
        panel.border = element_rect(colour = "black", fill = NA, linewidth = 1)) # Use linewidth instead of size


dev.off()

# Reorganize the dataframe

corr_plot_df_reorganized <- corr_plot_df %>%
  # Pivot longer to create a 'data_type' column for mRNA and protein
  pivot_longer(cols = c(`Mean Voom-normalised mRNA abundance`, `Mean Protein Abundance`),
               names_to = "data_type", values_to = "abundance") %>%
  # Rename the data_type values for clarity
  mutate(data_type = ifelse(data_type == "Mean Voom-normalised mRNA abundance", "mRNA abundance", "Protein abundance")) %>%
  # Pivot wider to create columns for each condition
  pivot_wider(names_from = CONDITION, values_from = abundance)

# Subset the dataframe for mRNA abundance and exclude ARID1A conditions
mRNA_df <- corr_plot_df_reorganized %>%
  dplyr::filter(data_type == "mRNA abundance") %>%
  dplyr::select(-data_type, -contains("ARID1A"))  # Remove the 'data_type' column and all columns containing 'ARID1A'

# Subset the dataframe for Protein abundance and exclude ARID1A conditions
protein_df <- corr_plot_df_reorganized %>%
  dplyr::filter(data_type == "Protein abundance") %>%
  dplyr::select(-data_type, -contains("ARID1A"))  # Remove the 'data_type' column and all columns containing 'ARID1A'

# Ensure no infinite or NA values in 'Mean peptide abundance'
molten_phospho_processed_clean <- molten_phospho_processed %>%
  dplyr::filter(is.finite(`Mean peptide abundance`)) %>%
  dplyr::filter(!is.na(`Mean peptide abundance`))

# Pivot the dataframe to a wide format
phospho_df <- molten_phospho_processed_clean %>%
  pivot_wider(names_from = CONDITION, values_from = `Mean peptide abundance`)

phospho_df <- phospho_df %>%
  dplyr::select(-contains("ARID1A"))

library(GGally)

pdf(
  width = 8,  # The width of the plot in inches
  height = 8, # The height of the plot in inches
  file = "./paper/Supplementary_plots/correlation_cond_vs_cond.pdf"
)

# Ensure no NA values
mRNA_df <- na.omit(mRNA_df)
protein_df <- na.omit(protein_df)
phospho_df <- na.omit(phospho_df)

# Ensure you have no infinite values in any column of protein_df
protein_df <- protein_df %>%
  dplyr::filter(across(everything(), is.finite))  # Keep rows where all values are finite

# Plot pairwise relationships for mRNA abundance excluding ARID1A conditions
ggpairs(
  mRNA_df,
  columns = 2:ncol(mRNA_df),
  title = "mRNA abundance plotted pairwise against all drug conditions",
  upper = list(continuous = wrap("cor", size = 3)),  # Add correlation coefficients
  diag = list(continuous = wrap("barDiag", fill = "lightblue", binwidth = 1)),  # Pale blue histograms with adjusted binwidth
  lower = list(continuous = wrap("points", alpha = 0.5))  # Slightly transparent points
) +
  cowplot::theme_cowplot() +
  theme(
    plot.title = element_text(size = 14, face = "bold"),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 1) # Use linewidth instead of size
  ) +
  grids(linetype = "dashed")

# Plot pairwise relationships for Protein abundance excluding ARID1A conditions
ggpairs(
  protein_df,
  columns = 2:ncol(protein_df),
  title = "Protein Abundance plotted pairwise against all drug conditions",
  upper = list(continuous = wrap("cor", size = 3)),  # Add correlation coefficients
  diag = list(continuous = wrap("barDiag", fill = "lightblue", binwidth = 1)),  # Pale blue histograms with adjusted binwidth
  lower = list(continuous = wrap("points", alpha = 0.5))  # Slightly transparent points
) +
  cowplot::theme_cowplot() +
  theme(
    plot.title = element_text(size = 14, face = "bold"),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 1) # Use linewidth instead of size
  ) +
  grids(linetype = "dashed")


# # Plot pairwise relationships for peptide abundance excluding ARID1A conditions
# ggpairs(
#   phospho_df,
#   columns = 2:ncol(phospho_df),
#   title = "Pairwise Plot: Peptide Abundance Conditions",
#   upper = list(continuous = wrap("cor", size = 3)),  # Add correlation coefficients
#   diag = list(continuous = wrap("barDiag", fill = "lightblue", binwidth = 1)),  # Pale blue histograms with adjusted binwidth
#   lower = list(continuous = wrap("points", alpha = 0.5))  # Slightly transparent points
# ) +
#   cowplot::theme_cowplot() +
#   theme(
#     plot.title = element_text(size = 14, face = "bold"),
#     panel.border = element_rect(colour = "black", fill = NA, linewidth = 1) # Use linewidth instead of size
#   ) +
#   grids(linetype = "dashed")


dev.off()


####Upset plots####
library(ComplexHeatmap)

phospho=data.frame(read.csv("./data/proteomic/processed/combat_peptide.csv"))
phosphorylated_uniprot_quantified <- gsub("^.*?__(.*?)\\s.*$", "\\1", phospho$X)
genename_df<-AnnotationDbi::select(EnsDb.Hsapiens.v86, keys = phosphorylated_uniprot_quantified, keytype = "UNIPROTID", columns = "GENENAME")
phosphorylated_proteins_quantified <- unique(genename_df$GENENAME)

lt<-list(`mRNA RNAseq/transcriptomics`= unique(list_of_inputs$`mRNA RNAseq/transcriptomics`$X),
         `Protein abundance`= unique(list_of_inputs$`Protein abundance`$X),
         `Phosphoproteomic abundance`=phosphorylated_proteins_quantified)
pdf(# The directory you want to save the file in
  width = 8.3, # The width of the plot in inches
  height = 11.7,
  file = "./paper/Supplementary_plots/upset_plot.pdf")
m1 = make_comb_mat(lt, mode = "distinct")
m2 = make_comb_mat(lt, mode = "intersect")
m3 = make_comb_mat(lt, mode = "union")
top_ha = HeatmapAnnotation(
  "Distict proteins/genes" = anno_barplot(comb_size(m1),
                           gp = gpar(fill = "black"), height = unit(4, "cm")),
  "Intersect of proteins/genes" = anno_barplot(comb_size(m2),
                             gp = gpar(fill = "black"), height = unit(4, "cm")),
  "Union of proteins/genes" = anno_barplot(comb_size(m3),
                         gp = gpar(fill = "black"), height = unit(4, "cm")),
  gap = unit(2, "mm"), annotation_name_side = "left", annotation_name_rot = 0)
# the same for using m2 or m3
UpSet(m1, top_annotation = top_ha)
dev.off()


####MOFA plot####

MOFAobject.trained<-load_model(file = "./results/mofa/mofa_object.hdf5")

weights <- get_weights(MOFAobject.trained,
                       views = "all",
                       as.data.frame = TRUE
)

# Load the dplyr package
library(dplyr)

# Assuming 'weights' is your dataframe

# Define a function to add label based on percentiles
add_label <- function(df, grouping_vars) {
  df %>%
    group_by(across(all_of(grouping_vars))) %>%
    mutate(label = case_when(
      value > quantile(value, 0.99) ~ feature,
      value < quantile(value, 0.01) ~ feature,
      TRUE ~ ""
    )) %>%
    ungroup()
}

# Add labels based on percentiles within factors
df <- weights %>%
  group_by(factor) %>%
  mutate(rank_within_factor = rank(value)) %>%
  add_label(c("factor"))

# Add labels based on percentiles within factors and views
df <- df %>%
  group_by(factor, view) %>%
  mutate(rank_within_view_and_factor = rank(value)) %>%
  add_label(c("factor", "view"))

ghibli_cs<-"YesterdayMedium"


pdf(# The directory you want to save the file in
  width = 8.3, # The width of the plot in inches
  height = 11.7,
  file = "./paper/Supplementary_plots/mofa_plot.pdf")
ggplot(df, aes(x = rank_within_factor,
               y = value,
               color= view,
               label=label)) +
  geom_point(alpha=.6) +
  facet_wrap(~factor, ncol = 3, scales = "free_y") + cowplot::theme_cowplot() +
  ggtitle("MOFA weights for the different factors") +
  theme(plot.title = element_text(size = 12, face = "bold"),
        panel.border = element_rect(colour = "black", fill=NA, size=1))+
  scale_color_ghibli_d(ghibli_cs, direction = -1) +
  geom_text_repel(min.segment.length = 0, seed = 42, box.padding = 0.5,
                  color = "black",
                  size=2, force_pull = 2)+
  geom_hline(yintercept = 0, linetype = 'dotted', col = 'darkred')
ggplot(df, aes(x = rank_within_factor,
               y = value,
               color= view,
               label=label)) +
  geom_point(alpha=.6) +
  facet_grid(view~factor, scales = "free_y") + cowplot::theme_cowplot() +
  ggtitle("MOFA weights for the different factors") +
  theme(plot.title = element_text(size = 12, face = "bold"),
        panel.border = element_rect(colour = "black", fill=NA, size=1))+
  scale_color_ghibli_d(ghibli_cs, direction = -1) +
  geom_text_repel(min.segment.length = 0, seed = 42, box.padding = 0.5,
                  color = "black",
                  size=2, force_pull = 2)+
  geom_hline(yintercept = 0, linetype = 'dotted', col = 'darkred')
dev.off()


