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

list_of_inputs<-list(#signalome_view=data.frame(read.csv("input_data/signalome_by_sample.csv")),
  phospho=data.frame(read.csv("./data/input_data/phosphosites.csv"))
  ,protein=data.frame(read.csv("./data/input_data/proteins.csv"))
  ,mRNA=data.frame(read.csv("./data/input_data/rna_expression.csv"))
)




plot_two_pca<-function(counts){
  
  #counts<-read.csv(file = "./data/input_data/proteins.csv")#log2(mofa_protein)
  # Use X column to make row names
  rownames(counts) <- counts$X
  counts <- counts[, -1]  # Remove the 'X' column
  # Remove the number after the full stop in the column names
  colnames(counts) <- gsub("\\.\\d+", "", colnames(counts))
  
  to_pca <- counts[!rowSums(is.infinite(as.matrix(counts))), ]
  dim(to_pca)
  pca <- prcomp(t(to_pca), scale.=TRUE )
  pca<-pca$x
  # Remove the number after the full stop in the column names
  rownames(pca) <- gsub("\\.\\d+", "", rownames(pca))
  # Extract drug names and genetic elements from column names
  col_names <- rownames(pca)
  drug_names <- sub("^(.*?)__.*", "\\1", col_names)
  gene_names <- sub(".*?__(.*)", "\\1", col_names)
  
  mds_df<-data.frame(drug=drug_names, gene=gene_names,
                     PC1=pca[,1], PC2=pca[,2], PC3=pca[,3])
  
  plot_mds<-ggpairs(mds_df,                 # Data frame
                    columns = 3:5,
                    lower = "blank", upper = list(continuous="points"),    # Columns
                    aes(color = gene_names, # Color by group (cat. variable)
                        shape = drug_names)) + cowplot::theme_cowplot() +
    scale_colour_manual(values=gene_colours) +
    scale_fill_manual(values=gene_colours) + 
    grids(linetype = "dashed")

  return(plot_mds)
}


pca_list<-lapply(list_of_inputs, plot_two_pca)
pdf(# The directory you want to save the file in
  width = 6, # The width of the plot in inches
  height = 6,
  file = "./paper/Supplementary_plots/raw_data.pdf")
pca_list
#genetic legend
plot(NULL ,xaxt='n',yaxt='n',bty='n',ylab='',xlab='', xlim=0:1, ylim=0:1)
legend("topleft", legend =names(gene_colours), pch=16, pt.cex=2, cex=1, bty='n',
       col = unname(gene_colours))
mtext("Genetic Condition", at=.1, cex=1)
#drug legend

plot(NULL ,xaxt='n',yaxt='n',bty='n',ylab='',xlab='', xlim=0:1, ylim=0:1)
legend("topleft", legend =names(drug_colours), pt.cex=2, cex=1, bty='n',
       pch = c(15,3,16,17))
mtext("Genetic Condition", at=.1, cex=1)
dev.off()

produce_lfc <- function(counts, contrast_matrix, replacement_Vec) {
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
  
  group <- interaction(new_drug_names, gene_names)
  
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

#for ARID1A
lfc_arid1a<-lapply(list_of_inputs, function(x){
  lfc<-produce_lfc(counts = x,
                        contrast_matrix = contr_arid1a,
                        replacement_Vec = replacement_Vec)
  return(lfc$top_table)
})
lfc_arid1a<-bind_rows(lfc_arid1a, .id = "column_label")

#for COMBINATION
lfc_combination<-lapply(list_of_inputs, function(x){
  lfc<-produce_lfc(counts = x,
                   contrast_matrix = contr_combination,
                   replacement_Vec = replacement_Vec)
  return(lfc$top_table)
})
lfc_combination<-bind_rows(lfc_combination, .id = "column_label")

#for TRAMETINIB
lfc_trametinib<-lapply(list_of_inputs, function(x){
  lfc<-produce_lfc(counts = x,
                   contrast_matrix = contr_trametinib,
                   replacement_Vec = replacement_Vec)
  return(lfc$top_table)
})
lfc_trametinib<-bind_rows(lfc_trametinib, .id = "column_label")

#for VEMURAFENIB
lfc_vemurafenib<-lapply(list_of_inputs, function(x){
  lfc<-produce_lfc(counts = x,
                   contrast_matrix = contr_vemurafenib,
                   replacement_Vec = replacement_Vec)
  return(lfc$top_table)
})
lfc_vemurafenib<-bind_rows(lfc_vemurafenib, .id = "column_label")


plot_volcano <- function(to_plot, title) {
  
  # Define the parameters
  a <- 2   # Horizontal asymptote
  b <- .5  # Vertical asymptote
  c <- 0   # Y-intercept
  
  # Internal function to define the mirrored function
  mirrored_asymptotic_function <- function(x) {
    y <- a / (abs(x) - b) + c
    return(y)
  }
  
  # Identify points above and below the mirrored function
  to_plot$below <- -log10(to_plot$adj.P.Val) < mirrored_asymptotic_function(to_plot$logFC)
  to_plot$below[abs(to_plot$logFC) < b] <- TRUE
  
  # Define alpha values based on 'below'
  alpha_values <- ifelse(to_plot$below, 0.1, 0.6)
  
  to_plot$label <- ""
  to_plot$label[!to_plot$below] <- to_plot$Gene[!to_plot$below]
  
  volcano_plot <- ggplot(to_plot, 
                         aes(logFC, -log10(adj.P.Val), color = below)) +
    geom_point(alpha = alpha_values) + 
    cowplot::theme_cowplot() + 
    geom_vline(xintercept = 0, linetype = 'dotted', col = 'darkred') +
    theme(legend.position = "none") +
    ylab("Log10(Adjusted P value)") + xlab("Log fold change") + 
    geom_function(fun = mirrored_asymptotic_function, colour = ghibli_palettes$YesterdayDark[4]) +
    scale_color_manual(values = c("TRUE" = "#FFECCC", "FALSE" = "#A4303F")) + 
    ylim(0, max(-log10(to_plot$adj.P.Val))) + 
    facet_wrap(~column_label, ncol = 3) +
    ggtitle(title)
  
  return(volcano_plot)
}


pdf(# The directory you want to save the file in
  width = 8.3, # The width of the plot in inches
  height = 11.7,
  file = "./paper/Supplementary_plots/volcano_plots.pdf")
ggpubr::ggarrange(plot_volcano(lfc_trametinib, title = "A.   Volcano plot of Untreated WT vs Trametinib-treated WT"),
                  plot_volcano(lfc_vemurafenib, title = "B.   Volcano plot of Untreated WT vs Vemurafenib-treated WT"),
                  plot_volcano(lfc_combination, title = "C.   Volcano plot of Untreated WT vs Combination-treated WT"),
                  plot_volcano(lfc_arid1a, title = "D.   Volcano plot of Untreated WT vs untreated ARID1A KO")
                  , nrow = 4) 

dev.off()


