
#find correct environment
packLib="/usr/lib/R"
if (file.exists(packLib)) {
  reticulate::use_condaenv("py37", required = T)
  setwd(dir = "~/Desktop/Melanoma_Resistance/")
}else {
  setwd(dir = "~/Desktop/Melanoma_Resistance/")
}


library(ggplot2)
library(ggpubr)
library(dplyr)
library(ggrepel)

log_pval_cutoff<-20
odds_ratio_cutoff<-60

enrichr_bioplanet2019<-read.csv(file = "./results/phuego/enrichr_dfs/factor3and1.csv")
enrichr_bioplanet2019$label<-
enrichr_bioplanet2019$label<-""
enrichr_bioplanet2019$label[enrichr_bioplanet2019$log_pval>log_pval_cutoff]<-enrichr_bioplanet2019$term[enrichr_bioplanet2019$log_pval>log_pval_cutoff]
enrichr_bioplanet2019$label[enrichr_bioplanet2019$odds_ratio>odds_ratio_cutoff]<-enrichr_bioplanet2019$term[enrichr_bioplanet2019$odds_ratio>odds_ratio_cutoff]

pdf(file = paste0("./paper/Supplementary_plots/enrichr_enrichment_factor3and1.pdf"),
    width = 8, height = 8)
# basic scatterplot
ggplot(enrichr_bioplanet2019, aes(x=odds_ratio, y=log_pval, label=label)) +
  geom_point() +
  ggrepel::geom_text_repel(size = 4) +
  cowplot::theme_cowplot() +
  theme(
    plot.title = element_text(size = 15, face = "bold"),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 1),
  ) +
  grids(linetype = "dashed") +
  xlab(paste("Odd Ratio")) +
  ylab("Log10(adjusted P value)") +
  ggtitle("Enrichment of terms in combined ARID1A/drug treatment network")
dev.off()

#check the immune in the combination.


log_pval_cutoff<-10
odds_ratio_cutoff<-30

list_factor2<-list(UpRegulated=read.csv(file = "./results/phuego/enrichr_dfs/factor2_up.csv"),
                   DownRegulated=read.csv(file = "./results/phuego/enrichr_dfs/factor2_down.csv"))
enrichr_bioplanet2019 <- bind_rows(list_factor2, .id = "direction")

enrichr_bioplanet2019$label<-
  enrichr_bioplanet2019$label<-""
enrichr_bioplanet2019$label[enrichr_bioplanet2019$log_pval>log_pval_cutoff]<-enrichr_bioplanet2019$term[enrichr_bioplanet2019$log_pval>log_pval_cutoff]
enrichr_bioplanet2019$label[enrichr_bioplanet2019$odds_ratio>odds_ratio_cutoff]<-enrichr_bioplanet2019$term[enrichr_bioplanet2019$odds_ratio>odds_ratio_cutoff]

pdf(file = paste0("./paper/Supplementary_plots/enrichr_enrichment_factor2.pdf"),
    width = 16, height = 6)
# basic scatterplot
ggplot(enrichr_bioplanet2019, aes(x=odds_ratio, y=log_pval, label=label)) +
  geom_point() +
  ggrepel::geom_text_repel(size = 4) +
  cowplot::theme_cowplot() +
  facet_wrap(~direction, ncol = 2) +
  theme(
    plot.title = element_text(size = 15, face = "bold"),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 1),
  ) +
  grids(linetype = "dashed") +
  xlab(paste("Odd Ratio")) +
  ylab("Log10(adjusted P value)") +
  ggtitle("Enrichment of terms in combined ARID1A/drug treatment network")
dev.off()

#check the immune in the combination.


log_pval_cutoff<-10
odds_ratio_cutoff<-35

enrichr_bioplanet2019<-read.csv(file = "./results/phuego/enrichr_dfs/factor2_down.csv")
enrichr_bioplanet2019$label<-""
enrichr_bioplanet2019$label[enrichr_bioplanet2019$log_pval>log_pval_cutoff]<-enrichr_bioplanet2019$term[enrichr_bioplanet2019$log_pval>log_pval_cutoff]
enrichr_bioplanet2019$label[enrichr_bioplanet2019$odds_ratio>odds_ratio_cutoff]<-enrichr_bioplanet2019$term[enrichr_bioplanet2019$odds_ratio>odds_ratio_cutoff]
list_to_keep <- c('ERBB signaling pathway', 'mTOR signaling pathway', 'IGF1 pathway',
                  'CBL−mediated ligand−induced downregulation of EGF receptors', 'Mammalian target of rapamycin complex 1 (mTORC1)−mediated signaling', 
                  'S6K1 signaling',
                  'Signaling events mediated by hepatocyte growth factor receptor (c−Met)', 
                  'PI3K events in ERBB2 signaling', 'Signaling by ERBB2')
#enrichr_bioplanet2019$label[!enrichr_bioplanet2019$label %in% list_to_keep]<- ''

pdf(file = paste0("./paper/Supplementary_plots/enrichr_enrichment_factor2.pdf"),
    width = 8, height = 8)
# Basic scatterplot with improved ggrepel
ggplot(enrichr_bioplanet2019, aes(x = odds_ratio, y = log_pval, label = label)) +
  geom_point() +
  ggrepel::geom_text_repel(
    min.segment.length = 0,
    force_pull = 0,
    size = 4, #should be 5
    segment.color = "black",  # Ensures lines are drawn to points
    segment.size = 0.5        # Adjust line thickness
  ) +
  cowplot::theme_cowplot(font_size = 20) +
  theme(
    plot.title = element_text(size = 15, face = "bold"),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 1)
  ) +
  grids(linetype = "dashed") +
  xlab("Odd Ratio") +
  ylab("Log10(adjusted P value)") +
  ggtitle("Enrichment of terms in downregulated combination-therapy network")

enrichr_bioplanet2019<-read.csv(file = "./results/phuego/enrichr_dfs/factor2_up.csv")
enrichr_bioplanet2019$label[enrichr_bioplanet2019$log_pval>log_pval_cutoff]<-enrichr_bioplanet2019$term[enrichr_bioplanet2019$log_pval>log_pval_cutoff]
enrichr_bioplanet2019$label[enrichr_bioplanet2019$odds_ratio>odds_ratio_cutoff]<-enrichr_bioplanet2019$term[enrichr_bioplanet2019$odds_ratio>odds_ratio_cutoff]

# Basic scatterplot with improved ggrepel
ggplot(enrichr_bioplanet2019, aes(x = odds_ratio, y = log_pval, label = label)) +
  geom_point() +
  ggrepel::geom_text_repel(
    min.segment.length = 0,
    size = 2, #should be 5
    segment.color = "black",  # Ensures lines are drawn to points
    segment.size = 0.5        # Adjust line thickness
  ) +
  cowplot::theme_cowplot() +
  theme(
    plot.title = element_text(size = 18, face = "bold"),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 1)
  ) +
  grids(linetype = "dashed") +
  xlab("Odd Ratio") +
  ylab("Log10(adjusted P value)") +
  ggtitle("Enrichment of terms in upregulated combination-therapy network")

dev.off()

pdf(file = paste0("./paper/Supplementary_plots/enrichr_enrichment_factor1.pdf"),
    width = 6)

enrichr_bioplanet2019<-read.csv(file = "./results/phuego/enrichr_dfs/factor1.csv")
enrichr_bioplanet2019$label[enrichr_bioplanet2019$log_pval>log_pval_cutoff]<-enrichr_bioplanet2019$term[enrichr_bioplanet2019$log_pval>log_pval_cutoff]
enrichr_bioplanet2019$label[enrichr_bioplanet2019$odds_ratio>odds_ratio_cutoff]<-enrichr_bioplanet2019$term[enrichr_bioplanet2019$odds_ratio>odds_ratio_cutoff]
list_to_keep <- c('Pathways in cancer', 'MAPK signaling pathway', 'Immune system',
                  'FGFR1b ligand binding and activation', 'MAP kinase pathway regulation through dual specificity phosphatases', 'ERK activation',
                  'Melanoma', 'Prolactin activation of MAPK signaling', 'FGF signaling pathway')
enrichr_bioplanet2019$label[!enrichr_bioplanet2019$label %in% list_to_keep]<- ''
# Basic scatterplot with improved ggrepel
ggplot(enrichr_bioplanet2019, aes(x = odds_ratio, y = log_pval, label = label)) +
  geom_point() +
  ggrepel::geom_text_repel(
    min.segment.length = 0,
    force_pull = 0,
    size = 5,
    segment.color = "black",  # Ensures lines are drawn to points
    segment.size = 0.5        # Adjust line thickness
  ) +
  cowplot::theme_cowplot(font_size = 20) +
  theme(
    plot.title = element_text(size = 15, face = "bold"),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 1)
  ) +
  grids(linetype = "dashed") +
  xlab("Odd Ratio") +
  ylab("Log10(adjusted P value)") +
  ggtitle("Enrichment of terms in combined drug-agnostic network")

dev.off()

