
#find correct environment
packLib="/usr/lib/R"
if (file.exists(packLib)) {
  reticulate::use_condaenv("py37", required = T)
  setwd(dir = "~/Desktop/Melanoma_Resistance/")
}else {
  setwd(dir = "~/Desktop/Melanoma_Resistance/")
}


library(ggplot2)
library(ggrepel)

log_pval_cutoff<-20
odds_ratio_cutoff<-60

enrichr_bioplanet2019<-read.csv(file = "./results/phuego/enrichr_dfs/factor3and1.csv")
enrichr_bioplanet2019$label<-
enrichr_bioplanet2019$label<-""
enrichr_bioplanet2019$label[enrichr_bioplanet2019$log_pval>log_pval_cutoff]<-enrichr_bioplanet2019$term[enrichr_bioplanet2019$log_pval>log_pval_cutoff]
enrichr_bioplanet2019$label[enrichr_bioplanet2019$odds_ratio>odds_ratio_cutoff]<-enrichr_bioplanet2019$term[enrichr_bioplanet2019$odds_ratio>odds_ratio_cutoff]

pdf(file = paste0("/Users/charliebarker/Desktop/Melanoma_Resistance/paper/Supplementary_plots/enrichr_enrichment_factor3and1.pdf"),
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

pdf(file = paste0("/Users/charliebarker/Desktop/Melanoma_Resistance/paper/Supplementary_plots/enrichr_enrichment_factor2.pdf"),
    width = 10, height = 16)
# basic scatterplot
ggplot(enrichr_bioplanet2019, aes(x=odds_ratio, y=log_pval, label=label)) +
  geom_point() +
  ggrepel::geom_text_repel(size = 4) +
  cowplot::theme_cowplot() +
  facet_wrap(~direction, ncol = 1) +
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

