
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


enrichr_bioplanet2019<-read.csv(file = "./results/phuego/enrichr_dfs/factor1.csv")
enrichr_bioplanet2019$label<-
enrichr_bioplanet2019$label<-""
enrichr_bioplanet2019$label[enrichr_bioplanet2019$log_pval>10]<-enrichr_bioplanet2019$term[enrichr_bioplanet2019$log_pval>10]
enrichr_bioplanet2019$label[enrichr_bioplanet2019$odds_ratio>50]<-enrichr_bioplanet2019$term[enrichr_bioplanet2019$odds_ratio>50]

pdf(file = paste0("/Users/charliebarker/Desktop/Melanoma_Resistance/paper/for_sumana/factor1/enrichr_enrichment.pdf"), 
    width = 8, height = 8)
# basic scatterplot
ggplot(enrichr_bioplanet2019, aes(x=odds_ratio, y=log_pval, label=label)) + 
  geom_point() +
  ggrepel::geom_text_repel(size = 2) + cowplot::theme_cowplot()
dev.off()

#check the immune in the combination. 