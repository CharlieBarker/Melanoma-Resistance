library(reshape2)
library(GGally)
library(dplyr)
library(ggrepel)
library(tidyr)
library(tidyverse)
library(purrr)
setwd("~/Desktop/Melanoma_Resistance/data/Zhao_etal_2020/")

a375r<-read.csv(file = "lfc.csv")
arid1a<-read.csv(file = "../../results/transcriptomics/arid1a_lfc.csv")


all_df<-rbind(data.frame(geneNames=a375r$genes,
                         lfc=a375r$logFC,
                         padj=a375r$adj.P.Val,
                         exp="Vemurafenib resistant A375"),
              data.frame(geneNames=arid1a$gene_symbol,
                 lfc=arid1a$log2FoldChange,
                 padj=arid1a$padj,
                 exp="ARID1A KO A375"))
all_df<-all_df[all_df$padj<0.1,]

# Add a unique identifier column
# Group by the shared columns (geneNames, padj, exp)
all_df <- all_df %>%
  group_by(geneNames, exp) %>%
  summarise(lfc = mean(lfc),
            padj = mean(padj)) %>%
  ungroup()

spread_df <- pivot_wider(all_df[!is.na(all_df$geneNames),],
                         names_from = exp, values_from = lfc, id_cols = geneNames)
load("../../results/maxflow/max_flow_graphs.Rdata")

spread_df$in_factor1 <- spread_df$geneNames %in% c(igraph::V(max_flow_graphs$Factor1$up)$name,
                                                   igraph::V(max_flow_graphs$Factor1$down)$name)
spread_df$in_factor2 <- spread_df$geneNames %in% c(igraph::V(max_flow_graphs$Factor2$up)$name,
                                                   igraph::V(max_flow_graphs$Factor2$down)$name)
spread_df$in_factor3 <- spread_df$geneNames %in% c(igraph::V(max_flow_graphs$Factor3$up)$name,
                                                   igraph::V(max_flow_graphs$Factor3$down)$name)


# Change the point size, and shape
ggplot(spread_df, aes(x=`ARID1A KO A375`, y=`Vemurafenib resistant A375`, colour=in_factor1)) +
  geom_point() + 
  geom_smooth(method='lm')

to_plot<-all_df[all_df$geneNames %in% intersect(igraph::V(max_flow_graphs$Factor3$up)$name, c(receptors, ligands)),]
# Barplot
ggplot(to_plot, aes(x=exp, y=lfc)) + 
  geom_bar(stat = "identity") + facet_wrap(~geneNames, scales = "free")