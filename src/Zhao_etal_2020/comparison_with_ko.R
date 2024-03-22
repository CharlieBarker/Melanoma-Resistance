library(fgsea)
library(reshape2)
library(GGally)
library(dplyr)
library(ggrepel)

setwd("~/MelanomaProject/A375r/data")

a375r<-read.csv(file = "lfc.csv")
arid1a<-read.csv(file = "~/MelanomaProject/Sumana_RNAseq_analysis/WT_Vemurafenib_vs_ARID1A_Vemurafenib/WT_Vemurafenib_vs_ARID1A_Vemurafenib_diff_seq.csv")
med12<-read.csv(file = "~/MelanomaProject/Sumana_RNAseq_analysis/WT_Vemurafenib_vs_MED12_Vemurafenib/WT_Vemurafenib_vs_MED12_Vemurafenib_diff_seq.csv")
list_of_dsets<-list(ARID1A_KO=a375r,MED12_KO=med12,A375r=a375r)

shared_ensg<-Reduce(intersect,  list(v1 = a375r$X, v2 = arid1a$X, v3 = med12$X)) 
all_df<-rbind(data.frame(ensg=a375r$X[a375r$X %in% shared_ensg],
                         lfc=a375r$logFC[a375r$X %in% shared_ensg],
                         exp="A375r"),
      data.frame(ensg=arid1a$X[arid1a$X %in% shared_ensg],
                 lfc=arid1a$logFC[arid1a$X %in% shared_ensg],
                 exp="ARID1A_KO"),
      data.frame(ensg=med12$X[med12$X %in% shared_ensg],
                 lfc=med12$logFC[med12$X %in% shared_ensg],
                 exp="MED12_KO"))

spread_all<-dcast(all_df, ensg ~ exp, value.var = "lfc")
plot(spread_all[,2:4])

ENSEMBL_convert<- read.table("~/MelanomaProject/Sumana_RNAseq_analysis/ENSEMBL_mapping_human.txt", header=T)
spread_all$name<-ENSEMBL_convert$symbol[match(spread_all$ensg, ENSEMBL_convert$ENSEMBL_ID)]
#visualise the correlation with high weight mofa mrna
mrna_weights<-mofa_weights[mofa_weights$view == "mRNA",]


pthways <- gmtPathways("~/MelanomaProject/mofa/pathway_sets/c2.cp.v2022.1.Hs.symbols.gmt") #Human_GO_AllPathways_with_GO_iea_April_01_2020_symbol.gmt
pathOut<-list()
for (exper in c("ARID1A_KO", "MED12_KO", "A375r")) {
  df_test<-list_of_dsets[[exper]]
  ranks<-df_test$logFC
  names(ranks)<-df_test$genes
  # Load the pathways into a named list
  fgseaRes <- fgsea(pathways=pthways, stats=ranks, nperm=10000)
  pathOut[[exper]]<-fgseaRes
  cat(sum(fgseaRes$padj<0.1), "\n")
}



pathOut_df<-bind_rows(pathOut, .id = "EXP")

ggplot(pathOut_df[pathOut_df$padj<0.1,], aes(x=EXP, y = reorder(pathway, -NES), color = NES, size = -log10(padj))) + 
  geom_point() + 
  scale_color_viridis_c(name = 'Normalised Enrichment Score (NES)') + 
  cowplot::theme_cowplot() + 
  theme(axis.line  = element_blank()) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  ylab('') +
  theme(axis.ticks = element_blank())  +
  theme(axis.ticks = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=2)) 
up_in_all<-spread_all$name[spread_all$A375r > 1 & spread_all$ARID1A_KO > 1 & spread_all$MED12_KO > 1]
p <- ggplot(spread_all[spread_all$name %in% pthways$WP_ALLOGRAFT_REJECTION,], 
            aes(A375r, ARID1A_KO, label = name)) +
  geom_point(color = "red")
p + geom_text() + 
  labs(title = "Proteins involved in ALLOGRAFT REJECTION (WP)") +
  cowplot::theme_cowplot()+ 
  geom_hline(yintercept=0, linetype="dashed", color = "red")+ 
  geom_vline(xintercept=0, linetype="dashed", color = "red")


#test which pathways are the best correlated
cor.out<-function(ko, pathway){
  pathway_ko<-spread_all[[ko]][spread_all$name %in% pathway]
  pathway_a375r<-spread_all$A375r[spread_all$name %in% pathway]
  if (length(pathway_ko)> 3 & length(pathway_a375r)> 3) {
    cor.result<-cor.test(pathway_ko, pathway_a375r)
  }
  else{return(NULL)}
  t.stat=cor.result$estimate
  p.value=cor.result$p.value
  return(data.frame(t.stat, p.value))
}

arid1a_pathays<-lapply(pthways, function(x){cor.out(pathway = x, ko = "ARID1A_KO")})
med12_pathays<-lapply(pthways, function(x){cor.out(pathway = x, ko = "MED12_KO")})

arid1a_pathays<-bind_rows(arid1a_pathays, .id = "PATHWAY")
med12_pathays<-bind_rows(med12_pathays, .id = "PATHWAY")

arid1a_pathays$fdr.adj.p<-p.adjust(arid1a_pathays$p.value)
med12_pathays$fdr.adj.p<-p.adjust(med12_pathays$p.value)

arid1a_pathays$label<-NULL
med12_pathays$label<-NULL

arid1a_pathays$label[arid1a_pathays$fdr.adj.p<0.01 & arid1a_pathays$t.stat > 0]<-arid1a_pathays$PATHWAY[arid1a_pathays$fdr.adj.p<0.05 & arid1a_pathays$t.stat > 0]
med12_pathays$label[med12_pathays$fdr.adj.p<0.01 & med12_pathays$t.stat > 0]<-med12_pathays$PATHWAY[med12_pathays$fdr.adj.p<0.05 & med12_pathays$t.stat > 0]

to_plot<-rbind(data.frame(arid1a_pathays, EXP="ARID1A_KO"), data.frame(med12_pathays, EXP="MED12_KO"))
to_plot$log10.fdr.adj.p<--log10(to_plot$fdr.adj.p)
ggplot(to_plot, aes(x = t.stat, y = log10.fdr.adj.p )) +
  geom_point(
    size = 1.5, 
    alpha = 0.8 # It's nice to add some transparency because there may be overlap.
  ) +
  geom_text_repel(
    aes(label = label),
    family = "Poppins",
    size = 3,
    min.segment.length = 0, 
    seed = 42, 
    box.padding = 0.5,
    max.overlaps = Inf,
    arrow = arrow(length = unit(0.010, "npc")),
    nudge_x = .15,
    nudge_y = .5,
    color = "grey50"
  ) + facet_wrap(~EXP) + cowplot::theme_cowplot()

cor.plot<-function(ko, pathway){
  pathway_ko<-
  pathway_a375r<-
  pathway_a375r<-spread_all$A375r[spread_all$name %in% pathway]
  cor_pathway<-data.frame(KO=spread_all[[ko]][spread_all$name %in% pathway],
             A375r=spread_all$A375r[spread_all$name %in% pathway],
             label=spread_all$name[spread_all$name %in% pathway])
  g1<-ggplot(cor_pathway, aes(x = KO, y = A375r )) +
    geom_point(
      size = 1.5, 
      alpha = 0.8 # It's nice to add some transparency because there may be overlap.
    ) +
    geom_text_repel(
      aes(label = label),
      family = "Poppins",
      size = 3,
      min.segment.length = 0, 
      seed = 42, 
      box.padding = 0.5,
      max.overlaps = Inf,
      arrow = arrow(length = unit(0.010, "npc")),
      nudge_x = .15,
      nudge_y = .5,
      color = "grey50"
    )  + cowplot::theme_cowplot()
  return(g1)
}
cor.plot("ARID1A_KO", pthways$PID_SHP2_PATHWAY)


#Random

#SOX2 transcription facotrs
spread_all[which(spread_all$name %in% c("KDR", "FLT1", "EGFR", "VEGFC")),]
#VEGF 
spread_all[which(spread_all$name %in% c("SOX2", "ZNF536", "ZEB2")),]
#intersect
spread_all[which(spread_all$name %in% c("EPHA1", "VEGFB", "SORBS1")),]
spread_all[which(spread_all$name %in% c("NRG1", "MET", "FGF2", "HGF")),]

to_plot<-reshape2::melt(spread_all[which(spread_all$name %in% c("VEGFC", "FGFR1", "FGF2", "KDR", "FLT1", "SOX2")),])
ggplot(data=to_plot, aes(x=variable, y=value)) +
  geom_bar(stat="identity", fill="steelblue") + facet_wrap(~name) + cowplot::theme_cowplot()