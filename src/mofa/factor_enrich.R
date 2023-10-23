library(MOFA2)
library(dplyr)
library(EnsDb.Hsapiens.v86)


setwd("~/Desktop/Melanoma_Resistance//")

MOFAobject.trained<-load_model(file = "./results/mofa/mofa_object.hdf5")

###enrichment of the weights #### 

weights <- get_weights(MOFAobject.trained, 
                       views = "all", 
                       as.data.frame = TRUE 
)

#aggregate factors 
weights_factors<-lapply(unique(weights$factor), function(x){weights[weights$factor==x,]})
names(weights_factors)<-unique(weights$factor)
weights_factors<-lapply(weights_factors, function(x){x$protein<-gsub("\\;.*","",x$feature); return(x)})
weights_factors<-lapply(weights_factors, function(x){x$protein<-gsub("\\_.*","",x$protein); return(x)})
#get value of most strongly changed psite 
weights_factors<-lapply(weights_factors, function(x){
  # Create an aggregated dataframe for "phospho" view
  aggregated_df <- x %>%
    dplyr::filter(view == "phospho") %>%  # Filter for "phospho" view
    group_by(protein) %>%  # Group by protein
    dplyr::slice(which.max(abs(value))) %>%  # Select the row with the highest absolute factor value
    ungroup()
  
  # Replace the original dataframe with the aggregated results
  x <- x %>%
    dplyr::filter(view != "phospho") %>%  # Remove "phospho" view rows
    bind_rows(aggregated_df)  # Append the aggregated results
  return(x)
})
weights_factors<-lapply(weights_factors, function(x){x[c("protein","value")]})
#now aggregate with everything else
weights_factors<-lapply(weights_factors, function(x){aggregate(.~protein, x, sum)})

library(fgsea)
library(gprofiler2)
pathways.sets <- gmtPathways("./data/pathway_sets/c2.cp.v2022.1.Hs.symbols.gmt")
cancer.sets <- gmtPathways("./data/pathway_sets/c6.all.v2022.1.Hs.symbols.gmt")

fgsea_enrich<-list()
for (factor in names(weights_factors)) {
  cat(factor)
  variable<-weights_factors[[factor]]
  ranks<-variable$value
  names(ranks)<-variable$protein
  ranks<-sort(ranks)
  fgseaRes <- fgsea(pathways = pathways.sets,
                    stats    = ranks,
                    minSize  = 15,
                    maxSize  = 500)
  fgseaRes_cancer <- fgsea(pathways = cancer.sets,
                           stats    = ranks,
                           minSize  = 15,
                           maxSize  = 500)
  fgseaRes<-fgseaRes[order(fgseaRes$padj),]
  fgseaRes$Rank<-seq(1:nrow(fgseaRes))
  fgseaRes$Significance<-fgseaRes$padj < 0.1

  fgseaRes_cancer<-fgseaRes_cancer[order(fgseaRes_cancer$padj),]
  fgseaRes_cancer$Rank<-seq(1:nrow(fgseaRes_cancer))
  fgseaRes_cancer$Significance<-fgseaRes_cancer$padj < 0.05

  fgsea_enrich[[factor]]<- rbind(data.frame(fgseaRes, pathway_set = "C2_Pathways"),
                                 data.frame(fgseaRes_cancer, pathway_set = "C6_Oncogenic_Signatures"))
}
fgsea_enrich_df<-bind_rows(fgsea_enrich, .id = "Factor")
dot_plot<-function(df){
  fig_out<-ggplot(df, aes(x=Factor, y = reorder(pathway, -NES), color = NES, size = -log10(padj), shape=Significance)) +
    geom_point() +
    scale_color_viridis_c(name = 'Combined.Score') +
    cowplot::theme_cowplot() +
    theme(axis.line  = element_blank()) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
    ylab('') +
    theme(axis.ticks = element_blank())  +
    theme(axis.ticks = element_blank(),
          panel.border = element_rect(colour = "black", fill=NA, size=2)) +
    scale_shape_manual(values=c(1,16) )
  return(fig_out)
}

C2_Pathways<-dot_plot(fgsea_enrich_df[fgsea_enrich_df$Rank<15 & fgsea_enrich_df$pathway_set == "C2_Pathways",])

#now how are these things are interacting based on network analysis?
#for netextract - or Pheugo or whatever... 

#get the absolute max logFC phosphosite for each protein
weights_factors_up_low<-lapply(weights_factors, function(x){
  x[x$value < unname(quantile(x$value, probs=c(0.05, 0.95))[1]) |
      x$value > unname(quantile(x$value, probs=c(0.05, 0.95))[2]),] 
})
weights_factors_up_low<-lapply(weights_factors_up_low, function(x){x[order(x$value),]})

weights_factors_up_low<-lapply(weights_factors_up_low, function(x){
  genename_df <- AnnotationDbi::select(EnsDb.Hsapiens.v86, keys = x$protein, keytype = "GENENAME", columns = "UNIPROTID");
  x$protein <- genename_df$UNIPROTID[match(x$protein, genename_df$GENENAME)]
  return(x)})
weights_factors_up_low<-lapply(weights_factors_up_low, function(x){x[!is.na(x$protein),]})
weights_factors_up_low<-lapply(weights_factors_up_low, function(x){x$padj<-0.01;return(x)})
weights_factors_up_low<-lapply(weights_factors_up_low, function(x){x[,c("protein", "value")]})

library(readr)
lapply(names(weights_factors_up_low),
       function(x){write_delim(x = weights_factors_up_low[[x]], file = paste0("./results/mofa/mofa_factors/", x, ".txt"), delim = "\t", col_names = F)})


