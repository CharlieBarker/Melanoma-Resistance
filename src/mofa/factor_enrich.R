
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
weights_factors<-lapply(weights_factors, function(x){x[c("protein","value")]})
weights_factors<-lapply(weights_factors, function(x){x$value<-abs(x$value); return(x)})
weights_factors<-lapply(weights_factors, function(x){aggregate(.~protein, x, sum)})

# library(fgsea)
# library(gprofiler2)
# pathways.sets <- gmtPathways("./data/pathway_sets/c2.cp.v2022.1.Hs.symbols.gmt") 
# cancer.sets <- gmtPathways("./data/pathway_sets/c6.all.v2022.1.Hs.symbols.gmt") 
# 
# fgsea_enrich<-list()
# for (factor in names(weights_factors)) {
#   cat(factor)
#   variable<-weights_factors[[factor]]
#   ranks<-variable$value
#   names(ranks)<-variable$protein
#   ranks<-sort(ranks)
#   fgseaRes <- fgsea(pathways = pathways.sets, 
#                     stats    = ranks,
#                     minSize  = 15,
#                     maxSize  = 500)
#   fgseaRes_cancer <- fgsea(pathways = cancer.sets, 
#                            stats    = ranks,
#                            minSize  = 15,
#                            maxSize  = 500)
#   fgseaRes<-fgseaRes[order(fgseaRes$padj),]
#   fgseaRes$Rank<-seq(1:nrow(fgseaRes))
#   fgseaRes$Significance<-fgseaRes$padj < 0.1
#   
#   fgseaRes_cancer<-fgseaRes_cancer[order(fgseaRes_cancer$padj),]
#   fgseaRes_cancer$Rank<-seq(1:nrow(fgseaRes_cancer))
#   fgseaRes_cancer$Significance<-fgseaRes_cancer$padj < 0.05
#   
#   fgsea_enrich[[factor]]<- rbind(data.frame(fgseaRes, pathway_set = "C2_Pathways"),
#                                  data.frame(fgseaRes_cancer, pathway_set = "C6_Oncogenic_Signatures"))
# }
# fgsea_enrich_df<-bind_rows(fgsea_enrich, .id = "Factor")
# dot_plot<-function(df){
#   fig_out<-ggplot(df, aes(x=Factor, y = reorder(pathway, -NES), color = NES, size = -log10(padj), shape=Significance)) +
#     geom_point() +
#     scale_color_viridis_c(name = 'Combined.Score') +
#     cowplot::theme_cowplot() +
#     theme(axis.line  = element_blank()) +
#     theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
#     ylab('') +
#     theme(axis.ticks = element_blank())  +
#     theme(axis.ticks = element_blank(),
#           panel.border = element_rect(colour = "black", fill=NA, size=2)) +
#     scale_shape_manual(values=c(1,16) )
#   return(fig_out)
# }
# 
# C2_Pathways<-dot_plot(fgsea_enrich_df[fgsea_enrich_df$Rank<5 & fgsea_enrich_df$pathway_set == "C2_Pathways",])

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

# library(pathview)
# MAPKview_input<-weights_factors_up_low
# factor_1<-MAPKview_input$Factor4
# kegg_input<-gprofiler2::gconvert(factor_1$protein)
# dme_input_df<-merge(x = kegg_input, y = factor_1,
#                     by.x= "input", by.y = "protein")
# dme_input_list<-dme_input_df$value
# names(dme_input_list)<-dme_input_df$target
# # Produce the native KEGG plot (PNG)
# dme <- pathview(gene.data=dme_input_list,
#                 pathway.id="hsa04010", species = kegg_organism, gene.idtype=gene.idtype.list[3])
