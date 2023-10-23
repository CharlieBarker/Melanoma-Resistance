home_dir<-"~/Desktop/Melanoma_Resistance/"
setwd(home_dir)

library(pathview)

###enrichment of the weights #### 

weights <- get_weights(MOFAobject.trained, 
                       views = "all", 
                       as.data.frame = TRUE 
)

#factor 1 has been times by -1 so that all factors are consistent in that control is down and test is up 
weights[weights$factor == "Factor1",]$value <- weights[weights$factor == "Factor1",]$value * -1


#aggregate factors 
weights_factors<-lapply(unique(weights$factor), function(x){weights[weights$factor==x,]})
names(weights_factors)<-unique(weights$factor)
weights_factors<-lapply(weights_factors, function(x){x$protein<-gsub("\\;.*","",x$feature); return(x)})
weights_factors<-lapply(weights_factors, function(x){x$protein<-gsub("\\_.*","",x$protein); return(x)})
weights_factors<-lapply(weights_factors, function(x){x[c("protein","value")]})
#weights_factors<-lapply(weights_factors, function(x){x$value<-abs(x$value); return(x)})
weights_factors<-lapply(weights_factors, function(x){aggregate(.~protein, x, sum)})


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

MAPKview_input<-weights_factors_up_low
factor_1<-MAPKview_input$Factor1
kegg_input<-gprofiler2::gconvert(factor_1$protein)
dme_input_df<-merge(x = kegg_input, y = factor_1,
                    by.x= "input", by.y = "protein")
dme_input_list<-dme_input_df$value
names(dme_input_list)<-dme_input_df$target
# Produce the native KEGG plot (PNG)
pathView<-"./results/mofa/pathView/"
setwd(pathView)
dme <- pathview(gene.data=dme_input_list,
                pathway.id="hsa01521", gene.idtype=gene.idtype.list[3])

setwd(home_dir)

