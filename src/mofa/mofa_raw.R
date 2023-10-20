reticulate::use_condaenv("py37", required = T)
setwd(dir = "~/phd/MelanomaProject/mofa")

library(data.table)
library(MOFA2)
library(dplyr)
library(purrr)
library(stringr)
library(ggplot2)

list_of_inputs<-list(#signalome_view=data.frame(read.csv("input_data/signalome_by_sample.csv")),
  phospho=data.frame(read.csv("input_data/phosphosites.csv"))
  ,protein=data.frame(read.csv("input_data/proteins.csv"))
  ,mRNA=data.frame(read.csv("input_data/rna_expression.csv"))
  #,kinase_view=data.frame(read.csv("input_data/kin_matrix.csv"))
  #,tf_view=data.frame(read.csv("input_data/tf_activitties.csv"))
  )


list_of_inputs<-lapply(list_of_inputs, function(x){rownames(x)<-x$X; return(x)})
list_of_inputs<-lapply(list_of_inputs, function(x){x$X<-NULL; return(x)})
list_of_inputs<-lapply(list_of_inputs, function(x){colnames(x)<-sub("*\\.[0-9]", "", colnames(x)); return(x)})

#anova to find the most important TF features
library(PhosR)
aov <- matANOVA(mat=list_of_inputs$mRNA, grps=colnames(list_of_inputs$mRNA))
idx <- (aov < 0.001) & (rowSums(abs(list_of_inputs$mRNA) > 1) > 0) #idx is the filter 0.5 and aov of 0.05
list_of_inputs$mRNA <- list_of_inputs$mRNA[idx, ,drop = FALSE]

list_of_inputs$protein<-list_of_inputs$protein[rowSums(is.infinite(data.matrix(list_of_inputs$protein))) == 0,]
aov <- matANOVA(mat=list_of_inputs$protein, grps=colnames(list_of_inputs$protein))
idx <- (aov < 0.001) & (rowSums(abs(list_of_inputs$protein) > 1) > 0) #idx is the filter 0.5 and aov of 0.05
list_of_inputs$protein <- list_of_inputs$protein[idx, ,drop = FALSE]

library(EnsDb.Hsapiens.v86)
library(dplyr)
genename_df <- AnnotationDbi::select(EnsDb.Hsapiens.v86, keys = rownames(list_of_inputs$protein), keytype = "UNIPROTID", columns = "GENENAME")
new_names<-make.unique(genename_df$GENENAME[match(rownames(list_of_inputs$protein), genename_df$UNIPROTID)])
rownames(list_of_inputs$protein)[!is.na(new_names)] <- new_names[!is.na(new_names)]

list_of_inputs_df<-lapply(list_of_inputs, function(dat){reshape2::melt(as.matrix(dat))})

#scale this?
#scale this and views
#list_of_inputs_df$proteome_view<-NULL

inputs_df<-bind_rows(list_of_inputs_df, .id = "view")
colnames(inputs_df)<-c("view", "feature", "sample", "value")
inputs_df$sample<-str_replace(string = inputs_df$sample, pattern = "vemurafenib.trametinib", replacement = "vemurafenib_and_trametinib")


#so we have to aggregate the samples together
#iwould suggest to do an annova
long_mofa_input<-aggregate(.~view+feature+sample, inputs_df, mean)

MOFAobject <- create_mofa(long_mofa_input)
data_overview_plot<-plot_data_overview(MOFAobject)

data_opts <- get_default_data_options(MOFAobject)
head(data_opts)

model_opts <- get_default_model_options(MOFAobject)
model_opts$num_factors <- 4
head(model_opts)

train_opts <- get_default_training_options(MOFAobject)
head(train_opts)
data_opts$scale_views<-T
data_opts$scale_groups<-T
MOFAobject <- prepare_mofa(
  object = MOFAobject,
  data_options = data_opts,
  model_options = model_opts,
  training_options = train_opts
)

#run MOFA
outfile = file.path("test.hdf5")
MOFAobject.trained <- run_mofa(MOFAobject, outfile, use_basilisk = F)

sample_metadata <- data.frame(
  sample = samples_names(MOFAobject.trained)[[1]],
  Drug = unlist(map(str_split(samples_names(MOFAobject.trained)[[1]], pattern = "__"), 1)),
  Genetic = unlist(map(str_split(samples_names(MOFAobject.trained)[[1]], pattern = "__"), 2))
)

samples_metadata(MOFAobject.trained) <- sample_metadata


head(MOFAobject.trained@cache$variance_explained$r2_total[[1]]) # group 1
head(MOFAobject.trained@cache$variance_explained$r2_per_factor[[1]]) # group 1

variance_per_view<-plot_variance_explained(MOFAobject.trained, x="view", y="factor")
variance_heat<-plot_variance_explained(MOFAobject.trained, x="group", y="factor", plot_total = T)[[2]]

facotrs_plot<-plot_factor(MOFAobject.trained, 
                          factor = 1:4,
                          color_by = "Drug",
                          shape_by = "Genetic"
)


#WRITE PDFS

# pdf(file = "plots/mofa_variance.pdf",width=6,height=10)
# data_overview_plot
# ggpubr::ggarrange(variance_heat, variance_per_view, nrow = 2,align = "hv") #density plots per batch
# dev.off()
# pdf(file = "plots/mofa_factors.pdf",width=7,height=2)
# facotrs_plot
# dev.off()


###enrichment of the weights #### 

weights <- get_weights(MOFAobject.trained, 
                       views = "all", 
                       as.data.frame = TRUE 
)
weights_factors<-lapply(unique(weights$factor), function(x){weights[weights$factor==x,]})
names(weights_factors)<-unique(weights$factor)
weights_factors<-lapply(weights_factors, function(x){x$protein<-gsub("\\;.*","",x$feature); return(x)})
weights_factors<-lapply(weights_factors, function(x){x$protein<-gsub("\\_.*","",x$protein); return(x)})
weights_factors<-lapply(weights_factors, function(x){x[c("protein","value")]})
weights_factors<-lapply(weights_factors, function(x){aggregate(.~protein, x, sum)})
library(fgsea)
library(gprofiler2)
pathways.sets <- gmtPathways("pathway_sets/c2.cp.v2022.1.Hs.symbols.gmt") 
cancer.sets <- gmtPathways("pathway_sets/c6.all.v2022.1.Hs.symbols.gmt") 

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

C2_Pathways<-dot_plot(fgsea_enrich_df[fgsea_enrich_df$Rank<10 & fgsea_enrich_df$pathway_set == "C2_Pathways",])
sig_more_once<-names(which(table(fgsea_enrich_df[fgsea_enrich_df$pathway_set == "C6_Oncogenic_Signatures" & fgsea_enrich_df$Significance == T,]$pathway) > 1))
C6_Oncogenic_Signatures<-dot_plot(fgsea_enrich_df[ fgsea_enrich_df$pathway_set == "C6_Oncogenic_Signatures" & fgsea_enrich_df$pathway %in% sig_more_once,])

word_cloud_weight<-function(weights ,df, factor, pathway){
  word_in<-weights[[factor]][weights[[factor]]$protein %in% unlist(df[df$pathway == pathway & df$Factor == factor,]$leadingEdge),]
  word_in$freq<-abs(word_in$value)
  word_in$freq <- round(scales::rescale(word_in$freq, to = c(1, 50)))
  return(wordcloud(words = word_in$protein, freq = word_in$freq,max.words=200, random.order=FALSE, rot.per=0.35, colors=brewer.pal(8, "Dark2")))
}
#word_cloud_weight(df = fgsea_enrich_df,weights = weights_factors, factor = "Factor1", pathway = "PID_CDC42_PATHWAY")

# pdf(file = "plots/mofa_raw.pdf",width=18,height=10)
# library(wordcloud)
# ggpubr::ggarrange(C2_Pathways,
#                   C6_Oncogenic_Signatures,
#                   nrow = 1,widths = c(3.5,2)) #density plots per batch
# facotrs_plot
# ggpubr::ggarrange(variance_heat, variance_per_view, nrow = 2,align = "hv") #density plots per batch
# dev.off()

#now how are these things are interacting based on network analysis?
#for netextract - or Pheugo or whatever... 

#get the absolute max logFC phosphosite for each protein
weights_factors_up_low<-lapply(weights_factors, function(x){
  x[x$value < unname(quantile(x$value, probs=c(0.05, 0.95))[1]) |
      x$value > unname(quantile(x$value, probs=c(0.05, 0.95))[2]),] 
})
weights_factors_up_low<-lapply(weights_factors_up_low, function(x){x[order(x$value),]})
library(igraph)
resnikbma<-read.delim("~/phd/MelanomaProject/data/NetExtract/semSimNets/resnikbma.txt", header = F)
resnikbma_g<-graph_from_data_frame(resnikbma)

weights_factors_up_low<-lapply(weights_factors_up_low, function(x){
  genename_df <- AnnotationDbi::select(EnsDb.Hsapiens.v86, keys = x$protein, keytype = "GENENAME", columns = "UNIPROTID");
  genename_df <- genename_df[genename_df$UNIPROTID %in% V(resnikbma_g)$name,]
  x$protein <- genename_df$UNIPROTID[match(x$protein, genename_df$GENENAME)]
  return(x)})
weights_factors_up_low<-lapply(weights_factors_up_low, function(x){x[!is.na(x$protein),]})
weights_factors_up_low<-lapply(weights_factors_up_low, function(x){x$padj<-0.01;return(x)})
weights_factors_up_low<-lapply(weights_factors_up_low, function(x){x[,c("value", "protein", "padj")]})

library(readr)
lapply(names(weights_factors_up_low),
       function(x){write_delim(x = weights_factors_up_low[[x]], file = paste0("./factor_raw_more/", x, ".txt"), delim = "\t", col_names = F)})

library(pathview)
MAPKview_input<-weights_factors_up_low
factor_1<-MAPKview_input$Factor4
kegg_input<-gprofiler2::gconvert(factor_1$protein)
dme_input_df<-merge(x = kegg_input, y = factor_1,
                 by.x= "input", by.y = "protein")
dme_input_list<-dme_input_df$value
names(dme_input_list)<-dme_input_df$target
# Produce the native KEGG plot (PNG)
dme <- pathview(gene.data=dme_input_list,
                pathway.id="hsa04010", species = kegg_organism, gene.idtype=gene.idtype.list[3])
