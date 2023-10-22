setwd(dir = "~/phd/MelanomaProject/mofa")

library(fgsea)
library(gprofiler2)
library(ggrepel)
library(MOFA2)
library(purrr)
library(stringr)
library(RColorBrewer)
library(tidyverse)

mofa<-load_model(file.path("~/phd/MelanomaProject/mofa/test.hdf5"))
#add meta data 
sample_metadata <- data.frame(
  sample = samples_names(mofa)[[1]],
  Drug = unlist(map(str_split(samples_names(mofa)[[1]], pattern = "__"), 1)),
  Genetic = unlist(map(str_split(samples_names(mofa)[[1]], pattern = "__"), 2))
)

samples_metadata(mofa) <- sample_metadata

plot_factors(mofa, factors = 1:4,
             color_by = "Drug",
             shape_by = "Genetic")


#get weights
weights <- get_weights(mofa, 
                       views = "all", 
                       as.data.frame = TRUE 
)
weights_factors<-lapply(unique(weights$factor), function(x){weights[weights$factor==x,]})
names(weights_factors)<-unique(weights$factor)
ordered_weights<-lapply(weights_factors, function(x){ordered_weights<-x[order(x$value),]
ordered_weights$rank<-c(1:nrow(ordered_weights)); return(ordered_weights)})
ordered_weights<-lapply(ordered_weights, function(x){x$protein<-gsub("\\;.*","",x$feature); return(x)})
ordered_weights<-lapply(ordered_weights, function(x){x$protein<-gsub("\\_.*","",x$protein); return(x)})
complete_weights<-bind_rows(ordered_weights)
#create files for omics analyser 
mofa_weights<-complete_weights
mofa_weights$psite<-map(str_split(mofa_weights$feature, ";"), 2)
mofa_weights$psite[unlist(lapply(mofa_weights$psite, function(x){length(x)==0}))]<-""
mofa_weights$psite<-unlist(mofa_weights$psite)

weight_df<-data.frame(name=mofa_weights$protein,
                      site=mofa_weights$psite, #if applicable
                      factor=mofa_weights$factor,
                      weights=mofa_weights$value,
                      view=mofa_weights$view)
#add the different lfcs
spread_weight <- weight_df %>% spread(view, weights)
library(dplyr)
spread_weight<-aggregate(. ~ name + factor, data=spread_weight, FUN=na.omit, na.action="na.pass")

factors_list<-lapply(unique(spread_weight$factor), function(x){spread_weight[spread_weight$factor == x,]})
names(factors_list)<-unique(spread_weight$factor)

#remove empty psites 
factors_list<-lapply(factors_list, function(x){x$site[x$site == ""]<-"-"; return(x)})
#remove empty psites in lists
factors_list<-lapply(factors_list, function(x){
  x$site<-lapply(x$site, function(x){x[x != ""]})
  return(x)
})
#concat the psites where there are more than one on a protein 
factors_list<-lapply(factors_list, function(x){
  x$site<-unlist(lapply(x$site, paste, collapse="-"))
  return(x)
})
#get the absolute sum of all the psites where there is more than one on a protein. 
factors_list<-lapply(factors_list, function(x){
  x$phospho<-lapply(x$phospho, as.numeric)
  return(x)
})
factors_list<-lapply(factors_list, function(x){
  x$phospho<-unlist(lapply(x$phospho, function(x){sum(abs(x))}))
  return(x)
})

#make sure the remaining attributes are as.numeric

factors_list<-lapply(factors_list, function(x){
  #replace elements without measurement with 0
  x$kinase_view[lengths(x$kinase_view) == 0] <- 0
  x$mRNA[lengths(x$mRNA) == 0] <- 0
  x$protein[lengths(x$protein) == 0] <- 0
  x$tf_view[lengths(x$tf_view) == 0] <- 0
  #ad.numeric
  x$kinase_view<-unlist(lapply(x$kinase_view, function(x){unlist(as.numeric(x))}))
  x$mRNA<-unlist(lapply(x$mRNA, function(x){unlist(as.numeric(x))}))
  x$protein<-unlist(lapply(x$protein, function(x){unlist(as.numeric(x))}))
  x$tf_view<-unlist(lapply(x$tf_view, function(x){unlist(as.numeric(x))}))
  #replace elements without measurement with NA
  x$kinase_view[x$kinase_view == 0] <- 0
  x$mRNA[x$mRNA == 0] <- 0
  x$protein[x$protein == 0] <- 0
  x$tf_view[x$tf_view == 0] <- 0
  x$phospho[x$phospho == 0] <- 0
  #round to 2 digits (indistinguishable to our eyes anyway in terms of plotting colour)
  x$kinase_view<-round(x$kinase_view, digits = 2)
  x$mRNA<-round(x$mRNA, digits = 2)
  x$protein<-round(x$protein, digits = 2)
  x$tf_view<-round(x$tf_view, digits = 2)
  x$phospho<-round(x$phospho, digits = 2)
  
  return(x)
})

library(readr)
lapply(names(factors_list),
       function(x){write_delim(x = factors_list[[x]], file = paste0("./omics_analyser_inpt/", x, ".txt"), delim = "\t", col_names = T)})

write.csv(complete_weights,file = "model_weights_raw_and_processed.csv")

plot_labelled_factor_weights<-function(factor_list, 
                                       factor_to_plot,
                                       quartile_vec=c(0.005, 0.995)){
  colour_vec= c("kinase_view" = "#D00000", 
                "phospho_view" = "#3F88C5",
                "proteome_view" = "#FFBA08",
                "transcriptome_view" = "#136F63")
  to_plot<-factor_list[[factor_to_plot]]
  quantile_range<-quantile(to_plot$value, quartile_vec)
  bottom_lim<-unname(quantile_range[1])
  top_lim<-unname(quantile_range[2])
  to_plot$label<-""
  to_label<- to_plot$value < bottom_lim | to_plot$value > top_lim
  to_plot$label[to_label]<-as.character(to_plot$feature[to_label])
  # Change the point size, and shape
  g1<-ggplot(to_plot, aes(x=value, y=rank,
                          fill=view, 
                          color=view, 
                          label=label)) +
    geom_point(size=2, shape=23) + cowplot::theme_cowplot()+
    geom_vline(xintercept=bottom_lim, linetype='dotted', col = 'red')+
    geom_vline(xintercept=top_lim, linetype='dotted', col = 'red')  +
    geom_text_repel(force_pull   = 20, # do not pull toward data points
                    nudge_y      = 0.6,
                    direction    = "x",
                    angle        = 90,
                    hjust        = 3,
                    vjust        =  0.8,
                    segment.size = 0.2,
                    max.iter = 1e4, 
                    max.time = 1,max.overlaps = 10,
                    size = 3)+
    scale_color_manual(values = colour_vec) +
    scale_fill_manual(values = colour_vec)
  return(g1)
}

pdf(file = "plots/Factor_weights.pdf",width=10,height=14)
ggpubr::ggarrange(#Factor 1
                  plot_factor(mofa, factor = 1, color_by = "Drug", shape_by = "Genetic") + 
                    theme(legend.position="none")+ coord_flip()+
                    theme(axis.text.x=element_blank(), #remove x axis labels
                          axis.ticks.x=element_blank(), #remove x axis ticks
                          axis.text.y=element_blank(),  #remove y axis labels
                          axis.ticks.y=element_blank()  #remove y axis ticks
                    ),         
                  plot_labelled_factor_weights(factor_list = ordered_weights, factor_to_plot = "Factor1") + 
                    theme(legend.position="none"),
                  #Factor 2 
                  plot_factor(mofa, factor = 2, color_by = "Drug", shape_by = "Genetic") + 
                    theme(legend.position="none")+ coord_flip()+
                    theme(axis.text.x=element_blank(), #remove x axis labels
                          axis.ticks.x=element_blank(), #remove x axis ticks
                          axis.text.y=element_blank(),  #remove y axis labels
                          axis.ticks.y=element_blank()  #remove y axis ticks
                    ),      
                  plot_labelled_factor_weights(factor_list = ordered_weights, factor_to_plot = "Factor2") + 
                    theme(legend.position="none"),
                  #Factor 3
                  plot_factor(mofa, factor = 3, color_by = "Drug", shape_by = "Genetic") + 
                    theme(legend.position="none")+ coord_flip()+
                    theme(axis.text.x=element_blank(), #remove x axis labels
                          axis.ticks.x=element_blank(), #remove x axis ticks
                          axis.text.y=element_blank(),  #remove y axis labels
                          axis.ticks.y=element_blank()  #remove y axis ticks
                    ),  
                  plot_labelled_factor_weights(factor_list = ordered_weights, factor_to_plot = "Factor3") + 
                    theme(legend.position="none"),
                  #Factor 4
                  plot_factor(mofa, factor = 4, color_by = "Drug", shape_by = "Genetic") + 
                    theme(legend.position="none")+ coord_flip()+
                    theme(axis.text.x=element_blank(), #remove x axis labels
                          axis.ticks.x=element_blank(), #remove x axis ticks
                          axis.text.y=element_blank(),  #remove y axis labels
                          axis.ticks.y=element_blank()  #remove y axis ticks
                    ),                  
                  plot_labelled_factor_weights(factor_list = ordered_weights, factor_to_plot = "Factor4") + 
                    theme(legend.position="none"),
                  #The overall parameters 
                  ncol = 1, heights = rep(c(1,5),4), align = "hv")

dev.off()

