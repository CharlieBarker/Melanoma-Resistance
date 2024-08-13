
#find correct environment 
packLib="/usr/lib/R"
if (file.exists(packLib)) {
  reticulate::use_condaenv("py37", required = T)
  setwd(dir = "~/Desktop/Melanoma_Resistance/")
}else {
  setwd(dir = "~/Desktop/Melanoma_Resistance/")
}
source("./src/functions/default_variables.R")

library(ggrepel)
library(ggpubr)
library(wesanderson)
library(gridExtra)

#empty plot for the networ
empty_plot <- ggplot() + geom_blank() + cowplot::theme_cowplot()

enrichr_df_process<-function(df, labels_to_show, title){
  df$labels <- df$term
  df$labels[!df$labels %in% labels_to_show] <- ""
  scale_0_to_1 <- function(x) {
    scaled <- (x - min(x)) / (max(x) - min(x))
    return(scaled)
  }
  df$shade <- scale_0_to_1(df$odds_ratio) * scale_0_to_1(df$log_pval) 
  pal <- wes_palette("Zissou1", 100, type = "continuous")
  volcano_plot<-ggplot(df, aes(x=odds_ratio, y=log_pval, label = labels, colour = shade)) +
    geom_point() + 
    cowplot::theme_cowplot() + 
    geom_text_repel(colour="black", size = 2.5, min.segment.length = 0, seed = 42, box.padding = 0.5) + 
    scale_color_gradientn(colours = pal) +
    labs(x = "Odds Ratio", y = "Log10(Adjusted P Value)") +
    theme(legend.position="none",
          plot.title=element_text(size=7)) +
    ggtitle(title)
  return(volcano_plot)
}

labels_to_show_up<-c("S6K1 signaling", "CBL-mediated ligand-induced downregulation of EGF receptors",
                     "Mammalian target of rapamycin complex 1 (mTORC1)-mediated signaling",
                     "ERBB signaling pathway", "Signaling events mediated by hepatocyte growth factor receptor (c-Met)",
                     "mTOR signaling pathway", 'PI3K events in ERBB2 signaling', 'IGF1 pathway')
labels_to_show_down<-c("Immune system", "Immune system signaling by interferons, interleukins, prolactin, and growth hormones",
                       "MAPK signaling pathway",
                       "TNF-alpha signaling pathway", "ATM-mediated phosphorylation of repair proteins",
                       "JNK (c-Jun kinases) phosphorylation and  activation mediated by activated human TAK1",
                       'RSK activation',"Recruitment of repair and signaling proteins to double-strand breaks")

g1 <- empty_plot
g2 <- grid.arrange(enrichr_df_process(read.csv(file = "./results/phuego/enrichr_dfs/factor2_up.csv"), labels_to_show_up,
                                      title="Enriched terms in changes associated with Combination therapy (Upregulated network)"), 
                   enrichr_df_process(read.csv(file = "./results/phuego/enrichr_dfs/factor2_down.csv"), labels_to_show_down, 
                                      title="Enriched terms in changes associated with Combination therapy (Downregulated network)"),
                   ncol = 1)

pdf(# The directory you want to save the file in
  width = 15, # The width of the plot in inches
  height = 10,
  file = "./paper/Figures/combination_figure.pdf")


grid.arrange(g1, g2, nrow = 1, widths = c(.66, .33))

dev.off()