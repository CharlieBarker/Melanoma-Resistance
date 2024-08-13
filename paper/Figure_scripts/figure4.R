
#find correct environment 
packLib="/usr/lib/R"
if (file.exists(packLib)) {
  reticulate::use_condaenv("py37", required = T)
  setwd(dir = "~/Desktop/Melanoma_Resistance/")
}else {
  setwd(dir = "~/Desktop/Melanoma_Resistance/")
}
source("./src/functions/default_variables.R")


library(igraph)
library(EnsDb.Hsapiens.v86)
library("ggVennDiagram")
library(MOFA2)
library(stringr)
library(purrr)
library(ggrepel)
library(OmnipathR)
library(magrittr)
library(dplyr)
library(tidyr)
library(ggpubr)
library(wesanderson)
library(ggplot2)
library(cowplot)

#globally improtant variables

Abs<-F #are we looking at up and down regulation seperately?



readr::local_edition(1) 

#get mofa weights 
MOFAobject.trained<-load_model(file = "./results/mofa/mofa_object.hdf5")

weights <- get_weights(MOFAobject.trained, 
                       views = "all", 
                       as.data.frame = TRUE 
)
weights$node<-unlist(map(str_split(weights$feature, pattern = "_"),1))
weights$node<-unlist(map(str_split(weights$node, pattern = ";"),1))

#get phuego graphs
KDE <- "0.5"
factorS <- c("Factor1", "Factor2", "Factor3")
results_dir <- "./results/phuego/results/"
factor_centrality<-list()

#provisionally we need to find a way to make the phuego networks directed. 

for (factor in factorS) {
  file_graphml_up<-paste0(results_dir, factor, "/increased/KDE_", KDE, "/networks/KDE.graphml")
  file_graphml_down<-paste0(results_dir, factor, "/decreased/KDE_", KDE, "/networks/KDE.graphml")
  up_graph<-read_graph(file=file_graphml_up,format = "graphml")
  centrality_df<-rbind(data.frame(centrality = degree(read_graph(file=file_graphml_up,format = "graphml")), 
                                  direction = "Up regulated"),
             data.frame(centrality = degree(read_graph(file=file_graphml_down,format = "graphml")), 
                        direction = "Down regulated"))
  centrality_df$uniprot<-rownames(centrality_df)
  factor_centrality[[factor]]<-centrality_df
}

factor_centrality_df<-bind_rows(factor_centrality, .id = "factor")

genename_df <- AnnotationDbi::select(EnsDb.Hsapiens.v86, keys = factor_centrality_df$uniprot, 
                                      keytype = "UNIPROTID", 
                                      columns = "GENENAME")
factor_centrality_df$genename<-genename_df$GENENAME[match(factor_centrality_df$uniprot, genename_df$UNIPROTID)]

# Summarize the weights by sum, accross the views 
summarized_weights <- weights %>%
  group_by(factor, node) %>%
  summarise(value = sum(value)) %>%
  ungroup()

to_plot<-merge(x = factor_centrality_df, y = summarized_weights, 
      by.x=c("factor", "genename"), by.y=c("factor", "node"))
to_plot<-to_plot[to_plot$centrality!=0,]
# Your ggplot code
panelA <- ggplot(to_plot[to_plot$factor=="Factor3",], aes(x=value, y=centrality, color=direction, label=genename)) +
  geom_point() + 
  geom_text_repel(data = to_plot[to_plot$factor=="Factor3",], color="black",
                  aes(label=ifelse(xor(abs(value) > 1, centrality > 3), genename, '')),
                  min.segment.length = 0, seed = 42, box.padding = 0.5, nudge_y = 0.1) +
  geom_vline(xintercept = 0, linetype = "dotted", color = "red") +  # Add dotted red line at x = 0
  labs(x = "Factor 3 weights (ARID1A KO)", y = "Degree in inferred signalling network") +  # Add axis titles
  theme_cowplot()+ scale_y_continuous(trans = 'log2') +
  scale_color_manual(values = wes_palette("Darjeeling1")) + grids()


list_of_inputs<-list(#signalome_view=data.frame(read.csv("input_data/signalome_by_sample.csv")),
  phospho=data.frame(read.csv("./data/input_data/phosphosites.csv"))
  ,protein=data.frame(read.csv("./data/input_data/proteins.csv"))
  ,mRNA=data.frame(read.csv("./data/input_data/rna_expression.csv"))
)


gene_ids <- list(
  FGFR1 = "P11362",
  FGFR2 = "P21802",
  FGF1 = "P05230",
  FGF2 = "P09038",
  FLT1 = "P17948",
  STK1 = "Q13188",
  ERBB3 = "P21860",
  KDR = "P35968",
  DUSP10 = "Q13115",
  MAP3K5 = "Q99683",
  PRKCD = "Q05655",
  BDNF = "P23560",
  FL1 = "test",
  EGFR="test"
  )


all_of_interest <- names(gene_ids)
rtk_mrna<-reshape2::melt(list_of_inputs$mRNA[list_of_inputs$mRNA$X %in% all_of_interest,])
rtk_mrna$X <- factor(rtk_mrna$X, levels=all_of_interest)
rtk_mrna$variable<-sub("*\\.[0-9]", "", rtk_mrna$variable)
rtk_mrna <- rtk_mrna |>
  separate_wider_delim(variable, delim = "__", names = c("drug", "ko"))


rtk_mrna$drug<-names(replacement_Vec)[match(rtk_mrna$drug, unname(replacement_Vec))]
# Define the desired order of levels
desired_order <- c("Untreated", "Vemurafenib", "Trametinib", "Combinations")
# Convert my_column to a factor with the specified order
rtk_mrna$drug <- factor(rtk_mrna$drug, levels = desired_order)

output_file<-"./paper/plots/phuego/factor3_rna_abundances.pdf"

panelB <- rtk_mrna %>%
  ggplot(aes(x=ko, y=value, fill=drug)) +
  geom_boxplot(position=position_dodge(width=0.75)) +  # Adjust width as needed
  geom_jitter(position=position_dodge(width=0.75), color="black", size=0.4, alpha=0.9) +  # Adjust width as needed
  scale_fill_manual(values = drug_colours) +
  cowplot::theme_cowplot() +
  theme(
    plot.title = element_text(size=11),
    axis.text.x = element_text(angle = 70, hjust = 1, size = rel(1)),
  ) +
  xlab("") +
  facet_wrap(~X, scales = "free_y", nrow = 3) + 
  grids(linetype = "dashed") +
  labs(
    x = "Genetic background",
    y = "RNA abundance"
  )

#get tf activities 
arid1a_tf_acts<-read.csv(file = "./results/transcriptomics/tf_activity/arid1a_tf_acts.csv")
arid1a_tf_pval<-read.csv(file = "./results/transcriptomics/tf_activity/arid1a_tf_pval.csv")
arid1a_tf<-merge(arid1a_tf_acts,arid1a_tf_pval,by="X")
colnames(arid1a_tf)<-c("TF", "Activity", "P_val")
arid1a_tf$diffexpressed <- arid1a_tf$Activity < 0
arid1a_tf$diffexpressed[arid1a_tf$P_val>0.1] <- "Not significant"
arid1a_tf$diffexpressed [arid1a_tf$diffexpressed == T] <- "Downregulated"
arid1a_tf$diffexpressed [arid1a_tf$diffexpressed == F] <- "Upregulated"
arid1a_tf$label<-arid1a_tf$TF
arid1a_tf$TF[arid1a_tf$diffexpressed=="Not significant"]<- ""

# Your ggplot code
panelC <- ggplot(arid1a_tf, aes(x=Activity, y=-log10(P_val), color=diffexpressed, label=label)) +
  geom_point() + 
  geom_text_repel(color="black",
                  min.segment.length = 0, seed = 42, box.padding = 0.5, nudge_y = 0.1) +
  geom_vline(xintercept = 0, linetype = "dotted", color = "red") +  # Add dotted red line at x = 0
  labs(x = "TF activity") +  # Add axis titles
  theme_cowplot() +
  scale_color_manual(values = wes_palette("Darjeeling1")) + grids()

#panel C - RFX5, RFXAP, CIITA and TWIST1
#collectri<-OmnipathR::collectri(organism=9606L, genesymbols=TRUE, loops=TRUE)
tfs_of_interest<-collectri[collectri$source_genesymbol %in% c("RFX5", "RFXAP", "CIITA", "TWIST1"),]

#print lfc for genes in tf regulons ARID1A - so that we can generate figure 4
arid1a_gene_lfc<-read.csv(file = "./results/transcriptomics/arid1a_lfc.csv")


vol_in<-rbind(data.frame(arid1a_gene_lfc, TF="RFX5", in_regulon=arid1a_gene_lfc$gene_symbol %in% tfs_of_interest$target_genesymbol[tfs_of_interest$source_genesymbol=="RFX5"]),
              data.frame(arid1a_gene_lfc, TF="RFXAP", in_regulon=arid1a_gene_lfc$gene_symbol %in% tfs_of_interest$target_genesymbol[tfs_of_interest$source_genesymbol=="RFXAP"]),
              data.frame(arid1a_gene_lfc, TF="CIITA", in_regulon=arid1a_gene_lfc$gene_symbol %in% tfs_of_interest$target_genesymbol[tfs_of_interest$source_genesymbol=="CIITA"]),
              data.frame(arid1a_gene_lfc, TF="TWIST1", in_regulon=arid1a_gene_lfc$gene_symbol %in% tfs_of_interest$target_genesymbol[tfs_of_interest$source_genesymbol=="TWIST1"]))
vol_in$X<-NULL
vol_in$label<-""
vol_in$label[vol_in$in_regulon]<-vol_in$gene_symbol[vol_in$in_regulon]
vol_in<-vol_in[vol_in$in_regulon,]

# Define alpha values based on 'below'
alpha_values <- ifelse(!vol_in$in_regulon, 0.05, 1)


# Create labels based on conditions
vol_in$label <- ifelse((vol_in$log2FoldChange >= -0.5 & vol_in$log2FoldChange <= 0.5) | vol_in$padj >= 0.1, "", vol_in$label)
vol_in$colour <- ifelse((vol_in$log2FoldChange >= -0.5 & vol_in$log2FoldChange <= 0.5) | vol_in$padj >= 0.1, "darkgrey", "darkred")

# Your existing ggplot code
panelD<- ggplot(vol_in, 
                       aes(log2FoldChange, -log10(padj), color = colour, label=label)) +
  geom_point() + 
  geom_vline(xintercept = c(-0.5, 0.5), linetype = "dashed") + # Vertical lines
  geom_hline(yintercept = 1, linetype = "dashed") +        # Horizontal line
  ylab("Log10(Adjusted P value)") + xlab("Log fold change") + facet_wrap(~TF, scales = "free") +
  cowplot::theme_cowplot()+
  geom_text_repel(min.segment.length = 0, seed = 42, box.padding = 0.5, colour="black") + 
  theme(legend.position="none") +
  scale_colour_manual(values = wes_palette("Royal1"))


pdf(# The directory you want to save the file in
  width = 20, # The width of the plot in inches
  height = 15,
  file = "./paper/Figures/ARID1a_figure.pdf")
# Showing plot
ggarrange(panelA, panelB, panelC,panelD, ncol = 2, nrow = 2, labels = "AUTO")
dev.off()
