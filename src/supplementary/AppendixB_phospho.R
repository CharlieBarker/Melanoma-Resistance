# load libraries (this gets ggplot2 and related libraries)
library(tidyverse)
library(data.table)
# these are from Bioconductor
library(limma) 
library(edgeR) 
library(sva)
# we use this for the pairwise correlation plots
library(psych)
suppressPackageStartupMessages({
  library(PhosR)
})
library(EnsDb.Hsapiens.v86)

#plot a 3d pca
pca_3D<-function(data, colour, out_file){
  dat_pc<-t(data.frame(na.omit(data)))
  pca <- prcomp(dat_pc)
  scores = as.data.frame(pca$x) 
  scores$colour <- as.factor(design[match(rownames(scores), design$code),colour])
  palette(rainbow(length(unique(scores$colour))))
  p1<-plot(PC1~PC3, scores, col = scores$colour)
  library(rgl)
  par3d(cex=.4)
  with(scores, plot3d(PC1,PC2,PC3, col = as.numeric(scores$colour), size = 5))
  legend3d("topright", legend = levels(scores$colour), col = as.factor(scores$colour), pch=19)
  rgl.postscript(out_file,fmt="pdf")
  return(p1)
}

####LOADING DATA####

path_data<-"/home/charlie/MelanomaProject/data/phosphoproteomics/main/swiss-prot/total/" #SWISS-plot
setwd(path_data)

combined_peptides<-read.csv(file = "./peptide_total.csv")
rownames(combined_peptides)<-combined_peptides$X
combined_peptides$X<-NULL
# read the design matrix
design<-read.csv(file = "~/MelanomaProject/data/phosphoproteomics/main/design.csv",header = T)
design$code<-paste0("X", design$TMT.labels,".plex", design$Plex)
design$Drug <- str_replace_all(design$Drug, pattern = "WT", replacement = "Untreated")
design$total<-paste0(design$Drug,"__", design$Genetic)
grps <- design$total[match(colnames(combined_peptides), design$code)]

####IMPUTATION####

#prep phosR data
normSLN<-function(x){
  # first basic normalization is to adjust each TMT experiment to equal signal per channel
  # figure out the global scaling value
  target <- mean(colSums(x, na.rm = T))
  # do the sample loading normalization before the IRS normalization
  # there is a different correction factor for each column
  norm_facs <- target / colSums(x, na.rm = T)
  data_sl <- sweep(x, 2, norm_facs, FUN = "*")
  return(data_sl)
}
raw_peptides<-combined_peptides
sln_peptides<-normSLN(combined_peptides)

phos_data<-sln_peptides
#extract phospho-site information
site = readr::parse_number(unlist(map(str_split(rownames(phos_data), pattern = "__"), last)))
uniprot_id = sub(" .*", "", (map(str_split(rownames(phos_data), pattern = "__"), 2))) 
genename_df<-AnnotationDbi::select(EnsDb.Hsapiens.v86, keys = uniprot_id, keytype = "UNIPROTID", columns = "GENENAME")
geneSymbol <- genename_df$GENENAME
residue = gsub('[[:digit:]]+', '', unlist(map(str_split(rownames(phos_data), pattern = "__"), last)))
sequence<-unlist(map(str_split(rownames(phos_data), "__"), 1))

counts <- table(design$total[match(colnames(phos_data), design$code)])
ppe <- PhosphoExperiment(assays = list(Quantification = as.matrix(phos_data)), 
                         Site = site,
                         GeneSymbol = geneSymbol,
                         Residue = residue,
                         Sequence = sequence)

#Next, we will perform some filtering of phosphosites so that only phosphosites with quantification for at least 50% of the replicates in at least one of the conditions are retained. 


dim(ppe)
ppe_filtered <- selectGrps(ppe, grps, 0.5, n=1) 
dim(ppe_filtered)
ppe_imputed_tmp <- scImpute(ppe_filtered, 0.5, grps)[,colnames(ppe_filtered)]
dim(na.omit(data.frame(SummarizedExperiment::assay(ppe_imputed_tmp,"imputed"))))

# Paired tail-based imputation
ppe_imputed <- ppe_imputed_tmp
#Impute between Trametinib treatments 
tram_ARID1A_log<-design$total[match(colnames(ppe_imputed), design$code)] == "Trametinib_10nM__ARID1A_KO"
tram_MED12_log<-design$total[match(colnames(ppe_imputed), design$code)] == "Trametinib_10nM__MED12_KO"
tram_WT_log<-design$total[match(colnames(ppe_imputed), design$code)] == "Trametinib_10nM__WT"

ver_ARID1A_log<-design$total[match(colnames(ppe_imputed), design$code)] == "Vermurafenib_1uM__ARID1A_KO"
var_MED12_log<-design$total[match(colnames(ppe_imputed), design$code)] == "Vermurafenib_1uM__MED12_KO"
ver_WT_log<-design$total[match(colnames(ppe_imputed), design$code)] == "Vermurafenib_1uM__WT"

WT_ARID1A_log<-design$total[match(colnames(ppe_imputed), design$code)] == "WT__ARID1A_KO"
WT_MED12_log<-design$total[match(colnames(ppe_imputed), design$code)] == "WT__MED12_KO"
WT_WT_log<-design$total[match(colnames(ppe_imputed), design$code)] == "WT__WT"

#impute is done within the drug treatments, as these are where the most variation occurs. 

#Impute between Trametinib treatments 
ppe_imputed[,tram_ARID1A_log] <- ptImpute(ppe_imputed[,tram_WT_log], #control
                                          ppe_imputed[,tram_ARID1A_log], #ARID1A_Trametinib
                                          percent1 = 0.6, percent2 = .2, paired = F)
ppe_imputed[,tram_MED12_log] <- ptImpute(ppe_imputed[,tram_WT_log], #control
                                         ppe_imputed[,tram_MED12_log], #MED12_Trametinib
                                         percent1 = 0.6, percent2 = .2, paired = F)
#Impute between Vermurafenib treatments 
ppe_imputed[,ver_ARID1A_log] <- ptImpute(ppe_imputed[,ver_WT_log], #control
                                         ppe_imputed[,ver_ARID1A_log], #ARID1A_Vermurafenib
                                         percent1 = 0.6, percent2 = .2, paired = F)
ppe_imputed[,var_MED12_log] <- ptImpute(ppe_imputed[,ver_WT_log], #control
                                        ppe_imputed[,var_MED12_log], #MED12_Vermurafenib
                                        percent1 = 0.6, percent2 = .2, paired = F)
dim(na.omit(data.frame(SummarizedExperiment::assay(ppe_imputed,"imputed"))))

#Impute between WT - care! we only have 2 replicates here. 

ppe_imputed[,WT_ARID1A_log] <- ptImpute(ppe_imputed[,WT_WT_log], #control
                                        ppe_imputed[,WT_ARID1A_log], #ARID1A_Vermurafenib
                                        percent1 = .6, percent2 = .2, paired = FALSE)
ppe_imputed[,WT_MED12_log] <- ptImpute(ppe_imputed[,WT_WT_log], #control
                                       ppe_imputed[,WT_MED12_log], #MED12_Vermurafenib
                                       percent1 = .6, percent2 = .2, paired = FALSE)

ppe_imputed_scaled <- medianScaling(ppe_imputed, scale = F, assay = "imputed")

raw_quant<-na.omit(SummarizedExperiment::assay(ppe_imputed_scaled,"Quantification"))
raw_imp<-na.omit(SummarizedExperiment::assay(ppe_imputed_scaled,"imputed"))
scaled_imp<-na.omit(SummarizedExperiment::assay(ppe_imputed_scaled,"scaled"))

####BATCH CORRECTION####


# run ComBat as alternative to IRS
# NOTE: SL norm is probably better to do before applying the Batch corection

combat_input<-data.frame(scaled_imp)
combat_input[,grep("131C", colnames(combat_input))]<-NULL 

design_combat<- design
design_combat<-design_combat[match(colnames(combat_input), design_combat$code),]
batch<-unlist(design_combat$Plex[match(colnames(combat_input), design_combat$code)])
design_combat$Drug <- factor(design_combat$Drug, levels = c("Untreated", "Vermurafenib_1uM", "Trametinib_10nM", "vemurafenib+trametinib"))
design_combat$Genetic <- factor(design_combat$Genetic, levels = c("WT", "ARID1A_KO", "MED12_KO"))

mod <- model.matrix(~ Drug + Genetic, data = design_combat)
data_combat <- ComBat(dat = as.matrix(combat_input), batch = batch, mod = mod, par.prior = TRUE)
data_combat <- as.data.frame(data_combat)
par(mfrow = c(1, 1)) # any plotting in the ComBat call leaves plots as 2x2

# ComBat introduces some negative corrected counts that we need to fix
data_combat <-  data_combat[apply(data_combat, 1, function(x) all(x >= 0)), ] 




#plot diagnostic plots 
colour_scheme<-c("#F94144", "#F3722C", "#F8961E", "#F9844A", "#F9C74F", "#90BE6D", "#43AA8B")
sample_colour_scheme<-c("#001219","055363", "#005F73", "#0A9396", "#94D2BD", "#E9D8A6", "#EE9B00", "#CA6702", "#BB3E03",
                        "#AE2012", "#9B2226", "#C61E24", "#5727A1")

fontSize = 10



#raw
# Change density plot line colors by groups
raw_peptides_density<-raw_peptides
colnames(raw_peptides_density)<-str_replace_all(colnames(raw_peptides_density),pattern = "plex", replacement = "B")
colnames(raw_peptides_density)<-str_remove_all(colnames(raw_peptides_density),pattern = "X")
raw_molt_density<-reshape2::melt(log2(raw_peptides_density))
raw_molt_density$Batch<-unlist(map(str_split(raw_molt_density$variable, pattern = "\\."), last))

d1<-ggplot(raw_molt_density, aes(x=value, fill=Batch)) +
  geom_density(alpha=0.8)+scale_fill_manual(values=colour_scheme)+cowplot::theme_cowplot(font_size = fontSize) +
  xlim(5, 25) + 
  ggtitle("A.   Distribution of raw PEPTIDE abundances per batch") + 
  labs(y = "Density", x = "Log2(Peptide Abundance)")
# Change box plot colors by groups
b1<-ggplot(raw_molt_density, aes(y=value, x=as.character(variable), fill=Batch)) +
  geom_boxplot()+scale_fill_manual(values=colour_scheme)+cowplot::theme_cowplot(font_size = fontSize) +
  ylim(5, 25) + ggtitle("A.   Raw PEPTIDE abundances per sample") + 
  labs(y = "Log2(PEPTIDE Abundance)", x = "Samples") + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))


#SLN
sln_peptides_density<-sln_peptides
colnames(sln_peptides_density)<-str_replace_all(colnames(sln_peptides_density),pattern = "plex", replacement = "B")
colnames(sln_peptides_density)<-str_remove_all(colnames(sln_peptides_density),pattern = "X")
sln_molt_density<-reshape2::melt(log2(sln_peptides_density))
sln_molt_density$Batch<-unlist(map(str_split(sln_molt_density$variable, pattern = "\\."), last))

d2<-ggplot(sln_molt_density, aes(x=value, fill=Batch)) +
  geom_density(alpha=0.8)+scale_fill_manual(values=colour_scheme)+cowplot::theme_cowplot(font_size = fontSize) +
  xlim(5, 25) + ggtitle("B.   Distribution of sample loading normalised PEPTIDE abundances per batch") + 
  labs(y = "Density", x = "Log2(PEPTIDE Abundance)")
# Change box plot colors by groups
b2<-ggplot(sln_molt_density, aes(y=value, x=as.character(variable), fill=Batch)) +
  geom_boxplot()+scale_fill_manual(values=colour_scheme)+cowplot::theme_cowplot(font_size = fontSize) +
  ylim(5, 25) + ggtitle("B.   Sample loading normalised PEPTIDE abundances per sample") + 
  labs(y = "Log2(PEPTIDE Abundance)", x = "Samples") + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))


#scaled imputed 
combat_peptides_density<-data_combat
colnames(combat_peptides_density)<-str_replace_all(colnames(combat_peptides_density),pattern = "plex", replacement = "B")
colnames(combat_peptides_density)<-str_remove_all(colnames(combat_peptides_density),pattern = "X")
combat_molt_density<-reshape2::melt(log2(combat_peptides_density))
combat_molt_density$Batch<-unlist(map(str_split(combat_molt_density$variable, pattern = "\\."), last))

d3<-ggplot(combat_molt_density, aes(x=value, fill=Batch)) +
  geom_density(alpha=0.8)+scale_fill_manual(values=colour_scheme)+cowplot::theme_cowplot(font_size = fontSize) +
  xlim(5, 25) + ggtitle("C.   Distribution of ComBat normalised PEPTIDE abundances per batch") + 
  labs(y = "Density", x = "Log2(PEPTIDE Abundance)")
# Change box plot colors by groups
b3<-ggplot(combat_molt_density, aes(y=value, x=as.character(variable), fill=Batch)) +
  geom_boxplot()+scale_fill_manual(values=colour_scheme)+cowplot::theme_cowplot(font_size = fontSize) +
  ylim(5, 25) + ggtitle("C.   ComBat normalised PEPTIDE abundances per sample") + 
  labs(y = "Log2(PEPTIDE Abundance)", x = "Samples") + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

#CVs


# function computes CVs per sample
make_CVs <- function(df, design) {
  list_dfs<-list()
  exp_conditions<-design$total[match(colnames(df), design$code)]
  for (variable in names(table(exp_conditions)[table(exp_conditions) != 1])) {
    cond_df<-data.frame(df[,exp_conditions==variable])
    cond_df$ave <- rowMeans(cond_df)
    cond_df$sd <- apply(data.frame(cond_df[,colnames(cond_df)!="ave"]), 1, sd)
    cond_df$cv <- 100 * cond_df$sd / cond_df$ave
    list_dfs[[variable]]<-cond_df
  }
  # save results in 3 data frames and put into a list
  ave_df <- data.frame(lapply(list_dfs, function(x){x$ave}))
  sd_df <- data.frame(lapply(list_dfs, function(x){x$sd}))
  cv_df <- data.frame(lapply(list_dfs, function(x){x$cv}))
  
  cv_molt<-reshape2::melt(cv_df)
  colnames(cv_molt)<-c("Sample", "CV")
  return(cv_molt)
}

# get CVs and averages
list_raw <- make_CVs(raw_peptides, design = design)
list_sln <- make_CVs(sln_peptides, design = design)
list_combat <- make_CVs(data_combat, design = design)

median(list_combat$CV)

# Change box plot colors by groups
cv1<-ggplot(list_raw, aes(y=CV, x=Sample, fill=Sample)) +
  geom_boxplot()+cowplot::theme_cowplot(font_size = fontSize) +
  ylim(0, 200) + ggtitle("A.   RAW Coefficient Of Variation for each PEPTIDE per sample (CVs)") + 
  labs(y = "CV (%)", x = "Samples") + 
  theme(axis.text.x = element_blank())+
  scale_fill_manual(values=sample_colour_scheme)
cv2<-ggplot(list_sln, aes(y=CV, x=Sample, fill=Sample)) +
  geom_boxplot()+cowplot::theme_cowplot(font_size = fontSize) +
  ylim(0, 200) + ggtitle("B.   SL Coefficient Of Variation for each PEPTIDE per sample (CVs)") + 
  labs(y = "CV (%)", x = "Samples") + 
  theme(axis.text.x = element_blank())+
  scale_fill_manual(values=sample_colour_scheme)
cv3<-ggplot(list_combat, aes(y=CV, x=Sample, fill=Sample)) +
  geom_boxplot()+cowplot::theme_cowplot(font_size = fontSize) +
  ylim(0, 200) + ggtitle("C.   ComBat Coefficient Of Variation for each PEPTIDE per sample (CVs)") + 
  labs(y = "CV (%)", x = "Samples") + 
  theme(axis.text.x = element_blank())+
  scale_fill_manual(values=sample_colour_scheme)




#prep phosR data
imput_quant<-function(x){
  numberMissing <- colSums(!is.na(x))
  numberRows<-dim(x)[1]
  out<-data.frame(numberMissing/numberRows)
  rownames(out)<-str_replace_all(rownames(out),pattern = "plex", replacement = "B")
  rownames(out)<-str_remove_all(rownames(out),pattern = "X")
  out$name<-rownames(out)
  out$Batch<-unlist(map(str_split(rownames(out), pattern = "\\."), last))
  rownames(out)<-NULL
  colnames(out)<-c("ProportionMissing", "sampleName", "Batch")
  return(out)
}

pre_filter = imput_quant(SummarizedExperiment::assay(ppe,"Quantification"))
# Change the colors manually
i1 <- ggplot(data=pre_filter, aes(x=sampleName, y=ProportionMissing, fill=Batch)) +
  geom_bar(stat="identity", color="black", position=position_dodge())+
  theme_minimal()+scale_fill_manual(values=colour_scheme)+cowplot::theme_cowplot(font_size = fontSize) +
  ylim(0, 1) + ggtitle("A.   Quantification of samples in raw PEPTIDE data") + 
  labs(y = "Proportion missing", x = "Samples") + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + 
  geom_hline(yintercept=median(pre_filter$ProportionMissing), linetype="dashed", color = "red", size = 1.2)
post_filter = imput_quant(SummarizedExperiment::assay(ppe_filtered,"Quantification"))
# Change the colors manually
i2 <- ggplot(data=post_filter, aes(x=sampleName, y=ProportionMissing, fill=Batch)) +
  geom_bar(stat="identity", color="black", position=position_dodge())+
  theme_minimal()+scale_fill_manual(values=colour_scheme)+cowplot::theme_cowplot(font_size = fontSize) +
  ylim(0, 1) + ggtitle("B.   Quantification of samples in filtered PEPTIDE data") + 
  labs(y = "Proportion missing", x = "Samples") + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + 
  geom_hline(yintercept=median(post_filter$ProportionMissing), linetype="dashed", color = "red", size = 1.2)
post_impute = imput_quant(SummarizedExperiment::assay(ppe_imputed_tmp,"imputed"))
# Change the colors manually
i3 <- ggplot(data=post_impute, aes(x=sampleName, y=ProportionMissing, fill=Batch)) +
  geom_bar(stat="identity", color="black", position=position_dodge())+
  theme_minimal()+scale_fill_manual(values=colour_scheme)+cowplot::theme_cowplot(font_size = fontSize) +
  ylim(0, 1) + ggtitle("C.   Quantification of samples in site/condition specific imputed PEPTIDE data") + 
  labs(y = "Proportion missing", x = "Samples") + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + 
  geom_hline(yintercept=median(post_impute$ProportionMissing), linetype="dashed", color = "red", size = 1.2)
post_impute_t = imput_quant(SummarizedExperiment::assay(ppe_imputed_scaled,"scaled"))
# Change the colors manually
i4 <- ggplot(data=post_impute_t, aes(x=sampleName, y=ProportionMissing, fill=Batch)) +
  geom_bar(stat="identity", color="black", position=position_dodge())+
  theme_minimal()+scale_fill_manual(values=colour_scheme)+cowplot::theme_cowplot(font_size = fontSize) +
  ylim(0, 1) + ggtitle("D.Quantification of samples in paired-tail imputed PEPTIDE data") + 
  labs(y = "Proportion missing", x = "Samples") + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + 
  geom_hline(yintercept=median(post_impute_t$ProportionMissing), linetype="dashed", color = "red", size = 1.2)


#PCAs showing Batch correction 


# see how things cluster now that we have nice boxplots and density plots
MDS_raw<-plotMDS(log2(raw_peptides), plot = F)
library(ComplexHeatmap)
binned_col<-colour_scheme
names(binned_col)<-unique(as.character(design$Plex))
PCA_out_raw<-data.frame(Sample=names(MDS_raw$x),
                        x=MDS_raw$x,
                        y=MDS_raw$y,
                        Batch=as.character(design$Plex[match(names(MDS_raw$x), design$code)]),
                        Drug=design$Drug[match(names(MDS_raw$x), design$code)],
                        Genetic=design$Genetic[match(names(MDS_raw$x), design$code)],
                        All_conditions=design$total[match(names(MDS_raw$x), design$code)])
ha = HeatmapAnnotation(Drug = design$Drug[match(colnames(MDS_raw$distance.matrix), design$code)], 
                       Genetic = design$Genetic[match(colnames(MDS_raw$distance.matrix), design$code)],
                       Batch = as.character(design$Plex[match(colnames(MDS_raw$distance.matrix), design$code)]),
                       col = list(Batch = binned_col))
colnames(MDS_raw$distance.matrix)<-NULL
rownames(MDS_raw$distance.matrix)<-NULL
heat_raw<-Heatmap(MDS_raw$distance.matrix, rect_gp = gpar(type = "none"), column_dend_side = "bottom", top_annotation = ha,
        cell_fun = function(j, i, x, y, w, h, fill) {
          if(as.numeric(x) <= 1 - as.numeric(y)) {
            grid.rect(x, y, w, h, gp = gpar(fill = fill, col = fill))
          }
        })
PCA_out_raw<-na.omit(PCA_out_raw)
# Change point shapes, colors and sizes
p1<-ggplot(PCA_out_raw, aes(x=x, y=y, shape=Drug, color=Batch)) +
  geom_point()+scale_color_manual(values=colour_scheme)+cowplot::theme_cowplot(font_size = fontSize) + 
  ggtitle("A.   Multidimensional scaling plot of distances between raw PEPTIDE samples")


#sample loading normalisation 
MDS_sln<-plotMDS(log2(sln_peptides), plot = F)
PCA_out_sln<-data.frame(Sample=names(MDS_sln$x),
                        x=MDS_sln$x,
                        y=MDS_sln$y,
                        Batch=as.character(design$Plex[match(names(MDS_sln$x), design$code)]),
                        Drug=design$Drug[match(names(MDS_sln$x), design$code)],
                        Genetic=design$Genetic[match(names(MDS_sln$x), design$code)],
                        All_conditions=design$total[match(names(MDS_sln$x), design$code)])
ha = HeatmapAnnotation(Drug = design$Drug[match(colnames(MDS_sln$distance.matrix), design$code)], 
                       Genetic = design$Genetic[match(colnames(MDS_sln$distance.matrix), design$code)],
                       Batch = as.character(design$Plex[match(colnames(MDS_sln$distance.matrix), design$code)]),
                       col = list(Batch = binned_col))
colnames(MDS_sln$distance.matrix)<-NULL
rownames(MDS_sln$distance.matrix)<-NULL
heat_sln<-Heatmap(MDS_sln$distance.matrix, rect_gp = gpar(type = "none"), column_dend_side = "bottom", top_annotation = ha,
                  cell_fun = function(j, i, x, y, w, h, fill) {
                    if(as.numeric(x) <= 1 - as.numeric(y)) {
                      grid.rect(x, y, w, h, gp = gpar(fill = fill, col = fill))
                    }
                  })
PCA_out_sln<-na.omit(PCA_out_sln)
# Change point shapes, colors and sizes
p2<-ggplot(PCA_out_sln, aes(x=x, y=y, shape=Drug, color=Batch)) +
  geom_point()+scale_color_manual(values=colour_scheme)+cowplot::theme_cowplot(font_size = fontSize) + 
  ggtitle("B.   Multidimensional scaling plot of distances between sample loading normalised PEPTIDE samples")

#imputation 
MDS_impute<-plotMDS(log2(SummarizedExperiment::assay(ppe_imputed_scaled,"scaled")), plot = F)
PCA_out_impute<-data.frame(Sample=names(MDS_impute$x),
                        x=MDS_impute$x,
                        y=MDS_impute$y,
                        Batch=as.character(design$Plex[match(names(MDS_impute$x), design$code)]),
                        Drug=design$Drug[match(names(MDS_impute$x), design$code)],
                        Genetic=design$Genetic[match(names(MDS_impute$x), design$code)],
                        All_conditions=design$total[match(names(MDS_impute$x), design$code)])

PCA_out_impute<-na.omit(PCA_out_impute)
# Change point shapes, colors and sizes
p3<-ggplot(PCA_out_impute, aes(x=x, y=y, shape=Drug, color=Batch)) +
  geom_point()+scale_color_manual(values=colour_scheme)+cowplot::theme_cowplot(font_size = fontSize) + 
  ggtitle("Multidimensional scaling plot of distances between imputed PEPTIDE samples")


 #ComBAT Batch correction
MDS_combat<-plotMDS(log2(data_combat), plot = F)
PCA_out_combat<-data.frame(Sample=names(MDS_combat$x),
                           x=MDS_combat$x,
                           y=MDS_combat$y,
                           Batch=as.character(design$Plex[match(names(MDS_combat$x), design$code)]),
                           Drug=design$Drug[match(names(MDS_combat$x), design$code)],
                           Genetic=design$Genetic[match(names(MDS_combat$x), design$code)],
                           All_conditions=design$total[match(names(MDS_combat$x), design$code)])

PCA_out_combat<-na.omit(PCA_out_combat)
# Change point shapes, colors and sizes
p4<-ggplot(PCA_out_combat, aes(x=x, y=y, shape=Drug, color=Batch)) +
  geom_point()+scale_color_manual(values=colour_scheme)+cowplot::theme_cowplot(font_size = fontSize)+ 
  ggtitle("C.   Multidimensional scaling plot of distances between ComBat corrected PEPTIDE samples")


#WRITE PDFS
pdf(file = "~/Thesis/figures/AppendixB/phospho_preprocessing.pdf",width=8.27,height=10)
ggpubr::ggarrange(d1, d2, d3, nrow = 3) #density plots per batch
ggpubr::ggarrange(b1, b2, b3, nrow = 3) #boxplots per sample 
ggpubr::ggarrange(i1, i2, i3, i4, nrow = 4) #imputation plots 
ggpubr::ggarrange(cv1, cv2, cv3, nrow = 3) #CVs 
dev.off()

pdf(file = "~/Thesis/figures/AppendixB/phospho_PCA.pdf",width=8,height=8)
p1
p2
p4
dev.off()
