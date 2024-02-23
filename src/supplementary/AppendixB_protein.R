# Analysis of IOVS mouse lens data (Supplemental Table S01):
# Khan, Shahid Y., et al. "Proteome Profiling of Developing Murine Lens Through Mass Spectrometry."
# Investigative Ophthalmology & Visual Science 59.1 (2018): 100-107.

# load libraries (this gets ggplot2 and related libraries)
library(tidyverse)
library(reshape2)
# these are from Bioconductor
library(limma) 
library(edgeR) 
library(sva)
library(psych)


#This looks at the proteomics from ALL the batches, since the drug treatment corresponds too well with the batches, are batch correction procedures will filter this out
#so we can only really look at the differences between the KOs. 

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

par(mfrow = c(3, 2))

path_data<-"/home/charlie/phd/MelanomaProject//data/phosphoproteomics/main/swiss-prot/total/" #SWISS-plot
setwd(path_data)
# read the Supplemental 01 file (saved as a CSV export from XLSX file)
data_start <- read_csv("Proteins_grouped_abundances.csv")
design<-read.csv(file = "~/phd/MelanomaProject/data/phosphoproteomics/main/design.csv",header = T)
design$code<-paste0(design$TMT.labels,".pl", design$Plex)
design$Drug <- str_replace_all(design$Drug, pattern = "WT", replacement = "Untreated")
design$total<-paste0(design$Drug,"__", design$Genetic)
# save the annotations (gene symbol and protein accession) and remove from data frame
annotate_df <- data_start[1:20] 
data_raw <- as.data.frame(data_start[grep("Abundances.Grouped", colnames(data_start))])
row_labels<-annotate_df$Accession
rownames(data_raw) <- row_labels

# filter out proteins not seen in all three runs
data_no_na <- na.omit(data_raw)
data_no_na[,grep("Count", colnames(data_no_na))]<-NULL
data_no_na[,grep("CV", colnames(data_no_na))]<-NULL
data_no_na[,grep("131N", colnames(data_no_na))]<-NULL #this is an empty channel that Fransezka used because they had one spare. 
data_no_na[,grep("131C", colnames(data_no_na))]<-NULL 

# fix the column headers
col_headers <- colnames(data_no_na) 
col_headers <- str_replace(col_headers, "Abundances.Grouped.", "") 
colnames(data_no_na) <- col_headers

data_in<-data_no_na
data_sl<-normSLN(data_in) #PERFECT! the function does work (minor panic there!!)

# run ComBat as alternative to IRS
# NOTE: SL norm is probably better to do before applying the batch corection

design_combat<- design
design_combat$code<-str_remove(string = design$code, pattern = "X")
design_combat<-design_combat[match(colnames(data_sl), design_combat$code),]
batch<-unlist(design_combat$Plex[match(colnames(data_sl), design_combat$code)])
design_combat$Drug <- factor(design_combat$Drug, levels = c("Untreated", "Vermurafenib_1uM", "Trametinib_10nM", "vemurafenib+trametinib"))
design_combat$Genetic <- factor(design_combat$Genetic, levels = c("WT", "ARID1A_KO", "MED12_KO"))

mod <- model.matrix(~ Drug + Genetic, data = design_combat)
data_combat <- ComBat(dat = as.matrix(data_sl), batch = batch, mod = mod, par.prior = TRUE)
data_combat <- as.data.frame(data_combat)


# apply TMM normalization to the ComBat-corrected data
# combat_tmm <- calcNormFactors(data_combat)
# data_combat <- sweep(data_combat, 2, combat_tmm, FUN = "/")

# ComBat introduces some negative corrected counts that we need to fix
data_combat <-  data_combat[apply(data_combat, 1, function(x) all(x >= 0)), ] 

#3D plot
dat_pc<-t(data_combat)
table(design$total[match(rownames(dat_pc), design_combat$code)])

pca <- prcomp(dat_pc, scale.=TRUE )
pca3d::pca3d(pca, group = design_combat$Drug[match(rownames(dat_pc), design_combat$code)], bg = "white",show.centroids = T) #unlist(map(str_split(rownames(dat_pc), "\\."), last))
pca3d::snapshotPCA3d("~/Desktop/Genetic.png")


#PLOTS

colour_scheme<-c("#F94144", "#F3722C", "#F8961E", "#F9844A", "#F9C74F", "#90BE6D", "#43AA8B")
sample_colour_scheme<-c("#001219","055363", "#005F73", "#0A9396", "#94D2BD",
                        "#E9D8A6", "#EE9B00", "#CA6702", "#BB3E03",
                        "#AE2012", "#9B2226", "#C61E24", "#5727A1")
fontSize = 10


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
list_raw <- make_CVs(data_in, design = design)
list_sln <- make_CVs(data_sl, design = design)
list_combat <- make_CVs(data_combat, design = design)

median(list_combat$CV)

# Change box plot colors by groups
cv1<-ggplot(list_raw, aes(y=CV, x=Sample, fill=Sample))+
  geom_boxplot()+cowplot::theme_cowplot(font_size = fontSize) +
  ylim(0, 200) + ggtitle("A.   RAW Coefficient Of Variation for each PROTEIN per sample (CVs)") + 
  labs(y = "CV (%)", x = "Samples") + 
  theme(axis.text.x = element_blank())+
  scale_fill_manual(values=sample_colour_scheme)
cv2<-ggplot(list_sln, aes(y=CV, x=Sample, fill=Sample))+
  geom_boxplot()+cowplot::theme_cowplot(font_size = fontSize) +
  ylim(0, 200) + ggtitle("A.   SLN Coefficient Of Variation for each PROTEIN per sample (CVs)") + 
  labs(y = "CV (%)", x = "Samples") + 
  theme(axis.text.x = element_blank())+
  scale_fill_manual(values=sample_colour_scheme)
cv3<-ggplot(list_combat, aes(y=CV, x=Sample, fill=Sample))+
  geom_boxplot()+cowplot::theme_cowplot(font_size = fontSize) +
  ylim(0, 200) + ggtitle("A.   ComBat Coefficient Of Variation for each PROTEIN per sample (CVs)") + 
  labs(y = "CV (%)", x = "Samples") + 
  theme(axis.text.x = element_blank())+
  scale_fill_manual(values=sample_colour_scheme)



#raw
# Change density plot line colors by groups
xlim=c(10,25)

raw_peptides_density<-data_in
colnames(raw_peptides_density)<-str_replace_all(colnames(raw_peptides_density),pattern = "pl", replacement = "B")
colnames(raw_peptides_density)<-str_remove_all(colnames(raw_peptides_density),pattern = "X")
raw_molt_density<-reshape2::melt(log2(raw_peptides_density))
raw_molt_density$Batch<-unlist(map(str_split(raw_molt_density$variable, pattern = "\\."), last))

d1<-ggplot(raw_molt_density, aes(x=value, fill=Batch)) +
  geom_density(alpha=0.8)+scale_fill_manual(values=colour_scheme)+cowplot::theme_cowplot(font_size = fontSize) +
  xlim(xlim[1], xlim[2]) + 
  ggtitle("A.   Distribution of raw PROTEIN abundances per batch") + 
  labs(y = "Density", x = "Log2(Peptide Abundance)")
# Change box plot colors by groups
b1<-ggplot(raw_molt_density, aes(y=value, x=as.character(variable), fill=Batch)) +
  geom_boxplot()+scale_fill_manual(values=colour_scheme)+cowplot::theme_cowplot(font_size = fontSize) +
  ylim(xlim[1], xlim[2]) + 
  ggtitle("A.   Raw PROTEIN abundances per sample") + 
  labs(y = "Log2(PROTEIN Abundance)", x = "Samples") + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))


#SLN
sln_peptides_density<-data_sl
colnames(sln_peptides_density)<-str_replace_all(colnames(sln_peptides_density),pattern = "pl", replacement = "B")
colnames(sln_peptides_density)<-str_remove_all(colnames(sln_peptides_density),pattern = "X")
sln_molt_density<-reshape2::melt(log2(sln_peptides_density))
sln_molt_density$Batch<-unlist(map(str_split(sln_molt_density$variable, pattern = "\\."), last))

d2<-ggplot(sln_molt_density, aes(x=value, fill=Batch)) +
  geom_density(alpha=0.8)+scale_fill_manual(values=colour_scheme)+cowplot::theme_cowplot(font_size = fontSize) +
  xlim(xlim[1], xlim[2]) + 
  ggtitle("B.   Distribution of sample loading normalised PROTEIN abundances per batch") + 
  labs(y = "Density", x = "Log2(PROTEIN Abundance)")
# Change box plot colors by groups
b2<-ggplot(sln_molt_density, aes(y=value, x=as.character(variable), fill=Batch)) +
  geom_boxplot()+scale_fill_manual(values=colour_scheme)+cowplot::theme_cowplot(font_size = fontSize) +
  ylim(xlim[1], xlim[2]) + 
  ggtitle("B.   Sample loading normalised PROTEIN abundances per sample") + 
  labs(y = "Log2(PROTEIN Abundance)", x = "Samples") + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))


#combat corrected  
combat_peptides_density<-data_combat
colnames(combat_peptides_density)<-str_replace_all(colnames(combat_peptides_density),pattern = "pl", replacement = "B")
colnames(combat_peptides_density)<-str_remove_all(colnames(combat_peptides_density),pattern = "X")
combat_molt_density<-reshape2::melt(log2(combat_peptides_density))
combat_molt_density$Batch<-unlist(map(str_split(combat_molt_density$variable, pattern = "\\."), last))

d3<-ggplot(combat_molt_density, aes(x=value, fill=Batch)) +
  geom_density(alpha=0.8)+scale_fill_manual(values=colour_scheme)+cowplot::theme_cowplot(font_size = fontSize) +
  xlim(xlim[1], xlim[2]) + 
  ggtitle("C.   Distribution of ComBat normalised PROTEIN abundances per batch") + 
  labs(y = "Density", x = "Log2(PROTEIN Abundance)")
# Change box plot colors by groups
b3<-ggplot(combat_molt_density, aes(y=value, x=as.character(variable), fill=Batch)) +
  geom_boxplot()+scale_fill_manual(values=colour_scheme)+cowplot::theme_cowplot(font_size = fontSize) +
  ylim(xlim[1], xlim[2]) + 
  ggtitle("C.   ComBat normalised PROTEIN abundances per sample") + 
  labs(y = "Log2(PROTEIN Abundance)", x = "Samples") + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))



#PCAs showing Batch correction 


# see how things cluster now that we have nice boxplots and density plots
# function computes CVs per sample

#raw
MDS_sln<-plotMDS(log2(data_in), plot = F)
PCA_out_sln<-data.frame(Sample=names(MDS_sln$x),
                        x=MDS_sln$x,
                        y=MDS_sln$y,
                        Batch=as.character(design$Plex[match(names(MDS_sln$x), design$code)]),
                        Drug=design$Drug[match(names(MDS_sln$x), design$code)],
                        Genetic=design$Genetic[match(names(MDS_sln$x), design$code)],
                        All_conditions=design$total[match(names(MDS_sln$x), design$code)])
PCA_out_sln<-na.omit(PCA_out_sln)
# Change point shapes, colors and sizes
p1<-ggplot(PCA_out_sln, aes(x=x, y=y, shape=Drug, color=Batch)) +
  geom_point()+scale_color_manual(values=colour_scheme)+cowplot::theme_cowplot(font_size = fontSize) + 
  ggtitle("A.   Multidimensional scaling plot of distances between raw PROTEIN samples")


#sample loading normalisation 
MDS_sln<-plotMDS(log2(data_sl), plot = F)
PCA_out_sln<-data.frame(Sample=names(MDS_sln$x),
                        x=MDS_sln$x,
                        y=MDS_sln$y,
                        Batch=as.character(design$Plex[match(names(MDS_sln$x), design$code)]),
                        Drug=design$Drug[match(names(MDS_sln$x), design$code)],
                        Genetic=design$Genetic[match(names(MDS_sln$x), design$code)],
                        All_conditions=design$total[match(names(MDS_sln$x), design$code)])
PCA_out_sln<-na.omit(PCA_out_sln)
# Change point shapes, colors and sizes
p2<-ggplot(PCA_out_sln, aes(x=x, y=y, shape=Drug, color=Batch)) +
  geom_point()+scale_color_manual(values=colour_scheme)+cowplot::theme_cowplot(font_size = fontSize) + 
  ggtitle("B.   Multidimensional scaling plot of distances between sample loading normalised PROTEIN samples")

#ComBAT
MDS_sln<-plotMDS(log2(data_combat), plot = F)
PCA_out_sln<-data.frame(Sample=names(MDS_sln$x),
                        x=MDS_sln$x,
                        y=MDS_sln$y,
                        Batch=as.character(design$Plex[match(names(MDS_sln$x), design$code)]),
                        Drug=design$Drug[match(names(MDS_sln$x), design$code)],
                        Genetic=design$Genetic[match(names(MDS_sln$x), design$code)],
                        All_conditions=design$total[match(names(MDS_sln$x), design$code)])
PCA_out_sln<-na.omit(PCA_out_sln)
# Change point shapes, colors and sizes
p3<-ggplot(PCA_out_sln, aes(x=x, y=y, shape=Drug, color=Batch)) +
  geom_point()+scale_color_manual(values=colour_scheme)+cowplot::theme_cowplot(font_size = fontSize) + 
  ggtitle("C.   Multidimensional scaling plot of distances between ComBat-corrected PROTEIN samples")

#Protein abundance and RNA transcript correlation. 

library(edgeR)
library(EnsDb.Hsapiens.v86)
library(ggpubr)
library(rstatix)

data_protein<-data_combat
All_counts<- read.csv("~/MelanomaProject/Sumana_RNAseq_analysis/geneCounts_fixed.csv", header=T, fill=T)
Sample_info<- read.csv("~/MelanomaProject/Sumana_RNAseq_analysis/Study_design.csv", header=T)
conv<-read.delim(file = "../../../rna_prot_design.txt", sep = ",", header = F)

Sample_info$Drug <- unlist(map(str_split(Sample_info$Sample_name, pattern = "__"), 1))
Sample_info$Genetic <- unlist(map(str_split(Sample_info$Sample_name, pattern = "__"), 2))

#Let's normalise the counts now:
mycounts <- All_counts
rownames(mycounts) <- All_counts$ENSEMBL_ID
mycounts<- mycounts[,2:37]
isexpr <- rowSums(cpm(mycounts)>2) >= 1
mycounts <- mycounts[isexpr,]
experiment_design<-Sample_info
rownames(experiment_design) <- experiment_design$Study_ID
group<-factor(experiment_design$Type)
design_mat <- model.matrix(~0+group)
colnames(design_mat)<- gsub("group","",colnames(design_mat))
nf <- calcNormFactors(mycounts)
y <- voom(mycounts,design_mat,lib.size=colSums(mycounts)*nf, plot = T)
counts.voom <- y$E
# save normalised expression data into output dir
Normalised_data_to_print<- as.data.frame(counts.voom)
data_rna<- melt(data.matrix(mycounts))

genename_df<-AnnotationDbi::select(EnsDb.Hsapiens.v86, keys = as.character(data_rna$Var1), keytype = "GENEID", columns = "UNIPROTID")
data_rna$UNIPROTID <- genename_df$UNIPROTID[match(as.character(data_rna$Var1) , genename_df$GENEID)]

Sample_info$Sample_name<-conv$V1[match(Sample_info$Sample_name, conv$V2)]
Sample_info$Drug<-unlist(map(str_split(Sample_info$Sample_name, pattern = "__"),1))
Sample_info$Genetic<-unlist(map(str_split(Sample_info$Sample_name, pattern = "__"),2))
Sample_info$Drug <- str_replace_all(Sample_info$Drug, pattern = "WT", replacement = "Untreated")
Sample_info$total<-paste0(Sample_info$Drug,"__", Sample_info$Genetic)
#select only WTs
cor_list<-list()
for (biological_condition in unique(Sample_info$total)) {
  WT_protein<-data_protein[colnames(data_protein) %in% design$code[design$total == biological_condition]]
  WT_protein<-reshape2::melt(data.matrix(WT_protein))
  
  WT_rna<-data_rna[data_rna$Var2 %in% Sample_info$Study_ID[Sample_info$total == biological_condition],]
  #aggregate counts of replicates 
  WT_rna <- WT_rna %>%
    group_by(UNIPROTID) %>%
    summarise_at(vars(value), funs(mean(., na.rm=TRUE)))
  WT_protein <- WT_protein %>%
    group_by(Var1) %>%
    summarise_at(vars(value), funs(mean(., na.rm=TRUE)))
  #merge
  prot_rna<-merge(x = WT_protein, y = WT_rna, by.x = "Var1", by.y = "UNIPROTID")
  cor_list[[biological_condition]]<-data.frame(prot_rna, biological_condition=biological_condition)
}
library(tidyverse)
prot_rna_toPlot<-bind_rows(cor_list)
#plot all proteins/protein coding genes 
c1<-ggplot(prot_rna_toPlot, aes(x=log2(value.x), y=log2(value.y))) + geom_point(size = 0.1) + geom_density_2d(bins=20)+stat_cor(method="pearson")+theme_classic()+ 
  facet_wrap(~ biological_condition) + ggtitle("Expression v Protein abundance = all genes") +
  xlab("Log2(Protein Abundance)") + ylab("Log2(Gene count)") + cowplot::theme_cowplot(font_size = fontSize)


#WRITE PDFS
pdf(file = "~/Thesis/figures/AppendixB/PROTEIN_preprocessing.pdf",width=8.27,height=10)
ggpubr::ggarrange(d1, d2, d3, nrow = 3) #density plots per batch
ggpubr::ggarrange(b1, b2, b3, nrow = 3) #boxplots per sample 
ggpubr::ggarrange(cv1, cv2, cv3, nrow = 3) #CVs 
c1
dev.off()

pdf(file = "~/Thesis/figures/AppendixB/PROTEIN_MDS.pdf",width=8,height=8)
p1
p2
p3
dev.off()
