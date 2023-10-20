
set.seed(1234)

path_data<-"~/Desktop/Melanoma_Resistance/" #SWISS-plot
setwd(path_data)

####LOAD FUNCTIONS####

source("./src/functions/SLN_normalisation.R")
source("./src/functions/plot_3D.R")

# load libraries (this gets ggplot2 and related libraries)
library(tidyverse)
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


####LOADING DATA####

# read the Supplemental 01 file (saved as a CSV export from XLSX file)
data_start <- read_csv("./data/proteomic/raw/Proteins_grouped_abundances.csv")
design<-read.csv(file = "./data/proteomic/design.csv",header = T)
design$code<-paste0(design$TMT.labels,".pl", design$Plex)
design$total<-paste0(design$Drug,"__", design$Genetic)
# save the annotations (gene symbol and protein accession) and remove from data frame
annotate_df <- data_start[1:20] 
data_raw <- as.data.frame(data_start[grep("Abundances.Grouped", colnames(data_start))])
row_labels<-annotate_df$Accession
rownames(data_raw) <- row_labels

# filter out proteins not seen in all three runs
data_no_na <- data_raw
data_no_na[,grep("Count", colnames(data_no_na))]<-NULL
data_no_na[,grep("CV", colnames(data_no_na))]<-NULL
data_no_na[,grep("131N", colnames(data_no_na))]<-NULL #this is an empty channel that Fransezka used because they had one spare. 
data_no_na[,grep("131C", colnames(data_no_na))]<-NULL 

# fix the column headers
col_headers <- colnames(data_no_na) 
col_headers <- str_replace(col_headers, "Abundances.Grouped.", "") 
colnames(data_no_na) <- col_headers
combided_proteins<-data_no_na

# read the design matrix
design<-read.csv(file = "./data/proteomic/design.csv",header = T)
design$code<-paste0(design$TMT.labels,".pl", design$Plex)
design$total<-paste0(design$Drug,"__", design$Genetic)
grps <- design$total[match(colnames(combided_proteins),design$code)]

####IMPUTATION####

combided_proteins<-normSLN(combided_proteins)
prot_data<-combided_proteins
#extract phospho-site information
site = readr::parse_number(unlist(map(str_split(rownames(prot_data), pattern = "__"), last)))
uniprot_id = rownames(prot_data)
genename_df<-AnnotationDbi::select(EnsDb.Hsapiens.v86, keys = uniprot_id, keytype = "UNIPROTID", columns = "GENENAME")
geneSymbol <- genename_df$GENENAME

counts <- table(design$total[match(colnames(prot_data), design$code)])
ppe <- PhosphoExperiment(assays = list(Quantification = as.matrix(prot_data)), 
                         GeneSymbol = geneSymbol)

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
#look at our data

p1 = plotQC(SummarizedExperiment::assay(ppe,"Quantification"), 
            labels=colnames(ppe), 
            panel = "quantify", grps = grps)
p2 = plotQC(SummarizedExperiment::assay(ppe_filtered,"Quantification"), 
            labels=colnames(ppe_filtered), 
            panel = "quantify", grps = grps)
p3 = plotQC(SummarizedExperiment::assay(ppe_imputed_scaled,"scaled"), 
            labels=colnames(ppe_imputed_scaled), panel = "quantify", grps = grps)
ggpubr::ggarrange(p1, p2, p3, nrow = 1)

raw_quant<-na.omit(SummarizedExperiment::assay(ppe_imputed_scaled,"Quantification"))
dim(raw_quant)

raw_imp<-na.omit(SummarizedExperiment::assay(ppe_imputed_scaled,"imputed"))
dim(raw_imp)

scaled_imp<-na.omit(SummarizedExperiment::assay(ppe_imputed_scaled,"scaled"))
dim(scaled_imp)

####BATCH CORRECTION####


# run ComBat as alternative to IRS
# NOTE: SL norm is probably better to do before applying the batch corection

combat_input<-data.frame(raw_imp)
combat_input[,grep("131C", colnames(combat_input))]<-NULL 

design_combat<- design
design_combat<-design_combat[match(colnames(combat_input), paste0("X", design_combat$code)),]
batch<-unlist(design_combat$Plex[match(colnames(combat_input),  paste0("X", design_combat$code))])
design_combat$Drug <- factor(design_combat$Drug, levels = c("WT", "Vermurafenib_1uM", "Trametinib_10nM", "vemurafenib+trametinib"))
design_combat$Genetic <- factor(design_combat$Genetic, levels = c("WT", "ARID1A_KO", "MED12_KO"))

mod <- model.matrix(~ Drug + Genetic, data = design_combat)
data_combat <- ComBat(dat = as.matrix(combat_input), batch = batch, mod = mod, par.prior = TRUE)
data_combat <- as.data.frame(data_combat)
par(mfrow = c(1, 1)) # any plotting in the ComBat call leaves plots as 2x2

# ComBat introduces some negative corrected counts that we need to fix
data_combat <-  data_combat[apply(data_combat, 1, function(x) all(x >= 0)), ] 


# usual box plot views
boxplot(log2(data_combat), notch = TRUE, col = rep(c("red", "green", "blue"), each = 9), 
        main = "ComBat batch correction of SL data\nExp1 (red), Exp2 (green), Exp3 (blue)",
        xlab = 'TMT Sample', ylab = 'log2 of Intensity')

# can also look at density plots (like a distribution histogram)
plotDensities(log2(data_combat), col = rep(c("red", "green", "blue"), 6), main = "ComBat data")
#3D plot
dat_pc<-t(data_combat)
table(design$total[match(rownames(dat_pc), design_combat$code)])

pca <- prcomp(dat_pc, scale.=TRUE )

#write end result 
write.csv(x = data_combat,file = "./data/proteomic/processed/protein__norm_combat_imp.csv")

