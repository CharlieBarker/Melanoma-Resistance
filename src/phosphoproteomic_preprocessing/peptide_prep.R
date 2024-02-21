library(tidyverse)
library(limma) 
library(edgeR) 
library(sva)
library(psych)
suppressPackageStartupMessages({
  library(PhosR)
})
library(EnsDb.Hsapiens.v86)

set.seed(1234)

path_data<-"~/Desktop/Melanoma_Resistance/" #SWISS-plot
setwd(path_data)

####LOAD FUNCTIONS####

source("./src/functions/SLN_normalisation.R")
####LOADING DATA####


combined_peptides<-read.csv(file = "./data/proteomic/raw/peptide_total.csv")
rownames(combined_peptides)<-combined_peptides$X
combined_peptides$X<-NULL
# read the design matrix
design<-read.csv(file = "./data/proteomic/design.csv",header = T)
grps <- design$total[match(colnames(combined_peptides),design$code)]

#remove MED12
combined_peptides<-combined_peptides[,!grepl(grps, pattern = "MED12")]
design<-design[!grepl(design$total, pattern = "MED12"),]
grps<-grps[!grepl(grps, pattern = "MED12")]

####IMPUTATION####

combined_peptides<-normSLN(combined_peptides)

phos_data<-combined_peptides
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


ppe_filtered <- selectGrps(ppe, grps, 0.5, n=1) 

#scImpute 

ppe_imputed_tmp <- scImpute(ppe_filtered, 0.6, grps)[,colnames(ppe_filtered)]

# Paired tail-based imputation
ppe_imputed <- ppe_imputed_tmp

ppe_imputed_scaled <- medianScaling(ppe_imputed, scale = F, assay = "imputed")
scaled_imp<-na.omit(SummarizedExperiment::assay(ppe_imputed_scaled,"scaled"))

####BATCH CORRECTION####


# run ComBat as alternative to IRS
# NOTE: SL norm is probably better to do before applying the batch corection

combat_input<-data.frame(scaled_imp)
combat_input[,grep("131C", colnames(combat_input))]<-NULL 

design_combat<- design
design_combat<-design_combat[match(colnames(combat_input), design_combat$code),]
batch<-unlist(design_combat$Plex[match(colnames(combat_input), design_combat$code)])
design_combat$Drug <- factor(design_combat$Drug, levels = c("WT", "Vermurafenib_1uM", "Trametinib_10nM", "vemurafenib+trametinib"))
design_combat$Genetic <- factor(design_combat$Genetic, levels = c("WT", "ARID1A_KO"))

mod <- model.matrix(~ Drug + Genetic, data = design_combat)
data_combat <- ComBat(dat = as.matrix(combat_input), batch = batch, mod = mod, par.prior = TRUE)
data_combat <- as.data.frame(data_combat)

# ComBat introduces some negative corrected counts that we need to fix
data_combat <-  data_combat[apply(data_combat, 1, function(x) all(x >= 0)), ] 


####WRITE RESULTS####
#write.csv(data_combat, file = "./data/proteomic/processed/combat_peptide.csv")
