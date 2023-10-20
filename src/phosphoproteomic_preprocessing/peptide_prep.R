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
ppe_imputed_tmp <- scImpute(ppe_filtered, 0.5, grps)[,colnames(ppe_filtered)]

# Paired tail-based imputation
ppe_imputed <- ppe_imputed_tmp
#Impute between Trametinib treatments 
tram_ARID1A_log<-design$total[match(colnames(ppe_imputed), design$code)] == "Trametinib_10nM__ARID1A_KO"
tram_MED12_log<-design$total[match(colnames(ppe_imputed), design$code)] == "Trametinib_10nM__MED12_KO"
tram_WT_log<-design$total[match(colnames(ppe_imputed), design$code)] == "Trametinib_10nM__WT"

ver_ARID1A_log<-design$total[match(colnames(ppe_imputed), design$code)] == "Vermurafenib_1uM__ARID1A_KO"
ver_MED12_log<-design$total[match(colnames(ppe_imputed), design$code)] == "Vermurafenib_1uM__MED12_KO"
ver_WT_log<-design$total[match(colnames(ppe_imputed), design$code)] == "Vermurafenib_1uM__WT"

combo_ARID1A_log<-design$total[match(colnames(ppe_imputed), design$code)] == "vemurafenib+trametinib__ARID1A_KO"
combo_MED12_log<-design$total[match(colnames(ppe_imputed), design$code)] == "vemurafenib+trametinib__MED12_KO"
combo_WT_log<-design$total[match(colnames(ppe_imputed), design$code)] == "vemurafenib+trametinib__WT"

WT_ARID1A_log<-design$total[match(colnames(ppe_imputed), design$code)] == "WT__ARID1A_KO"
WT_MED12_log<-design$total[match(colnames(ppe_imputed), design$code)] == "WT__MED12_KO"
WT_WT_log<-design$total[match(colnames(ppe_imputed), design$code)] == "WT__WT"

#impute is done within the drug treatments, as these are where the most variation occurs. 
pc1=.6
pc2=.2
paired=F
# #Impute between Trametinib treatments
ppe_imputed[,tram_ARID1A_log] <- ptImpute(ppe_imputed[,tram_WT_log], #control
                                          ppe_imputed[,tram_ARID1A_log], #ARID1A_Trametinib
                                          percent1 = pc1, percent2 = pc2, paired = paired)
ppe_imputed[,tram_MED12_log] <- ptImpute(ppe_imputed[,tram_WT_log], #control
                                         ppe_imputed[,tram_MED12_log], #MED12_Trametinib
                                         percent1 = pc1, percent2 = pc2, paired = paired)
#Impute between Vermurafenib treatments
ppe_imputed[,ver_ARID1A_log] <- ptImpute(ppe_imputed[,ver_WT_log], #control
                                         ppe_imputed[,ver_ARID1A_log], #ARID1A_Vermurafenib
                                         percent1 = pc1, percent2 = pc2, paired = paired)
ppe_imputed[,ver_MED12_log] <- ptImpute(ppe_imputed[,ver_WT_log], #control
                                        ppe_imputed[,ver_MED12_log], #MED12_Vermurafenib
                                        percent1 = pc1, percent2 = pc2, paired = paired)
# # #Impute between Combination treatments
ppe_imputed[,combo_ARID1A_log] <- ptImpute(ppe_imputed[,combo_WT_log], #control
                                         ppe_imputed[,combo_ARID1A_log], #ARID1A_Vermurafenib
                                         percent1 = pc1, percent2 = pc2, paired = paired)
ppe_imputed[,combo_MED12_log] <- ptImpute(ppe_imputed[,combo_WT_log], #control
                                        ppe_imputed[,combo_MED12_log], #MED12_Vermurafenib
                                        percent1 = pc1, percent2 = pc2, paired = paired)

# #Impute between WT - care! we only have 2 replicates here.

ppe_imputed[,WT_ARID1A_log] <- ptImpute(ppe_imputed[,WT_WT_log], #control
                                        ppe_imputed[,WT_ARID1A_log], #ARID1A_Vermurafenib
                                        percent1 = pc1, percent2 = pc2, paired = paired)
ppe_imputed[,WT_MED12_log] <- ptImpute(ppe_imputed[,WT_WT_log], #control
                                       ppe_imputed[,WT_MED12_log], #MED12_Vermurafenib
                                       percent1 = pc1, percent2 = pc2, paired = paired)



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
design_combat$Genetic <- factor(design_combat$Genetic, levels = c("WT", "ARID1A_KO", "MED12_KO"))

mod <- model.matrix(~ Drug + Genetic, data = design_combat)
data_combat <- ComBat(dat = as.matrix(combat_input), batch = batch, mod = mod, par.prior = TRUE)
data_combat <- as.data.frame(data_combat)

# ComBat introduces some negative corrected counts that we need to fix
data_combat <-  data_combat[apply(data_combat, 1, function(x) all(x >= 0)), ] 


####WRITE RESULTS####
write.csv(data_combat, file = "./data/proteomic/processed/combat_peptide.csv")
