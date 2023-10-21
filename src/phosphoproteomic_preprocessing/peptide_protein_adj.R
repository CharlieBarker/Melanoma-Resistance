library(tidyverse)
library(limma) 
library(edgeR) 
library(sva)
library(psych)


path_data<-"~/Desktop/Melanoma_Resistance/" #SWISS-plot
setwd(path_data)

####LOAD FUNCTIONS####

source("./src/functions/SLN_normalisation.R")

# regress phospho to total proteome to estimate the "NET" phosphorylation. Written by Lourdes. 
lm_fun<-function(x){
  y<-as.numeric(unlist(x[grep(pattern = "__PHOS", names(x))]))
  y1<-as.numeric(unlist(x[grep(pattern = "__PROT", names(x))]))
  y<-log2(y)
  y1<-log2(y1)
  if(all(is.na(y)) | all(is.na(y1)))
    fill<-c(rep(NA,length(y)))
  else(fill<-residuals(lm(y~y1)))
  return(fill)
}


proteome_start<-read.csv(file = "./data/proteomic/processed/protein__norm_combat_imp.csv")
phosphoproteome_start<-read.csv(file = "./data/proteomic/processed/combat_peptide.csv")
#Phosphopeptide pre processing 
#get all modifications 

#sum the peptides with the slightly different things that we are not interested in.

prot_raw <-proteome_start
rownames(prot_raw)<-prot_raw$X
prot_raw$X<-NULL
phos_raw <- data.frame(phosphoproteome_start)
rownames(phos_raw) <- phos_raw$X
phos_raw$X<-NULL
colnames(phos_raw)<-str_replace(colnames(phos_raw), "plex", "pl")

#Get design matrix
design_file = "./data/proteomic/design.csv"
design<-read.csv(file = design_file,header = T)
design$code<-paste0("X", design$TMT.labels, ".pl", design$Plex)
design$total<-paste0(design$Genetic, "___" ,design$Drug)

df_info = S4Vectors::DataFrame(design[match(colnames(phos_raw), design$code), ])
all(rownames(df_info) == df_info$code) #write exception to bring up and error if this is false
rownames(df_info) = colnames(phos_raw)
prot_raw<-prot_raw[,colnames(prot_raw) %in% colnames(phos_raw)] #only for imputed
all(colnames(phos_raw) == colnames(prot_raw)) #write exception to bring up and error if this is false
#add labels so we can distinguish the colomns 
colnames(prot_raw)<- paste0(colnames(prot_raw), "__PROT")
colnames(phos_raw) <- paste0(colnames(phos_raw), "__PHOS")

#next we normalise the proteome and the peptide fragments using identical normalisation procedures. 
cut_off<-50
#proteome
prot_naless<-na.omit(prot_raw)
prot_naless[prot_naless<cut_off]<-NA
prot_sl_tmm<-na.omit(prot_naless)
#prot_sl_tmm<-normSLN(prot_sl_tmm)



#peptide fragments
phos_naless<-na.omit(phos_raw)
phos_naless[phos_naless<cut_off]<-NA
phos_sl_tmm<-na.omit(phos_naless)
#phos_sl_tmm<-normSLN(phos_sl_tmm)

#check that the column names are in the correct order after all of this 
all(unlist(map(str_split(colnames(phos_sl_tmm), pattern = "__"), 1)) == unlist(map(str_split(colnames(prot_sl_tmm), pattern = "__"), 1)))

#combine the phosphoproteome and proteome 
fragment_in_master<-unlist(map(str_split(rownames(phos_sl_tmm), pattern = "__"), 2))
phospho_accession<-unlist(map(str_split(fragment_in_master, pattern = "\\s"), 1))
match_indices<-match(phospho_accession,
                     rownames(prot_sl_tmm))
combine<-cbind(phos_sl_tmm, prot_sl_tmm[match_indices,])
combine<-as.data.frame(combine)

#test lm works 
test<-lm(as.numeric(unlist(combine[200,grep(pattern = "__PHOS", colnames(combine))])) ~ as.numeric(unlist(combine[200,grep(pattern = "__PROT", colnames(combine))])))
#combine the proteome and the phosphoproteome, by normalising the two (logfold change)  
#we can use simple log change or we can use the function we tested above 
residual_list<-apply(combine,1,lm_fun)
regressed_phospho<-data.frame(t(residual_list))

#the apply hasnt affected the order of our numbers so long as our sanity checks have passed 
colnames(regressed_phospho)<-df_info$code
colnames(prot_sl_tmm)<-df_info$code
colnames(phos_sl_tmm)<-df_info$code
#you should also check the correlation of these seperately - i used sln_phos.R script for this. 
#also the same script checks for the abundance of MED12 and ARID1A... ARID1A is ko'd but MED12 isnt really.

write.csv(regressed_phospho, file = "./data/proteomic/processed/phosphopeptide__protein_adjusted.csv")