library(limma) 
library(ggrepel)
library(EnsDb.Hsapiens.v86)
library(ggpubr)
library(enrichR)
library(readr)
library(stringr)
library(purrr)
library(data.table)
library(dplyr)
rename_rows<-function(to_plot, edb){
  uniprot_accession<-unlist(map(str_split(sub(" .*", "", rownames(to_plot)), "__"), 2))
  frag_sequence<-unlist(map(str_split(sub(" .*", "", rownames(to_plot)), "__"), 1))
  residue_in_frag<-parse_number(unlist(map(str_split(rownames(to_plot), "__"), last)))
  frag_sequence_first_residue<-parse_number(unlist(map(str_split(str_extract_all(rownames(to_plot), "\\[[^()]+\\]"), "-"), 1)))  
  residue<-gsub('[[:digit:]]+', '', unlist(map(str_split(rownames(to_plot), "__"), last)))
  site<-frag_sequence_first_residue+residue_in_frag-1
  to_plot$future_rownames  <- paste0(uniprot_accession, ";",
                                     residue,
                                     site, ";")
  return(to_plot)
}

#path_data<-"~/MelanomaProject/data/phosphoproteomics/main/" #original
path_data<-"/home/charlie/MelanomaProject/data/phosphoproteomics/main/swiss-prot/total/" #SWISS-plot
setwd(path_data)

#Get experimental design 
design<-read.csv(file = "~/MelanomaProject/data/phosphoproteomics/main/design.csv",header = T)
design$code<-paste0("X", design$TMT.labels,".pl", design$Plex)
design$total<-paste0(design$Drug,"__", design$Genetic)

adjusted_phospho<-read.csv(file = "./phosphopeptide__protein_adjusted.csv")
unadjusted_phospho<-read.csv(file = "./phosphopeptide__norm_combat.csv")
colnames(unadjusted_phospho)<-str_replace(string = colnames(unadjusted_phospho), pattern = "plex", replacement = "pl")

proteins<-read.csv(file = "./protein__norm_combat_imp.csv")

prepare_mofa_input<-function(input, is.phos=T){
  rownames(input)<-input$X
  input$X<-NULL
  to_Study<-na.omit(input)
  
  if (is.phos) {
    #aggregate fragments covering the same psite 
    to_Study<-rename_rows(to_Study)
    to_Study<-aggregate(. ~ future_rownames, data=to_Study, FUN=sum)
    rownames(to_Study)<-to_Study$future_rownames
    to_Study$future_rownames<-NULL
  }

  colnames(to_Study)<-design$total[match(colnames(to_Study), design$code)]
  #make names for mofa
  mofa_in<-to_Study
  colnames(mofa_in)<-str_replace(string = colnames(mofa_in), pattern = "WT__", replacement = "Untreated__")
  colnames(mofa_in)<-str_replace(string = colnames(mofa_in), pattern = "vermurafenib.trametinib", replacement = "vermurafenib_and_trametinib")
  colnames(mofa_in) <-str_replace(colnames(mofa_in), pattern = "\\+", replacement = "_and_")
  colnames(mofa_in)<-make.unique(colnames(mofa_in))
  return(mofa_in)
}

mofa_adjusted_phos<-prepare_mofa_input(adjusted_phospho)
mofa_unadjusted_phos<-prepare_mofa_input(unadjusted_phospho)
mofa_protein<-prepare_mofa_input(proteins, is.phos = F)

write.csv(mofa_adjusted_phos, file = "~/MelanomaProject/mofa/input_data/phosphosites_adjusted.csv")
write.csv(log2(mofa_unadjusted_phos), file = "~/MelanomaProject/mofa/input_data/phosphosites_unadjusted.csv")
write.csv(log2(mofa_protein), file = "~/MelanomaProject/mofa/input_data/proteins_unadjusted.csv")

