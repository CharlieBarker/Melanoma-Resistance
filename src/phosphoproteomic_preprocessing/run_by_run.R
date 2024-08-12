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

plex_to_study<- 3

combined_peptides<-read.csv(file = "./data/proteomic/raw/peptide_total.csv")
rownames(combined_peptides)<-combined_peptides$X
combined_peptides$X<-NULL
combined_peptides[,grep("131C", colnames(combined_peptides))]<-NULL 

# read the design matrix
design<-read.csv(file = "./data/proteomic/design.csv",header = T)
grps <- design$total[match(colnames(combined_peptides),design$code)]

#remove MED12
combined_peptides<-combined_peptides[,!grepl(grps, pattern = "MED12")]
design<-design[!grepl(design$total, pattern = "MED12"),]
grps<-grps[!grepl(grps, pattern = "MED12")]

#seperate plexes 
combined_peptides<-combined_peptides[,grepl(colnames(combined_peptides), pattern = paste0("plex", plex_to_study))]
design<-design[design$Plex == plex_to_study,]


####Sample loading normalisation####

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
phos_raw<-na.omit(SummarizedExperiment::assay(ppe,"Quantification"))


#protein time 

# read the Supplemental 01 file (saved as a CSV export from XLSX file)
data_start <- read_csv("./data/proteomic/raw/Proteins_grouped_abundances.csv")
design_protein<-read.csv(file = "./data/proteomic/design.csv",header = T)
design_protein$code<-paste0(design_protein$TMT.labels,".pl", design_protein$Plex)
design_protein$total<-paste0(design_protein$Drug,"__", design_protein$Genetic)
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

grps <- design_protein$total[match(colnames(combided_proteins),design_protein$code)]

#remove MED12
combided_proteins<-combided_proteins[,!grepl(grps, pattern = "MED12")]
design_protein<-design_protein[!grepl(design_protein$total, pattern = "MED12"),]
grps<-grps[!grepl(grps, pattern = "MED12")]
#seperate plexes 
combided_proteins<-combided_proteins[,grepl(colnames(combided_proteins), pattern = paste0("pl", plex_to_study))]
design_protein<-design_protein[design_protein$Plex == plex_to_study,]


####Sample loading normalisation####

combided_proteins<-normSLN(combided_proteins)
prot_data<-combided_proteins
#extract phospho-site information
site = readr::parse_number(unlist(map(str_split(rownames(prot_data), pattern = "__"), last)))
uniprot_id = rownames(prot_data)
genename_df<-AnnotationDbi::select(EnsDb.Hsapiens.v86, keys = uniprot_id, keytype = "UNIPROTID", columns = "GENENAME")
geneSymbol <- genename_df$GENENAME

counts <- table(design$total[match(colnames(prot_data), design$code)])
ppe_protein <- PhosphoExperiment(assays = list(Quantification = as.matrix(prot_data)), 
                         GeneSymbol = geneSymbol)
prot_raw<-na.omit(SummarizedExperiment::assay(ppe_protein,"Quantification"))


####Peptide/Protein adjustment####


#sum the peptides with the slightly different things that we are not interested in.
phos_raw<-data.frame(phos_raw)
prot_raw<-data.frame(prot_raw)

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

residual_list<-apply(combine,1,lm_fun)
regressed_phospho<-data.frame(t(residual_list))

#the apply hasnt affected the order of our numbers so long as our sanity checks have passed 
colnames(regressed_phospho)<-df_info$code
colnames(prot_sl_tmm)<-df_info$code
colnames(phos_sl_tmm)<-df_info$code

colnames(regressed_phospho) <- design$Genetic[match(colnames(regressed_phospho), design$code)]
colnames(prot_sl_tmm) <- design$Genetic[match(colnames(prot_sl_tmm), design$code)]
colnames(phos_sl_tmm) <- design$Genetic[match(colnames(phos_sl_tmm), design$code)]

#Differential abundance calculation 



counts<-phos_sl_tmm
fit_levels<-c("groupARID1A_KO", "groupWT") 
contr_arid1a <- makeContrasts(groupARID1A_KO - groupWT, levels = fit_levels)

#for ARID1A
contrast_matrix = contr_arid1a
group <- colnames(counts)

counts <- counts[!rowSums(is.infinite(as.matrix(counts))), ]

mm <- model.matrix(~0 + group)
fit <- lmFit(counts, mm)

contr <- contrast_matrix
tmp <- contrasts.fit(fit, contr)
tmp <- eBayes(tmp)
top.table <- topTable(tmp, sort.by = "P", n = Inf)

nsig<-length(which(top.table$adj.P.Val < 0.05))
nsig
top.table$Gene <- rownames(top.table)
top.table <- top.table[,c("Gene", names(top.table)[1:6])]



