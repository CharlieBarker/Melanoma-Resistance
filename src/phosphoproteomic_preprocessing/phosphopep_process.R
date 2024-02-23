suppressPackageStartupMessages({
  library(PhosR)
  library(dplyr)
  library(ggplot2)
  library(GGally)
  library(ggpubr)
  library(calibrate)
  library(network)
  library(EnsDb.Hsapiens.v86)
  library(purrr)
  library(stringr)
  library(readr)
  library(gdata)
  library(data.table)
})
data("KinaseMotifs")
data("KinaseFamily")
data('PhosphoSitePlus')
data('PhosphoELM')

rename_rows<-function(to_plot, edb){
  uniprot_accession<-unlist(map(str_split(sub(" .*", "", rownames(to_plot)), "__"), 2))
  frag_sequence<-unlist(map(str_split(sub(" .*", "", rownames(to_plot)), "__"), 1))
  residue_in_frag<-parse_number(unlist(map(str_split(rownames(to_plot), "__"), last)))
  frag_sequence_first_residue<-parse_number(unlist(map(str_split(str_extract_all(rownames(to_plot), "\\[[^()]+\\]"), "-"), first)))  
  residue<-gsub('[[:digit:]]+', '', unlist(map(str_split(rownames(to_plot), "__"), last)))
  site<-frag_sequence_first_residue+residue_in_frag-1
  full_genename_df<-AnnotationDbi::select(edb, keys = uniprot_accession, keytype = "UNIPROTID", columns = "GENENAME")
  to_plot$future_rownames  <- paste0(full_genename_df$GENENAME[match(uniprot_accession , full_genename_df$UNIPROTID)], ";",
                                     residue,
                                     site, ";")
  return(to_plot)
}


path_data<-"~/Desktop/Melanoma_Resistance/" #SWISS-plot
setwd(path_data)


#Get experimental design 
design<-read.csv(file = "./data/proteomic/design.csv",header = T)
design$code<-paste0("X", design$TMT.labels,".pl", design$Plex)
design$Drug[design$Drug == "WT"] = "Untreated"
design$total<-paste0(design$Drug,"__", design$Genetic)

regressed_phospho<-read.csv(file = "./data/proteomic/processed/phosphopeptide__protein_adjusted.csv")
rownames(regressed_phospho)<-regressed_phospho$X

regressed_phospho$X<-NULL
#mark the phosphorylated residue as an X on the sequence
#this is necessary so we can centre the sequence on the phosphorylated protein 
fragment_rownames<-list()
for (fragment in rownames(regressed_phospho)) {
  index<-readr::parse_number(unlist(map(str_split(fragment, pattern = "__"), last)))
  fragment_ind<-str_remove(fragment, pattern = "__.+")
  substr(fragment_ind, index, index) <- "X"
  fragment_info<-unlist(str_split(fragment, pattern = "__"))[-1]
  new_fragment<-paste0(fragment_ind, "__",fragment_info[1], "__", fragment_info[2])
  fragment_rownames[[fragment]]<-new_fragment
}
all(names(fragment_rownames) == rownames(regressed_phospho) )
rownames(regressed_phospho) <- unname(fragment_rownames)
to_study<-rename_rows(regressed_phospho, EnsDb.Hsapiens.v86)

sequence_mapping<-data.frame(sequence = unlist(map(str_split(rownames(to_study), "__"), 1)),
                             psite = to_study$future_rownames)
uniq_psites<-unique(sequence_mapping$psite)
seq_psite_list<-lapply(uniq_psites, function(x){sequence_mapping$sequence[sequence_mapping$psite==x]})
names(seq_psite_list)<-uniq_psites
#get the overlap between the psites that have more than one sequence associated with them. 
sequence_to_psite<-lapply(seq_psite_list, function(x){paste(Reduce(union, sapply(x, strsplit, "")), collapse = "")})
#very slight differences in the fragments - we should be doing this earlier -- really. 

# Extract site sequences from the list, replace 'X' with the correct residue, and create a new dataframe
SiteSequence = sapply(names(sequence_to_psite), function(x) {
  parts <- strsplit(x, ";")[[1]]
  protein <- parts[1]
  residue <- parts[2]
  fragment <- gsub("X", paste0(tolower(substr(residue, 1, 1)), "*"), sequence_to_psite[[x]])
  fragment
})


prot_mapper<-data.frame(ID=unlist(map(str_split(unlist(map(str_split(string = rownames(to_study),
                                                                     pattern = "__"), 2)), pattern = "\\s"), 1)),
                        type="uniprot",
                        AA=gsub('[[:digit:]]+', '', unlist(map(str_split(rownames(to_study), pattern = "__"), 3))),
                        residue=as.integer(parse_number(unlist(map(str_split(to_study$future_rownames, pattern = ";"), 2)))))
# write_delim(x = prot_mapper,
#             file = "~/MelanomaProject/data/prot_mapper_input.csv",
#             col_names = F, 
#             delim = ",")

to_study<-aggregate(. ~ future_rownames, data=to_study, FUN=sum)
rownames(to_study)<-to_study$future_rownames
to_study$future_rownames<-NULL

#name columns after biological conditions, 
colnames(to_study)<-design$total[match(colnames(to_study), design$code)]
to_study<-to_study#[!grepl(colnames(to_study), pattern = "MED12") & !grepl(colnames(to_study), pattern = "ARID1A")]
colnames(to_study)<-str_remove(colnames(to_study), pattern = "\\..+")

#subset for relevant genetic background. 
to_study<-to_study#[grep(colnames(to_study), pattern = "__WT")]
colnames(to_study)<-str_remove(colnames(to_study), pattern = "\\..+")

# filter for up-regulated phosphosites
phosphoL6.mean <- meanAbundance(to_study, grps = colnames(to_study))
aov <- matANOVA(mat=to_study, grps=colnames(to_study))
idx <- (aov < 0.1) & (rowSums(abs(phosphoL6.mean) > .5) > 0) #idx is the filter 0.5 and aov of 0.05
without_dup.reg <- to_study[idx, ,drop = FALSE]

dim(without_dup.reg)

without_dup.reg.std <- standardise(without_dup.reg)
kinase_substrate_input<-without_dup.reg.std
write.csv(kinase_substrate_input, file = "./data/input_data/phosphosites.csv")
