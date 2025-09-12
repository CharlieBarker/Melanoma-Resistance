library(purrr)
library(stringr)

format_peptide <- function(entry) {
  # Split entry into parts
  parts <- str_split(entry, "__")[[1]]
  seq <- parts[1]                        # peptide sequence
  mod <- parts[length(parts)]            # last part (e.g., "S11" or "T13" or "Y7")
  
  # Extract AA and position
  aa <- str_extract(mod, "^[A-Z]")                 # phosphorylated residue
  pos <- as.integer(str_extract(mod, "[0-9]+"))    # numeric position
  
  if (is.na(pos)) return(NA)   # safeguard
  
  # 10-mer window with site at position 6 (1-based)
  start <- pos - 5
  end   <- pos + 4
  
  # Build 10-mer, padding with "_" if out of bounds
  window <- vapply(start:end, function(i) {
    if (i < 1 || i > nchar(seq)) {
      "_"
    } else {
      substr(seq, i, i)
    }
  }, character(1))
  
  paste(window, collapse = "")
}


path_data<-"~/Desktop/Melanoma_Resistance/" #SWISS-plot
setwd(path_data)

phospho<-read.csv(file = "./data/proteomic/processed/combat_peptide.csv")

#Get experimental design 
design<-read.csv(file = "./data/proteomic/design.csv",header = T)
design$code<-paste0("X", design$TMT.labels,".pl", design$Plex)
design$Drug[design$Drug == "WT"] = "Untreated"
design$total<-paste0(design$Drug,"__", design$Genetic)

rownames(phospho)<-phospho$X
phosphorylated_peptides<-sapply(phospho$X, format_peptide)

phospho$X<- unname(phosphorylated_peptides)
sty_phospho<-phospho[!is.na(phospho$X),]              #removes NAs (peptides where we are unsure of the residue)
sty_phospho<-sty_phospho[!duplicated(sty_phospho$X),] #removes duplicates (condensing to 10 peptides makes some duplocated)
rownames(sty_phospho)<-sty_phospho$X
sty_phospho$X<-NULL
colnames(sty_phospho)<-gsub(colnames(sty_phospho), pattern = 'plex', replacement = 'pl')

#name columns after biological conditions, 
colnames(sty_phospho)<-design$total[match(colnames(sty_phospho), design$code)]
colnames(sty_phospho)<-str_remove(colnames(sty_phospho), pattern = "\\..+")


process_columns<-function(counts, replacement_Vec){
  # Extract drug names and genetic elements from column names
  col_names <- colnames(counts)
  drug_names <- sub("^(.*?)__.*", "\\1", col_names)
  gene_names <- sub(".*?__(.*)", "\\1", col_names)
  # Replace drug names with the names in replacement vector
  new_drug_names <- names(replacement_Vec)[match(drug_names, unname(replacement_Vec))]
  # Combine drug names and genetic elements to form new column names
  snames <- paste(new_drug_names, gene_names, sep = ".")
  # Assign new column names
  colnames(counts) <- snames
  out<-list()
  out$counts <- counts
  out$drug_names <- new_drug_names
  out$gene_names <- gene_names
  return(out)
}

####Produce volcano plots####


out<-process_columns(sty_phospho, replacement_Vec)
counts<-out$counts
counts <- log2(counts + 1)


drug_names <- out$drug_names
gene_names <- out$gene_names

drug <- factor(drug_names)       # "Untreated", "Vemurafenib", "Trametinib", "Combinations"
genotype <- factor(gene_names)   # "WT", "ARID1A_KO"

# Make sure Untreated and WT are the reference levels
genotype <- relevel(genotype, ref = "WT")
drug <- relevel(drug, ref = "Untreated")

design <- model.matrix(~ genotype * drug)
colnames(design) <- make.names(colnames(design))
# [1] "Intercept" "genotypeARID1A_KO" "drugCombinations" "drugTrametinib" "drugVemurafenib" "genotypeARID1A_KO.drugCombinations" ...

# Fit model with eBayes
fit <- lmFit(counts, design)
fit <- eBayes(fit)

cont.matrix <- makeContrasts(
  ARID1A_KO_comb_vs_untreated = drugCombinations + genotypeARID1A_KO.drugCombinations,
  WT_comb_vs_untreated = drugCombinations,
  levels = design
)

fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2)

# Extract results
ARID1A_KO_comb_vs_untreated<-topTable(fit2, coef="ARID1A_KO_comb_vs_untreated", sort.by = "P", n = Inf)
WT_comb_vs_untreated<-topTable(fit2, coef="WT_comb_vs_untreated", sort.by = "P", n = Inf)

# Make sure your object is a data.frame
df <- ARID1A_KO_comb_vs_untreated
df$Gene <- rownames(df)
# Extract peptide sequences + logFC
seqrnk <- df[, c("Gene", "logFC")]

# Sort by logFC, decreasing
seqrnk <- seqrnk[order(-seqrnk$logFC), ]
seqrnk$logFC <- round(seqrnk$logFC, digits = 5)
# Write to file (tab-separated, no headers/quotes, rownames suppressed)
write.table(
  seqrnk,
  file = "results/lfc_phosX/ARID1A_KO_comb_vs_untreated.seqrnk",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  col.names = FALSE
)