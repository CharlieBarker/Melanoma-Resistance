library(reshape2)
library(plyr)

RTKs<-c("EGFR", "FGFR2", "FGFR1", "IGF1R", "INSR", "FYN", "NTRK2", "MET",  #receptor tyrosine.
        "ABL2", "TXK", "FLT1") #non receptor tyrosine.
STKs<-c("MAPK3", "MAPK1", "PRKD1", "JUN", "DAPK1", "MAPK8", "MAPK9", "MAPK10") #serine/threonine kinases.

kinase_list<-rbind(data.frame(gene_name=RTKs, type="tyrosine kinase"),
                   data.frame(gene_name=STKs, type="serine/threonine kinase"))
conv_nodes<-data.frame(uniprt=V(g)$name,
                       gene_name=V(g)$Gene_name)
kinase_list$uniprt <- conv_nodes$uniprt[match(kinase_list$gene_name, conv_nodes$gene_name)]

#get lfc omics
mRNA_files<-list.files(path = "./results/lfc/mRNA", full.names = T, recursive = T)
protein_files<-list.files(path = "./results/lfc/protein", full.names = T, recursive = T)

mRNA_list<-list()
for (UKA_file in mRNA_files) {
  mRNA_results<-read.csv(file =UKA_file)
  mRNA_list[[UKA_file]]<-mRNA_results
}
complete_mRNA <- bind_rows(mRNA_list, .id = "file")
subset_mRNA<-complete_mRNA[complete_mRNA$X %in% kinase_list$gene_name,]
# Create matrices for logFC and P.Value
mRNA_logFC_matrix <- dcast(subset_mRNA, file ~ Gene, value.var = "logFC")
mRNA_pvalue_matrix <- dcast(subset_mRNA, file ~ Gene, value.var = "adj.P.Val")

protein_list<-list()
for (UKA_file in protein_files) {
  protein_results<-read.csv(file =UKA_file)
  protein_list[[UKA_file]]<-protein_results
}
complete_protein <- bind_rows(protein_list, .id = "file")
subset_protein<-complete_protein[complete_protein$X %in% kinase_list$gene_name,]
# Create matrices for logFC and P.Value
protein_logFC_matrix <- dcast(subset_protein, file ~ Gene, value.var = "logFC")
protein_pvalue_matrix <- dcast(subset_protein, file ~ Gene, value.var = "adj.P.Val")

library(ComplexHeatmap)
rownames(protein_logFC_matrix)<-protein_logFC_matrix$file
protein_logFC_matrix$file<-NULL

rownames(mRNA_logFC_matrix)<-mRNA_logFC_matrix$file
mRNA_logFC_matrix$file<-NULL

# Merge the data frames with full outer join
merged_protein_logFC <- rbind.fill(data.frame(protein_logFC_matrix), 
                                   data.frame(mRNA_logFC_matrix))
rownames(merged_protein_logFC)<-c(rownames(protein_logFC_matrix), rownames(mRNA_logFC_matrix))
# Split the merged data frame back into two matrices
protein_logFC_matrix <- merged_protein_logFC[rownames(protein_logFC_matrix),]
mRNA_logFC_matrix <- merged_protein_logFC[rownames(mRNA_logFC_matrix),]


#get kinomics
kinomics_files<-list.files(path = "./data/kinomics//", full.names = T, recursive = T)

results_list<-list()
for (UKA_file in kinomics_files[grep(kinomics_files, pattern = "vs")]) {
  kinomics_results<-read_xlsx(path =UKA_file)
  results_list[[UKA_file]]<-kinomics_results
}
# Apply the conversion to each data frame in the list
results_list <- lapply(results_list, function(x) {
  x$`SD Kinase Statitistic` <- as.numeric(x$`SD Kinase Statitistic`)
  return(x)
})

# Bind the rows
complete_results <- bind_rows(results_list, .id = "file")
complete_results$experiment<-unlist(map(str_split(complete_results$file, pattern = "/"), last))
complete_results$background <- unlist(map(str_split(complete_results$file, pattern = "/"), 6))

# Assuming `complete_results` is your data frame containing the required columns

# Reorder Kinase Name by Mean Kinase Statistic
complete_results$`Kinase Name` <- reorder(complete_results$`Kinase Name`, complete_results$`Mean Kinase Statistic`)
# Retrieve gene names
genename_df <- AnnotationDbi::select(EnsDb.Hsapiens.v86, 
                                     keys = as.character(complete_results$`Kinase Uniprot ID`), 
                                     keytype = "UNIPROTID", 
                                     columns = "GENENAME")

complete_results$Gene <- genename_df$GENENAME[match(complete_results$`Kinase Uniprot ID`, genename_df$UNIPROTID)]
kinomics_statistic_matrix <- dcast(complete_results, experiment ~ Gene, value.var = "Median Kinase Statistic")
kinomics_statistic_matrix<-kinomics_statistic_matrix[grepl(kinomics_statistic_matrix$experiment, pattern = " vs Untreated"),]
kinomics_statistic_matrix<-kinomics_statistic_matrix[!grepl(kinomics_statistic_matrix$experiment, pattern = "ARID1A"),]
rownames(kinomics_statistic_matrix)<-kinomics_statistic_matrix$experiment
kinomics_statistic_matrix$experiment<-NULL
kinomics_statistic_matrix<-kinomics_statistic_matrix[,colnames(kinomics_statistic_matrix) %in% kinase_list$gene_name]

# Merge the data frames with full outer join
double_merged_protein_logFC <- rbind.fill(data.frame(merged_protein_logFC), 
                                   data.frame(kinomics_statistic_matrix))
rownames(double_merged_protein_logFC)<-c(rownames(merged_protein_logFC), rownames(kinomics_statistic_matrix))
# Split the merged data frame back into two matrices
kinomics_statistic_matrix <- double_merged_protein_logFC[rownames(kinomics_statistic_matrix),]
groups<-c("Other kinases",  # ABL2
  "Other kinases",  # DAPK1
  "Receptor Tyrosine Kinases (RTKs)",  # EGFR
  "Receptor Tyrosine Kinases (RTKs)",  # FGFR1
  "Receptor Tyrosine Kinases (RTKs)",  # FGFR2
  "Receptor Tyrosine Kinases (RTKs)",  # FLT1
  "Other kinases",  # FYN
  "Receptor Tyrosine Kinases (RTKs)",  # IGF1R
  "Receptor Tyrosine Kinases (RTKs)",  # INSR
  "Mitogen-Activated Protein Kinases (MAPKs)",  # MAPK1
  "Mitogen-Activated Protein Kinases (MAPKs)",  # MAPK10
  "Mitogen-Activated Protein Kinases (MAPKs)",  # MAPK3
  "Mitogen-Activated Protein Kinases (MAPKs)",  # MAPK8
  "Mitogen-Activated Protein Kinases (MAPKs)",  # MAPK9
  "Receptor Tyrosine Kinases (RTKs)",  # MET
  "Other kinases",  # PRKD1
  "Other kinases",  # JUN
  "Receptor Tyrosine Kinases (RTKs)",  # NTRK2
  "Other kinases")  # TXK


# Assuming `protein_logFC_matrix` is already defined

# Remove rows related to ARID1A and rename remaining rows
protein_logFC_matrix <- protein_logFC_matrix[!grepl("arid1a_lfc.csv", rownames(protein_logFC_matrix)), ]
rownames(protein_logFC_matrix) <- c("Combined vs Untreated", "Trametinib vs Untreated", "Vemurafenib vs Untreated")

# Assuming `mRNA_logFC_matrix` is already defined

# Remove rows related to ARID1A and rename remaining rows
mRNA_logFC_matrix <- mRNA_logFC_matrix[!grepl("arid1a_lfc.csv", rownames(mRNA_logFC_matrix)), ]
rownames(mRNA_logFC_matrix) <- c("Combined vs Untreated", "Trametinib vs Untreated", "Vemurafenib vs Untreated")
# Assuming `kinomics_statistic_matrix` is already defined

# Remove rows related to ARID1A and rename remaining rows
kinomics_statistic_matrix <- kinomics_statistic_matrix[!grepl("arid1a_lfc.csv", rownames(kinomics_statistic_matrix)), ]
rownames(kinomics_statistic_matrix) <- c("Combined vs Untreated", "Trametinib vs Untreated", "Vemurafenib vs Untreated")



ht1 = Heatmap(data.matrix(protein_logFC_matrix), name = "Protein Abundance (LFC)", na_col = "darkgrey", 
              rect_gp = gpar(col = "white", lwd = 2), cluster_columns = F, cluster_rows = F,
              column_split = groups, 
              column_gap=unit(.05, "npc"),
              row_title = "Protein", row_title_gp = gpar(fontsize = 10),column_title_gp = gpar(fontsize = 10))
ht2 = Heatmap(data.matrix(mRNA_logFC_matrix), name = "mRNA abundance (LFC)", na_col = "darkgrey", 
              rect_gp = gpar(col = "white", lwd = 2), cluster_columns = F, cluster_rows = F,
              column_split = groups, 
              column_gap=unit(.05, "npc"),
              row_title = "mRNA", row_title_gp = gpar(fontsize = 10),column_title_gp = gpar(fontsize = 10))
ht3 = Heatmap(data.matrix(kinomics_statistic_matrix), name = "Kinase activity (LFC)", na_col = "darkgrey", 
              rect_gp = gpar(col = "white", lwd = 2), cluster_columns = F, cluster_rows = F,
              column_split = groups, 
              column_gap=unit(.05, "npc"),
              row_title = "Kinase activity", row_title_gp = gpar(fontsize = 10),column_title_gp = gpar(fontsize = 10))
pdf(# The directory you want to save the file in
  width = 12, # The width of the plot in inches
  height = 4,
  file = "~/Desktop/test.pdf")

ht_list = ht1 %v% ht2 %v% ht3
draw(ht_list)


dev.off()

