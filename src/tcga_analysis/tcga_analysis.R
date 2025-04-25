# Load necessary libraries
library(dplyr)
library(ggplot2)
library(cowplot)
library(tidyr)
library(decoupleR)
library(TCGAretriever)
library(edgeR)
library(limma)
library(reshape2)
library(ggrepel)
library(wesanderson)
library(tibble)

# Set working directory
home_dir <- "~/Desktop/Melanoma_Resistance/"
setwd(home_dir)

# Fetch and prepare data
fetch_and_prepare_data <- function(extract_new=F) {
  if (extract_new==T) {
    all_skcm_RNA <- fetch_all_tcgadata(case_list_id = 'skcm_tcga_all',
                                       gprofile_id = 'skcm_tcga_gdc_mrna_seq_fpkm',
                                       mutations = FALSE)
    counts <- all_skcm_RNA %>%
      # Remove unnecessary columns
      select(-entrezGeneId, -type) %>%
      # Set 'hugoGeneSymbol' as row names
      column_to_rownames(var = "hugoGeneSymbol")
    
    # Define the file path
    file_path <- "./data/tcga_analysis/skcm_tcga_all_mrna_seq_fpkm.csv"
    
    # Create directories if they don't exist
    dir.create(dirname(file_path), recursive = TRUE, showWarnings = FALSE)
    
    # Write the dataframe `df` to the CSV file
    write.csv(counts, file = file_path, row.names = FALSE)
    
    return(counts)
  }else{
    skcm_tcga <- get_case_lists("skcm_tcga")
    counts <- read.csv(file = "./data/tcga_analysis/skcm_tcga_all_mrna_seq_fpkm.csv", row.names = 1)
    return(counts)
  }
}

# Function to analyze ARID1A mutation status
get_arid1a_mutation_status <- function(gene_list) {
  arid1a_mut_data <- get_mutation_data(case_list_id = 'skcm_tcga_all',
                                       gprofile_id = 'skcm_tcga_mutations',
                                       glist = gene_list)
  arid1a_cnv_data <- get_molecular_data(case_list_id = 'skcm_tcga_all',
                                        gprofile_id = 'skcm_tcga_gistic',
                                        glist = gene_list)
  arid1a_cnv_data <- reshape2::melt(arid1a_cnv_data)
  arid1a_del<-as.character(arid1a_cnv_data[arid1a_cnv_data$value < 0, ]$variable)
  samples_arid1a_affected <- unique(c(arid1a_del,
                                      as.character(arid1a_mut_data$sampleId)))
  patients_out <- data.frame(samples = unique(gsub(samples_arid1a_affected, pattern = "-", replacement = ".")),
                             cn_deletion = unique(samples_arid1a_affected) %in% arid1a_del,
                             mutation = unique(samples_arid1a_affected) %in% as.character(arid1a_mut_data$sampleId))
  patients_out
}

# Perform TF Activity Analysis
perform_tf_analysis <- function(counts, net, samples_arid1a_affected) {
  sample_acts <- decoupleR::run_ulm(mat = counts, net = net, .source = 'source', .target = 'target',
                                    .mor = 'mor', minsize = 5)
  sample_acts$is_arid1a_mutant <- sample_acts$condition %in% samples_arid1a_affected
  sample_acts
}

# Main script execution
counts <- fetch_and_prepare_data()
samples_arid1a_affected <- get_arid1a_mutation_status(c("ARID1A"))
net <- get_collectri(organism = 'human', split_complexes = FALSE)
sample_acts <- perform_tf_analysis(counts, net, samples_arid1a_affected$samples)

#first venn diagram of patients
library(ggvenn)

x <- list(
  `ARID1A deletion` = samples_arid1a_affected$samples[samples_arid1a_affected$cn_deletion],
  `ARID1A mutation` = samples_arid1a_affected$samples[samples_arid1a_affected$mutation],
  `ARID1A WT` = colnames(counts)[!colnames(counts) %in% samples_arid1a_affected$samples]
)



pdf(# The directory you want to save the file in
  width = 6, # The width of the plot in inches
  height = 6,
  file = "./results/tcga_analysis/patient_venn.pdf")
ggvenn(
  x,
  fill_color = c("#0073C2FF", "#EFC000FF", "#868686FF"),
  stroke_size = 0.5, set_name_size = 4
)
dev.off()

#arid1a violin plot
data_subset <- sample_acts[sample_acts$source %in% c("JUN", "ARID1A", "MYC", "RFX5", "RFXANK", "TWIST1", "SMAD4"), ]
data_subset$is_arid1a_mutant <- paste0(data_subset$condition %in% samples_arid1a_affected$samples[samples_arid1a_affected$mutation],
                                       "__",
                                       data_subset$condition %in% samples_arid1a_affected$samples[samples_arid1a_affected$cn_deletion])
# Modify `is_arid1a_mutant` to have informative labels
data_subset$is_arid1a_mutant <- factor(
  paste0(data_subset$condition %in% samples_arid1a_affected$samples[samples_arid1a_affected$mutation],
         "__",
         data_subset$condition %in% samples_arid1a_affected$samples[samples_arid1a_affected$cn_deletion]),
  levels = c("FALSE__FALSE", "TRUE__FALSE", "FALSE__TRUE", "TRUE__TRUE"),
  labels = c("WT ARID1A", "ARID1A mutation", "ARID1A mutation or CN loss", "ARID1A mutation or CN loss")
)

# Create the violin plot
violin_plot<-ggplot(data_subset[data_subset$source %in% c("ARID1A", "RFX5", "RFXANK", "TWIST1"),],
                    aes(x = is_arid1a_mutant, y = score, color = is_arid1a_mutant)) +
  geom_violin(alpha = 1, width = 0.8, colour = "black") +                         # Violin plot
  geom_jitter(width = 0.2, alpha = 0.4, size = 1) +               # Add jitter points for individual values
  facet_wrap(~source, scales = "free", ncol = 2) +                # Facet by source, free scales, 4 columns
  cowplot::theme_cowplot() +
  theme(
    plot.title = element_text(size = 15, face = "bold"),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 1),
    axis.text.x = element_text(angle = 45, hjust = 1)             # Rotate x-axis labels to 45 degrees
  ) +
  grids(linetype = "dashed") +
  scale_color_manual(values=wes_palette("Darjeeling1")) +
  xlab("ARID1A Mutation Status") +                                # X-axis label
  ylab("Transcription Factor Activity Score") +                                                 # Y-axis label
  ggtitle("TF activity by melanoma patient") +
  geom_signif(comparisons = list(c("WT ARID1A", "ARID1A mutation")),
              map_signif_level=c("***"=0.01, "**"=0.1, "*"=0.25),
              test = "t.test", colour = "black") +
  geom_signif(comparisons = list(c("WT ARID1A", "ARID1A mutation or CN loss")),
            map_signif_level=c("***"=0.01, "**"=0.1, "*"=0.25),
            test = "t.test", colour = "black")


pdf(# The directory you want to save the file in
  width = 10, # The width of the plot in inches
  height = 12,
  file = "./results/tcga_analysis/tf_activity_by_melanoma_patient.pdf")
# Showing plot

violin_plot
dev.off()



#see the levels of these in ARID1A mutant patients.
library(edgeR)
# Reshape the data to have genes as columns and conditions as rows
wide_data <- data_subset %>%
  dplyr::select(source, condition, score) %>%
  pivot_wider(names_from = source, values_from = score)
wide_data$is_ARID1a_mutant <- wide_data$condition %in% samples_arid1a_affected$samples

meta_data<-data.frame(patient = wide_data$condition, is_ARID1A_affected = wide_data$is_ARID1a_mutant)
mm <- model.matrix(~0 + is_ARID1A_affected, data = meta_data)
dge <- DGEList(counts)
dge <- calcNormFactors(dge)
keep <- filterByExpr(counts, mm)
counts_to_keep <- dge[keep,]
y <- voom(counts_to_keep, mm, plot = T)
fit <- lmFit(y, mm)

contr <- makeContrasts(is_ARID1A_affectedTRUE - is_ARID1A_affectedFALSE, levels = colnames(coef(fit)))
tmp <- contrasts.fit(fit, contr)
tmp <- eBayes(tmp)
top.table <- topTable(tmp, adjust.method = "BH", sort.by = "P", n = Inf)


# Extract t-values per gene
deg <- top.table %>%
  dplyr::select(logFC, t, adj.P.Val) %>%
  dplyr::filter(!is.na(t)) %>%
  as.matrix()

write.csv(x = deg, file = "./results/tcga_lfc/arid1a_mutant_vs_wt.csv")
# Run ulm
contrast_acts <- decoupleR::run_ulm(mat = deg[, 't', drop = FALSE],
                                    net = net,
                                    .source = 'source',
                                    .target = 'target',
                                    .mor='mor',
                                    minsize = 10)

# Filter top TFs in both signs
f_contrast_acts <- contrast_acts %>%
  dplyr::mutate(rnk = NA)
msk <- f_contrast_acts$score > 0
f_contrast_acts[msk, 'rnk'] <- rank(-f_contrast_acts[msk, 'score'])
f_contrast_acts[!msk, 'rnk'] <- rank(-abs(f_contrast_acts[!msk, 'score']))
tfs <- f_contrast_acts %>%
  dplyr::arrange(rnk) %>%
  dplyr::pull(source)
f_contrast_acts <- f_contrast_acts %>%
  dplyr::filter(source %in% tfs)
colors <- rev(RColorBrewer::brewer.pal(n = 11, name = "RdBu")[c(2, 10)])


p <- ggplot2::ggplot(data = f_contrast_acts[f_contrast_acts$p_value < 0.05,],
                     mapping = ggplot2::aes(x = stats::reorder(source, score),
                                            y = score)) +
  ggplot2::geom_bar(mapping = ggplot2::aes(fill = score),
                    color = "black",
                    stat = "identity") +
  ggplot2::scale_fill_gradient2(low = colors[1],  # Blue (downregulated)
                                mid = "whitesmoke",
                                high = colors[2], # Red (upregulated)
                                midpoint = 0) +
  cowplot::theme_cowplot() +
  theme(
    plot.title = element_text(size = 15, face = "bold"),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 1),
    axis.text.y = element_text(angle = 45, hjust = 1), # Rotating y-axis labels
    axis.text.x = element_text(angle = 45, hjust = 1)  # Rotating x-axis labels
  ) +
  grids(linetype = "dashed") +
  ggplot2::xlab("Transcription Factors (TFs)") +
  ggplot2::ylab("Score (TF Activity)") +
  ggplot2::ggtitle("Top Transcription Factors in ARID1A Mutant vs. Non-Mutant Melanoma Patients") +
  ggplot2::labs(fill = "Activity Score")




pdf(# The directory you want to save the file in
  width = 10, # The width of the plot in inches
  height = 4,
  file = "./results/tcga_analysis/tf_activity_true_vs_false.pdf")
# Showing plot

p
dev.off()



#panel C - RFX5, RFXAP, CIITA and TWIST1

vol_in<-rbind(data.frame(deg, TF="RFX5", in_regulon=rownames(deg) %in% net$target[net$source=="RFX5"]),
              data.frame(deg, TF="RFXAP", in_regulon=rownames(deg) %in% net$target[net$source=="RFXAP"]),
              data.frame(deg, TF="CIITA", in_regulon=rownames(deg) %in% net$target[net$source=="CIITA"]),
              data.frame(deg, TF="TWIST1", in_regulon=rownames(deg) %in% net$target[net$source=="TWIST1"]),
              data.frame(deg, TF="JUN", in_regulon=rownames(deg) %in% net$target[net$source=="JUN"]),
              data.frame(deg, TF="ARID1A", in_regulon=rownames(deg) %in% net$target[net$source=="ARID1A"]))

vol_in$X<-NULL
vol_in$label<-""
vol_in$label[vol_in$in_regulon]<-rownames(vol_in)[vol_in$in_regulon]
vol_in<-vol_in[vol_in$in_regulon,]

# Define alpha values based on 'below'
alpha_values <- ifelse(!vol_in$in_regulon, 0.05, 1)

# Create labels based on conditions
vol_in$label <- ifelse((vol_in$logFC >= -0.5 & vol_in$logFC <= 0.5) | vol_in$adj.P.Val >= 0.01, "", vol_in$label)
vol_in$colour <- ifelse((vol_in$logFC >= -0.5 & vol_in$logFC <= 0.5) | vol_in$adj.P.Val >= 0.01, "darkgrey", "darkred")
library(wesanderson)
# Your existing ggplot code
ggplot(vol_in, aes(logFC, -log10(adj.P.Val), color = colour, label=label)) +
  geom_point() +
  geom_vline(xintercept = c(-0.5, 0.5), linetype = "dashed") + # Vertical lines
  geom_hline(yintercept = 1, linetype = "dashed") +        # Horizontal line
  ylab("Log10(Adjusted P value)") + xlab("Log fold change") + facet_wrap(~TF, scales = "free") +
  cowplot::theme_cowplot()+
  geom_text_repel(min.segment.length = 0, seed = 42, box.padding = 0.5, colour="black") +
  theme(legend.position="none") +
  scale_colour_manual(values = wes_palette("Royal1"))





