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
