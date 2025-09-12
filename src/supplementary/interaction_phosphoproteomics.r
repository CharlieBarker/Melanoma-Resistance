library(dplyr)

setwd("~/Desktop/Melanoma_Resistance//")

# Read your psp dataset
kinase_data <- read.delim("data/psp/Kinase_Substrate_Dataset", stringsAsFactors = FALSE)
kinase_data <- kinase_data[kinase_data$KIN_ORGANISM == 'human' & kinase_data$SUB_ORGANISM == 'human',]

# Remove rows with missing kinase or substrate gene or mod site
kinase_data <- subset(kinase_data, !is.na(KINASE) & !is.na(SUB_GENE) & !is.na(SUB_MOD_RSD))
all_ids <- unique(kinase_data$KIN_ACC_ID, kinase_data$SUB_ACC_ID)

# Construct target as "SUB_GENE;SUB_MOD_RSD;"
kinase_data$target <- paste0(kinase_data$SUB_GENE, ";", kinase_data$SUB_MOD_RSD, ";")

# Construct network data.frame
kinase_net <- data.frame(
  source = kinase_data$KINASE,
  target = toupper(kinase_data$target),
  mor = 1
)
kinase_net <- kinase_net %>%
  group_by(source, target) %>%
  summarise(mor = mean(mor), .groups = "drop")


combo_wt_phospho<-read.csv(file = "./results/lfc/phospho/combination_lfc.csv")
combo_arid1a_phospho<-read.csv(file = "./results/lfc/phospho/combination_arid1a_lfc.csv")


# Extract relevant columns
combo_wt <- combo_wt_phospho[, c("X", "logFC", 'adj.P.Val')]
combo_arid1a <- combo_arid1a_phospho[, c("X", "logFC", 'adj.P.Val')]

# Rename logFC columns to condition names
colnames(combo_wt) <- c("site", "combination_wt", 'wt_adj.P.Val')
colnames(combo_arid1a) <- c("site", "combination_arid1a", 'arid1a_adj.P.Val')

# Merge datasets on 'site'
merged <- Reduce(function(x, y) merge(x, y, by = "site", all = TRUE), list(combo_wt, combo_arid1a))

pdf(file = 'paper/Supplementary_plots/arid1a_vs_wt_drug_response.pdf')

# Set site as rown

# basic scatterplot
ggplot(merged, aes(x=combination_wt, y=combination_arid1a)) + 
  geom_point() + ggtitle('Comparison of phosphosite drug responses between WT and ARID1A-KO (Combination therapy)')


merged$delta_logFC <- merged$combination_arid1a - merged$combination_wt
threshold <- 1  # change threshold, e.g. 1 log2 fold
diff_sites <- merged[abs(merged$delta_logFC) > threshold, ]

sig_diff_sites <- merged %>%
  dplyr::filter((wt_adj.P.Val < 0.05 | arid1a_adj.P.Val < 0.05) & 
           abs(delta_logFC) > 1)
merged$sig_diff <- merged$site %in% sig_diff_sites$site 

ggplot(merged, aes(x = delta_logFC, y = -log10(pmin(wt_adj.P.Val, arid1a_adj.P.Val)),
                   colour = sig_diff)) +
  geom_point() +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  labs(x = "Δ logFC (ARID1A - WT)", y = "-log10(adj.P.Val)") + ggtitle('Phosphosites with altered drug response in ARID1A-KO (Combination therapy)')




home_dir<-"~/Desktop/Melanoma_Resistance/"
setwd(home_dir)

source("./src/functions/default_variables.R")

library(tidyr)
library(ggplot2)
library(limma)
library(dplyr)
library(ggpubr)
library(scatterplot3d)
library(ghibli)
library(ggrepel)
library(EnsDb.Hsapiens.v86)

####Preprocessing####


phos_abundance=data.frame(read.csv("./data/input_data/phosphosites_unadjusted.csv"))

process_columns<-function(counts, replacement_Vec){
  # Use X column to make row names
  rownames(counts) <- counts$X
  counts <- counts[, -1]  # Remove the 'X' column
  # Remove the number after the full stop in the column names
  colnames(counts) <- gsub("\\.\\d+", "", colnames(counts))
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


out<-process_columns(phos_abundance, replacement_Vec)
counts<-out$counts
drug_names <- out$drug_names
gene_names <- out$gene_names

drug <- factor(drug_names)       # "Untreated", "Vemurafenib", "Trametinib", "Combinations"
genotype <- factor(gene_names)   # "WT", "ARID1A_KO"

# Make sure Untreated and WT are the reference levels
genotype <- relevel(genotype, ref = "WT")
drug <- relevel(drug, ref = "Untreated")

design <- model.matrix(~ genotype * drug)

fit <- lmFit(counts, design)
fit <- eBayes(fit)

ARID1A_KOxdrugCombinations<-topTable(fit, coef = "genotypeARID1A_KO:drugCombinations", sort.by = "P", n = Inf)
ARID1A_KOxdrugTrametinib<-topTable(fit, coef = "genotypeARID1A_KO:drugTrametinib", sort.by = "P", n = Inf)
ARID1A_KOxdrugVemurafenib<-topTable(fit, coef = "genotypeARID1A_KO:drugVemurafenib", sort.by = "P", n = Inf)
ARID1A_KOxdrugCombinations$Gene <- rownames(ARID1A_KOxdrugCombinations)
ARID1A_KOxdrugTrametinib$Gene <- rownames(ARID1A_KOxdrugTrametinib)
ARID1A_KOxdrugVemurafenib$Gene <- rownames(ARID1A_KOxdrugVemurafenib)

plot_volcano <- function(to_plot, title='',
                         significant_parameters = list("a"=1,   #horizontal asymptote.
                                                       "b"=.1), #vertical asymptote.
                         labelling_parameters = NULL         #seperate parameters for labelling.
) {
  
  # Define the parameters
  c <- 0   # Y-intercept
  
  # Internal function to define the mirrored function
  mirrored_asymptotic_function <- function(x,
                                           a=significant_parameters$a,
                                           b=significant_parameters$b
  ) {
    c=0
    y <- a / (abs(x) - b) + c
    return(y)
  }
  
  # Identify points above and below the mirrored function
  which_significant<- mirrored_asymptotic_function(to_plot$logFC)
  to_plot$below <- -log10(to_plot$adj.P.Val) < which_significant
  to_plot$below[abs(to_plot$logFC) < significant_parameters$b] <- TRUE
  
  # Define alpha values based on 'below'
  alpha_values <- ifelse(to_plot$below, 0.1, 0.5)
  
  to_plot$label <- ""
  
  #if labelling parameters are supplied, add labells
  if(!is.null(labelling_parameters)){
    which_labelled<- mirrored_asymptotic_function(to_plot$logFC, a=labelling_parameters$a, b=labelling_parameters$b)
    to_plot$to_label <- -log10(to_plot$adj.P.Val) < which_labelled
    to_plot$below[abs(to_plot$logFC) < labelling_parameters$b] <- TRUE
    to_plot$label[!to_plot$below] <- to_plot$Gene[!to_plot$below]
  }
  
  volcano_plot <- ggplot(to_plot,
                         aes(logFC, -log10(adj.P.Val), color = below, label = label)) +
    geom_point(alpha = alpha_values) +
    cowplot::theme_cowplot() +
    geom_text_repel(min.segment.length = 0, seed = 42, box.padding = 0.5,
                    color = "black",
                    size=2) +
    geom_vline(xintercept = 0, linetype = 'dotted', col = 'darkred') +
    theme(legend.position = "none") +
    ylab("Log10(Adjusted P value)") + xlab("Log fold change") +
    geom_function(fun = mirrored_asymptotic_function,
                  colour = ghibli_palettes$YesterdayDark[4], alpha = 0.5) +
    scale_color_manual(values = c("TRUE" = "lightgrey", "FALSE" = "#A42820")) +
    ylim(0, max(-log10(to_plot$adj.P.Val))) +
    ggtitle(title) +
    grids(linetype = "dashed") +
    theme(plot.title = element_text(size = 8, face = "bold"),
          panel.border = element_rect(colour = "black", fill = NA, linewidth = 1)) # Use linewidth instead of size
  
  
  return(volcano_plot)
}

plot_volcano(ARID1A_KOxdrugCombinations, title = 'Differential Combination therapy response: ARID1A-KO vs WT', labelling_parameters = list("a"=1,  "b"=.1))
plot_volcano(ARID1A_KOxdrugTrametinib, title = 'Differential Trametinib response: ARID1A-KO vs WT')
plot_volcano(ARID1A_KOxdrugVemurafenib, title = 'Differential Vemurafenib response: ARID1A-KO vs WT')

dev.off()

