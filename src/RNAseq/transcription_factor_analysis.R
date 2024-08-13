
#find correct environment 
packLib="/usr/lib/R"
if (file.exists(packLib)) {
  reticulate::use_condaenv("py37", required = T)
  setwd(dir = "~/Desktop/Melanoma_Resistance/")
}else {
  setwd(dir = "~/Desktop/Melanoma_Resistance/")
}

# Load necessary libraries
library(ggplot2)
library(dplyr)
library(tidyr)
library(tibble)  # Ensure tibble is loaded for rownames_to_column
library(wesanderson)
library(tidyverse)

# Define directory with the files
tf_activity_dir <- '/Users/charliebarker/Desktop/Melanoma_Resistance/results/tf_activity/'

# List all CSV files in the directory
tf_files <- list.files(tf_activity_dir, pattern = '*_tf_acts.csv', full.names = TRUE)

# Function to create volcano plot
create_volcano_plot <- function(acts_file, pvals_file, output_dir) {
  file_name <- basename(acts_file)
  
  # Remove the suffix "_tf_acts.csv" to get the experiment name
  experiment_name <- sub("_tf_acts.csv$", "", file_name)
  
  # Read the data
  tf_acts <- read.csv(acts_file, row.names = 1)
  tf_pvals <- read.csv(pvals_file, row.names = 1)
  # Reset row names as a column for both data frames
  tf_acts <- tf_acts %>%
    rownames_to_column(var = "TF") 
  colnames(tf_acts)<-c("TF", "logFC")
  tf_pvals <- tf_pvals %>%
    rownames_to_column(var = "TF") 
  colnames(tf_pvals)<-c("TF", "padj")
  
  # Merge the data
  df <- inner_join(tf_acts, tf_pvals, by = "TF")

  df$exp<-experiment_name
  return(df)
}

# Output directory for plots
output_dir <- tf_activity_dir

vol_plot_list<-list()
# Loop over each file and create plots
for (acts_file in tf_files) {
  # Construct the corresponding p-values file name
  pvals_file <- gsub("_tf_acts.csv", "_tf_pval.csv", acts_file)
  
  # Check if the corresponding p-values file exists
  if (file.exists(pvals_file)) {
    # Create volcano plot
    vol_plot_list[[acts_file]]<-create_volcano_plot(acts_file, pvals_file, output_dir)
  } else {
    warning(paste("No corresponding p-values file found for", acts_file))
  }
}

df <- do.call("rbind", vol_plot_list)
rownames(df)<-NULL

jun_tfs<-df[grepl(df$exp, pattern = "Untreated_WT_vs_") | 
              grepl(df$exp, pattern = "Untreated_ARID1A_KO_vs_"),]
pal <- wes_palette("Zissou1", 100, type = "continuous")

# Add 'genetics' and 'drug_treatment' columns
jun_tfs <- jun_tfs %>%
  mutate(
    genetics = case_when(
      grepl("ARID1A_KO", exp) ~ "ARID1A",
      grepl("WT", exp) ~ "WT",
      TRUE ~ NA_character_
    ),
    drug_treatment = case_when(
      grepl("Trametinib", exp) ~ "Untreated vs Trametinib",
      grepl("vemurafenib_and_trametinib", exp) ~ "Untreated vs Combination therapy",
      grepl("Vermurafenib", exp) ~ "Untreated vs Vemurafenib",
      grepl("Untreated_WT_vs_Untreated_ARID1A_KO", exp) ~ "Untreated vs ARID1A KO",
      TRUE ~ NA_character_
    )
  )


# Prepare the data for plotting
plot_data <- jun_tfs %>%
  dplyr::select(-exp)  %>% 
  pivot_wider(names_from = genetics, values_from = c(logFC, padj), names_sep = "_") %>%
  dplyr::filter(complete.cases(.))  # Ensure that only rows with no missing values are included


myc_tfs<-jun_tfs[jun_tfs$TF=="MYC",]

# Create a ggplot
myc_plot <- ggplot(myc_tfs, aes(x = drug_treatment, y = logFC, fill = genetics)) +
  geom_bar(stat = "identity", position = position_dodge(), width = 0.7) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  theme_minimal() +
  labs(title = "Log Fold Change of MYC Gene Under Different Conditions",
       x = "Drug Treatment",
       y = "logFC",
       fill = "Genetics") +
  theme_cowplot() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
pdf(file = paste0("~/Desktop/Melanoma_Resistance/paper/plots/tf_activity_myc.pdf"), 
    width = 12, height = 6)
myc_plot
dev.off()

plot_data$label = ""
both_sig<-plot_data$padj_ARID1A < 0.001 | plot_data$padj_WT < 0.001
plot_data$label[both_sig]<- plot_data$TF[both_sig]
# Create the scatter plot with x = y line
library(ggrepel)

# Create the linear model
model <- lm(logFC_WT ~ logFC_ARID1A, data = plot_data)
# Extract residuals
plot_data$residuals <- resid(model)

tf_plot<-ggplot(plot_data, aes(x = logFC_ARID1A, y = logFC_WT, label=label, colour = residuals)) +
  geom_point() +  # Plot the points
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "red") +  # Add x = y line
  facet_wrap(~drug_treatment) +  # Facet by drug treatment
  labs(x = "logFC ARID1A", y = "logFC WT", title = "LogFC Comparison of TF activation Between ARID1A and WT") +
  theme_minimal()  + # Use a minimal theme
  geom_text_repel(colour="black")

plot_residuals_tf <- plot_data[plot_data$padj_ARID1A < 0.1 | plot_data$padj_WT < 0.1,] %>%
  dplyr::select(-logFC_ARID1A, -logFC_WT, -padj_ARID1A, -padj_WT, -label)  %>% 
  pivot_wider(names_from = drug_treatment, values_from = residuals, names_sep = "_") %>%
  dplyr::filter(complete.cases(.))  # Ensure that only rows with no missing values are included

residual_plot_tf<-ggplot(plot_residuals_tf, aes(x = `Untreated vs Trametinib`, y = `Untreated vs Vemurafenib`, label=TF)) +
  geom_point() +  # Plot the points
  geom_vline(xintercept = 0, linetype = "dashed", color = "red") + 
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +  
  labs(x = "ARID1A ~ WT Residuals for Trametinib", y = "ARID1A ~ WT Residuals for Vermurafenib", title = "TF activity Comparison Between ARID1A and WT") +
  theme_minimal() + # Use a minimal theme
  geom_text_repel(colour="black")

#visualise the interesting transcription factors 

tfs_of_interest<-read.csv(file = "./results/collectri/tfs_of_interest.csv")

# Define the file paths
file_paths <- list.files(path = "./results/transcriptomics/", pattern = "*.csv", full.names = T)

# Create a function to extract information from the filename
extract_info <- function(filename) {
  parts <- strsplit(filename, "_vs_|_lfc.csv")[[1]]
  list(
    control = parts[1],
    test = parts[2],
    genetic_ko = ifelse(grepl("ARID1A", parts[1]), "ARID1A_KO", "WT"),
    drug = gsub("_ARID1A_KO|_WT", "", parts[2])
  )
}

# Initialize an empty list to store data frames
dfs <- list()

# Loop through each file and read the data
for (file_path in file_paths) {
  data <- read_csv(file_path, show_col_types = F)
  info <- extract_info(file_path)
  data <- data %>%
    mutate(
      control = basename(info$control),
      test = info$test,
      genetic_ko = info$genetic_ko,
      drug = info$drug
    )
  dfs[[file_path]] <- data
}

# Combine all data frames into one
combined_df <- bind_rows(dfs)
#look at HLAs
hla_df<-combined_df[grep(combined_df$gene_symbol, pattern="HLA"),]
hla_df$diffexpressed <- "Not significant"
hla_df$label <- ""

hla_df$diffexpressed[hla_df$padj < 0.1 & hla_df$log2FoldChange < 0] <- "Downregulated"
hla_df$diffexpressed[hla_df$padj < 0.1 & hla_df$log2FoldChange > 0] <- "Upregulated"
hla_df$label[hla_df$diffexpressed != "Not significant"] <- hla_df$gene_symbol[hla_df$diffexpressed != "Not significant"]

hla_df$control <- str_remove_all(hla_df$control, pattern = "_WT")
hla_df$control <- str_remove_all(hla_df$control, pattern = "_ARID1A_KO")

hla_df$test <- str_remove_all(hla_df$test, pattern = "_WT")
hla_df$test <- str_remove_all(hla_df$test, pattern = "_ARID1A_KO")

ggplot(data = hla_df, aes(x = log2FoldChange, y = -log10(padj), col = diffexpressed, label = label)) + 
  geom_vline(xintercept = c(-0.6, 0.6), col = "gray", linetype = 'dashed') + 
  geom_hline(yintercept = -log10(0.05), col = "gray", linetype = 'dashed') + 
  geom_point(size = 2) + 
  scale_color_manual(values = c("#00AFBB", "grey", "#bb0c00"), # to set the colours of our variable 
                     labels = c("Downregulated", "Not significant", "Upregulated")) + # to set the labels in case we want to overwrite the categories from the dataframe (UP, DOWN, NO) 
  labs(color = 'Severe', #legend_title, 
       x = expression("log"[2]*"FC"), y = expression("-log"[10]*"p-value")) + 
  scale_x_continuous(breaks = seq(-10, 10, 2)) + # to customise the breaks in the x axis 
  ggtitle('Expression of HLA proteins across different conditions') + # Plot title 
  geom_text_repel(max.overlaps = Inf) + # To show all labels
  facet_grid(genetic_ko~control+test)


plot_expression <- combined_df %>%
  dplyr::select(-baseMean, -lfcSE, -stat, -pvalue)  %>% 
  pivot_wider(names_from = genetic_ko, values_from = c(log2FoldChange, padj), names_sep = "_") %>%
  dplyr::filter(complete.cases(.))  # Ensure that only rows with no missing values are included
# Create the linear model
model <- lm(log2FoldChange_ARID1A_KO ~ log2FoldChange_WT, data = plot_expression)
# Extract residuals
plot_expression$residuals <- resid(model)

# Join the data frames to add the transcription factor (tf) information
plot_expression <- plot_expression %>%
  left_join(tfs_of_interest, by = c("gene_symbol" = "target")) 
plot_expression<-plot_expression[!is.na(plot_expression$source),]

expression_plot<-ggplot(plot_expression, aes(x = log2FoldChange_ARID1A_KO, y = log2FoldChange_WT, label=gene_symbol, colour =weight)) +
  geom_point() +  # Plot the points
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "red") +  # Add x = y line
  facet_grid(source~drug) +  # Facet by drug treatment
  labs(x = "logFC ARID1A", y = "logFC WT", title = "LogFC Comparison of expression Between ARID1A and WT") +
  theme_minimal() + # Use a minimal theme
  geom_text_repel(colour="black")


plot_residuals <- plot_expression[plot_expression$padj_ARID1A_KO < 0.1 | plot_expression$padj_WT < 0.1,] %>%
  dplyr::select(-log2FoldChange_ARID1A_KO, -log2FoldChange_WT, -padj_ARID1A_KO, -padj_WT, -PMID, -weight)  %>% 
  pivot_wider(names_from = drug, values_from = residuals, names_sep = "_") %>%
  dplyr::filter(complete.cases(.))  # Ensure that only rows with no missing values are included

residual_plot<-ggplot(plot_residuals, aes(x = Trametinib_10nM, y = Vermurafenib_1uM, label=gene_symbol)) +
  geom_point() +  # Plot the points
  geom_vline(xintercept = 0, linetype = "dashed", color = "red") + 
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +  
  facet_grid(~source) +  # Facet by drug treatment
  labs(x = "ARID1A ~ WT Residuals for Trametinib", y = "ARID1A ~ WT Residuals for Vermurafenib", title = "LogFC Comparison Between ARID1A and WT") +
  theme_minimal() + # Use a minimal theme
  geom_text_repel(colour="black")



pdf(file = paste0("~/Desktop/Melanoma_Resistance/paper/plots/tf_activity_deepdive.pdf"), 
    width = 12, height = 6)
tf_plot
expression_plot
residual_plot
residual_plot_tf
dev.off()

