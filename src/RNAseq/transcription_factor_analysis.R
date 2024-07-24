
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

jun_tfs<-df# [grepl(df$TF, pattern = "JUN"),]
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
      TRUE ~ NA_character_
    )
  )
# ggplot(jun_tfs, aes(x=TF, y=logFC, colour=logFC)) +
#   scale_colour_gradientn(colours = pal) + 
#   geom_hline(yintercept = 0, color = "black") + # Add horizontal line at y=0
#   geom_segment(aes(x=TF, xend=TF, y=0, yend=logFC), color="grey") +
#   geom_point(aes(size = -log10(padj))) +  # Adjust the size of points here
#   theme_light() +
#   theme(
#     panel.grid.major.x = element_blank(),
#     panel.border = element_rect(color = "lightgrey", fill = NA, size = 1),
#     axis.ticks.x = element_blank()
#   ) +
#   xlab("") +
#   ylab("TF activity logFC") +
#   facet_grid(genetics ~ drug_treatment)


# Prepare the data for plotting
plot_data <- jun_tfs %>%
  dplyr::select(-exp)  %>% 
  pivot_wider(names_from = genetics, values_from = c(logFC, padj), names_sep = "_") %>%
  dplyr::filter(complete.cases(.))  # Ensure that only rows with no missing values are included

plot_data$label = ""
both_sig<-plot_data$padj_ARID1A < 0.001 | plot_data$padj_WT < 0.001
plot_data$label[both_sig]<- plot_data$TF[both_sig]
# Create the scatter plot with x = y line
library(ggrepel)
ggplot(plot_data, aes(x = logFC_ARID1A, y = logFC_WT, label=label)) +
  geom_point() +  # Plot the points
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "red") +  # Add x = y line
  facet_wrap(~drug_treatment) +  # Facet by drug treatment
  labs(x = "logFC ARID1A", y = "logFC WT", title = "LogFC Comparison Between ARID1A and WT") +
  theme_minimal()  + # Use a minimal theme
  geom_text_repel()


library(dplyr)
library(ggplot2)
library(ggfortify)

wide_df <- df %>%
  pivot_wider(names_from = exp, values_from = logFC, id_cols = TF)

# Transpose the matrix: Experiments as rows, TFs as columns
data_matrix <- wide_df %>%
  dplyr::select(-TF) %>%        # Remove the TF column
  t()                    # Transpose the data

# Convert transposed data matrix to a data frame
data_matrix_df <- as.data.frame(data_matrix)
colnames(data_matrix_df) <- wide_df$TF  # Set column names as TFs
rownames(data_matrix_df) <- colnames(wide_df)[-1]  # Set row names as experiments

# Standardize the data
data_matrix_scaled <- scale(data_matrix_df)

# Perform PCA
pca_result <- prcomp(data_matrix_scaled, center = TRUE, scale. = TRUE)

# Plot with ggplot2
autoplot(pca_result, data = data.frame(Experiment = rownames(data_matrix_df)), label = TRUE, colour = 'Experiment')


# Extract PCA loadings
pca_loadings <- pca_result$rotation

# Get the loadings for PC2
pc2_loadings <- pca_loadings[, "PC1"]

# Identify TFs that explain the most variance in PC2
# Sort by absolute value to get the most contributing TFs
top_tf_pc2 <- tibble(TF = names(pc2_loadings), Loading_PC2 = pc2_loadings) %>%
  arrange(desc(abs(Loading_PC2)))

# Print the top TFs explaining the most variance in PC2
print(top_tf_pc2)