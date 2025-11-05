library(MOFA2)
library(dplyr)
library(EnsDb.Hsapiens.v86)
library(tidyr)
library(ggplot2)
library(wesanderson)
library(tidyverse)

setwd("~/Desktop/Melanoma_Resistance//")
source("./src/functions/default_variables.R")

MOFAobject.trained<-load_model(file = "./results/mofa/mofa_object.hdf5")
factors_melanoma<-data.frame(MOFA2::get_factors(MOFAobject.trained)$single_group)
# Convert to dataframe and keep rownames
factors_melanoma <- factors_melanoma %>%
  as.data.frame() %>%
  tibble::rownames_to_column("sample")

# Extract ARID1A status and drug treatment from sample names
factors_melanoma <- factors_melanoma %>%
  dplyr::mutate(
    ARID1A_status = ifelse(grepl("KO", sample), "KO", "WT"),
    Treatment = gsub("__.*", "", sample) # everything before __
  )

# Reshape to long format for ggplot
factors_melanoma_long <- factors_melanoma %>%
  pivot_longer(cols = starts_with("Factor"), names_to = "Factor", values_to = "Value")

factors_melanoma_long <- factors_melanoma_long %>%
  mutate(Treatment = dplyr::recode(Treatment,
                                   "Trametinib_10nM" = "Trametinib",
                                   "Untreated" = "Untreated",
                                   "vemurafenib_and_trametinib" = "Combination",
                                   "Vermurafenib_1uM" = "Vemurafenib"
  ))

pdf(file = "./paper/Figures/mofa_Factors.pdf",width=8,height=5.5)

# Plot
write.csv(factors_melanoma_long, file = 'paper/all_figure_data/1D.csv')
ggplot(factors_melanoma_long, aes(x = 1, y = Value, 
                    shape = ARID1A_status, color = Treatment)) +
  facet_wrap(~ Factor, scales = 'free_y') +
  geom_jitter(width = 0.1, size = 4, alpha = 0.8) +
  labs(y = "Factor value", x = "Factor", 
       shape = "ARID1A status", color = "Treatment") + 
  cowplot::theme_cowplot(font_size = 16) +   # increase overall font size
  theme(
    plot.title = element_text(size = 20, face = "bold"),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 1),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    strip.text = element_text(size = 16, face = "bold"),  # facet labels bigger
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 15)
  ) 
dev.off()

variance_explained<-data.frame(head(MOFAobject.trained@cache$variance_explained$r2_per_factor[[1]]))

# Convert rownames to a column
variance_explained <- variance_explained %>%
  rownames_to_column(var = "Factor")

# Pivot the dataframe to a long format
variance_explained_long <- variance_explained %>%
  pivot_longer(cols = -Factor, names_to = "Type", values_to = "Variance") %>%
  mutate(Type = factor(Type, levels = c("mRNA", "phospho", "protein")),
         Factor = factor(Factor, levels = c("Factor1", "Factor2", "Factor3")))

#WRITE PDFS
pdf(file = "./paper/Supplementary_plots/mofa_supp.pdf",width=8.27,height=8)

# Plot with ggplot2
ggplot(variance_explained_long, aes(x = Factor, y = Variance, fill = Type)) +
  geom_bar(stat = "identity", position = "dodge", color = "black") +
  labs(title = "Variance explained by each factor in each data modality.",
       x = "Factor",
       y = "Variance Explained (%)",
       fill = "Data Type") +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 14, face = "bold"),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 1)
  ) +
  scale_fill_manual(values = wes_palette("Royal1")) +
  geom_text(aes(label = sprintf("%.1f%%", Variance)),
            position = position_dodge(width = 0.9),
            vjust = -0.5, size = 3.5)  # Add text labels above bars

dev.off()



weights <- get_weights(MOFAobject.trained,
                       views = "all",
                       as.data.frame = TRUE
)

output_dir <- "./paper/Supplementary_tables"
write.csv(weights, file = paste0(output_dir, "/Supplementary_Table_S4.csv"), row.names = FALSE)
