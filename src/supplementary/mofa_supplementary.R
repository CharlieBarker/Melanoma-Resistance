library(MOFA2)
library(dplyr)
library(EnsDb.Hsapiens.v86)
library(tidyr)
library(ggplot2)
library(wesanderson)
library(tidyverse)

setwd("~/Desktop/Melanoma_Resistance//")

MOFAobject.trained<-load_model(file = "./results/mofa/mofa_object.hdf5")

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
