
jun_tfs<-df[grepl(df$exp, pattern = "Vermurafenib_1uM_WT_vs_") | 
              grepl(df$exp, pattern = "Trametinib_10nM_WT_vs_"),]
pal <- wes_palette("Zissou1", 100, type = "continuous")

# Add 'genetics' and 'drug_treatment' columns
jun_tfs <- jun_tfs %>%
  mutate(
    drug_treatment = case_when(
      grepl("Trametinib_10nM_WT_vs_vemurafenib_and_trametinib_WT", exp) ~ "Trametinib vs Combination therapy",
      grepl("Vermurafenib_1uM_WT_vs_vemurafenib_and_trametinib_WT", exp) ~ "Vermurafenib vs Combination therapy",
      TRUE ~ NA_character_
    )
  )

# Prepare the data for plotting
plot_data <- jun_tfs %>%
  dplyr::select(-exp)  %>% 
  pivot_wider(names_from = drug_treatment, values_from = c(logFC, padj), names_sep = "_") %>%
  dplyr::filter(complete.cases(.))  # Ensure that only rows with no missing values are included



plot_data$label = ""
both_sig<-plot_data$`padj_Trametinib vs Combination therapy` < 0.1 & plot_data$`padj_Vermurafenib vs Combination therapy` < 0.1
plot_data$label[both_sig]<- plot_data$TF[both_sig]
# Create the scatter plot with x = y line
library(ggrepel)

# Create the linear model
model <- lm(`logFC_Trametinib vs Combination therapy` ~ `logFC_Vermurafenib vs Combination therapy`, 
            data = plot_data)
# Extract residuals
plot_data$residuals <- resid(model)

tf_plot<-ggplot(plot_data, aes(x = `logFC_Trametinib vs Combination therapy`, 
                               y = `logFC_Vermurafenib vs Combination therapy`, 
                               label=label, colour = residuals)) +
  geom_point() +  # Plot the points
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "red") +  # Add x = y line
  labs(x = "logFC_Trametinib vs Combination therapy",
       y = "logFC_Vermurafenib vs Combination therapy", 
       title = "LogFC Comparison of TF activation between combination drug therapies") +
  theme_minimal()  + # Use a minimal theme
  geom_text_repel(colour="black")


pdf(file = paste0("~/Desktop/Melanoma_Resistance/paper/plots/tf_activity_combo.pdf"), 
    width = 12, height = 6)
tf_plot
dev.off()

