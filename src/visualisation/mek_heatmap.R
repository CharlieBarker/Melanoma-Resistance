library(dplyr)
library(ComplexHeatmap)
library(readxl)
library(tidyr)
library(colorRamp2)
library(wesanderson)
library(ggpubr)

setwd("~/Desktop/Melanoma_Resistance//")

data<-read.csv('data/sumana/Plate_measurements_ARID1A_luminex.csv')
psites <- read_xlsx('data/sumana/media-7.xlsx',skip = 3)


# Compute baseline (mean across replicates for Cells in no serum)
ref_means <- colMeans(data[data$Gene == "Cells in no serum", c("ERK", "MEK")])

# Add LFC columns for ERK and MEK relative to serum-free baseline
data$LFC_ERK <- log2(data$ERK / ref_means["ERK"])
data$LFC_MEK <- log2(data$MEK / ref_means["MEK"])

summary_df <- data %>%
  dplyr::filter(Gene %in% c("Control", "ARID1A", "MED12")) %>% 
  group_by(Gene, Control) %>%
  summarise(LFC_ERK = mean(LFC_ERK),
            LFC_MEK = mean(LFC_MEK), .groups = "drop")

heatmap_df <- summary_df %>%
  pivot_longer(cols = starts_with("LFC"),
               names_to = "Protein",
               values_to = "LFC") %>%
  mutate(Protein = ifelse(Protein == "LFC_ERK", "pERK", "pMEK"),
         Column = paste(Gene, Control, sep = "_")) %>%
  dplyr::select(Protein, Column, LFC) %>%
  pivot_wider(names_from = Column, values_from = LFC)

heatmap_mat <- as.matrix(heatmap_df[,-1])
rownames(heatmap_mat) <- heatmap_df$Protein

# Make a clean mapping table: Protein -> Label with phospho site
psite_labels <- psites %>%
  mutate(Label = paste(`Protein Name`, `Phospho Residue`, sep = ";")) %>%
  dplyr::select(`Protein Name`, Label)

# Match your data rows (pERK, pMEK) to the mapping
# First, create a lookup for ERK and MEK only
lookup <- c("pERK" = "ERK1/2", "pMEK" = "MEK1")

# Replace rownames with psite labels
rownames(heatmap_mat) <- sapply(rownames(heatmap_mat), function(x) {
  prot <- lookup[[x]]
  label <- psite_labels$Label[psite_labels$`Protein Name` == prot]
  if(length(label) == 0) return(x) else return(label)
})

# Make a tidy annotation dataframe from colnames
col_annot <- data.frame(
  Sample = colnames(heatmap_mat),
  stringsAsFactors = FALSE
)

# Split Gene and Treatment
col_annot <- col_annot %>%
  tidyr::separate(Sample, into = c("Gene", "Treatment"), sep = "_")
col_annot$Sample = colnames(heatmap_mat)
# Rename Gene labels
col_annot$Gene <- dplyr::recode(col_annot$Gene,
                                "Control" = "Empty vector",
                                "ARID1A" = "ARID1A-KO",
                                "MED12"  = "MED12-KO")

# Convert Treatment into Trametinib +/-
col_annot$Trametinib <- ifelse(col_annot$Treatment == "T", "+", "-")


# Desired column group order
desired_order <- c("Empty vector", "ARID1A-KO", "MED12-KO")
desired_order_treatment <- c("-", "+")

# Reorder heatmap_mat columns
col_annot$Gene <- factor(col_annot$Gene, levels = desired_order)
col_annot$Trametinib <- factor(col_annot$Trametinib, levels = desired_order_treatment)

col_annot <- col_annot[order(col_annot$Gene), ]
col_annot <- col_annot[order(col_annot$Trametinib), ]

# Reorder the matrix to group by Gene
heatmap_mat <- heatmap_mat[, col_annot$Sample]

# Create colour scale centered at 0
col_fun <- colorRamp2(c(min(heatmap_mat), 0, max(heatmap_mat)),
                      c("blue", "white", "red"))

# Create a column annotation showing Trametinib treatment as text
ha <- HeatmapAnnotation(
  `Trametinib treatment` = anno_text(
    col_annot$Trametinib,      # text to display
    rot = 1,
    location = 0.5,             # center text in cell
    just = "center",            # justify
    gp = gpar(fontsize = 20)    # text appearance
  )
)

pdf(file = 'paper/Figures/luminx_mek.pdf', width = 7, height = 1.5)
Heatmap(heatmap_mat,
        name = "log2 Fold Change vs Serum-free",
        cluster_rows = FALSE,
        cluster_columns = FALSE,
        top_annotation = ha,
        column_split = col_annot$Gene,  # splits by perturbation
        na_col = "darkgrey",
        rect_gp = gpar(col = "white", lwd = 12),
        column_gap = unit(.02, "npc"),
        row_title_gp = gpar(fontsize = 5),
        column_title_gp = gpar(fontsize = 10),
        col = col_fun,
        show_column_names = FALSE,          # remove column names
        heatmap_legend_param = list(direction = "horizontal"))  # legend at bottom


dev.off()

#or alternatively 

library(dplyr)
library(tidyr)
library(ggplot2)

# Compute baseline for normalization (Cells in no serum, untreated)
ref_means <- data %>%
  filter(Gene == "Cells in no serum", Control == "UT") %>%
  summarise(
    ERK_ref = mean(ERK),
    P38_ref = mean(P38.MAPK),
    MEK_ref = mean(MEK)
  )

# Add log2 fold change vs serum-free
data <- data %>%
  mutate(
    LFC_ERK = log2(ERK / ref_means$ERK_ref),
    LFC_P38 = log2(P38.MAPK / ref_means$P38_ref),
    LFC_MEK = log2(MEK / ref_means$MEK_ref)
  )

# Pivot longer for ggplot
plot_data <- data %>%
  select(Gene, Replicate, Control, LFC_ERK, LFC_P38, LFC_MEK) %>%
  pivot_longer(cols = starts_with("LFC"),
               names_to = "Protein",
               values_to = "LFC") %>%
  mutate(Protein = case_when(
    Protein == "LFC_ERK" ~ "ERK1/2;T202/Y204",
    Protein == "LFC_P38" ~ "P38 MAPK",
    Protein == "LFC_MEK" ~ "MEK1;S218/222"
  ))

plot_data <- plot_data %>%
  mutate(KO = case_when(
    Gene == 'ARID1A' ~ "ARID1A-KO",
    Gene == 'MED12' ~ "MED12-KO",
    Gene == 'Control' ~ "Empty vector",
    TRUE ~ NA_character_
  ))

plot_data<- plot_data[plot_data$Protein != 'P38 MAPK' & plot_data$Gene != 'Cells in no serum',]

# Reorder heatmap_mat columns
desired_order_treatment <- c("UT", "T")
plot_data$Control <- factor(plot_data$Control, levels = desired_order_treatment)
plot_data <- plot_data[order(plot_data$Control), ]

# Add sample size per group for labeling
n_table <- plot_data %>%
  group_by(Gene, Control, Protein) %>%
  summarise(n = n(), .groups = "drop")

plot_data <- plot_data %>%
  left_join(n_table, by = c("KO", "Control", "Protein"))

# Compute mean per Protein x Control x KO
means_df <- plot_data %>%
  group_by(Protein, Control, KO) %>%
  summarise(mean_LFC = mean(LFC), .groups = "drop")

# Ensure Control is a factor in the correct order (UT before T)
means_df$Control <- factor(means_df$Control, levels = c("UT", "T"))

# Base plot with points and SD
ggplot(plot_data, aes(x = Control, y = LFC, color = KO, group = KO)) +
  geom_jitter(width = 0.1, size = 4) +                       # individual points
  stat_summary(fun = mean, geom = "point", shape = 18, size = 4) +   # mean points
  stat_summary(fun.data = mean_sdl, geom = "errorbar", width = 0.2) + # SD bars
  # Add line connecting mean LFC per KO
  geom_line(data = means_df, aes(x = Control, y = mean_LFC, group = KO, color = KO),
            size = 2, linetype = "solid") +
  facet_wrap(~ Protein) +
  theme_bw() +
  labs(
    x = "Treatment",
    y = "log2 Fold Change vs Serum-free",
    title = "Phospho-protein levels per gene perturbation",
    color = "Genetic background"
  ) +
  cowplot::theme_cowplot() +
  theme(
    plot.title = element_text(size = 15, face = "bold"),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 1),
    axis.text.x = element_text(angle = 90, hjust = 1) # Rotate x-axis labels
  ) +
  grids(linetype = "dashed") +
    scale_colour_manual(values = wes_palette("BottleRocket2"))