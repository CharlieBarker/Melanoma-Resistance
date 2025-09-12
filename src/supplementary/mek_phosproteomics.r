library(benchmarKIN)
library(dplyr)
library(biomaRt)

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


combo_phospho<-read.csv(file = "./results/lfc/phospho/combination_lfc.csv")
tram_phospho<-read.csv(file = "./results/lfc/phospho/trametinib_lfc.csv")
vem_phospho<-read.csv(file = "./results/lfc/phospho/vemurafenib_lfc.csv")


# Extract relevant columns
combo <- combo_phospho[, c("X", "logFC")]
tram <- tram_phospho[, c("X", "logFC")]
vem <- vem_phospho[, c("X", "logFC")]

# Rename logFC columns to condition names
colnames(combo) <- c("site", "combination")
colnames(tram) <- c("site", "trametinib")
colnames(vem) <- c("site", "vemurafenib")

# Merge datasets on 'site'
merged <- Reduce(function(x, y) merge(x, y, by = "site", all = TRUE), list(combo, tram, vem))

# Set site as rownames
rownames(merged) <- merged$site
merged$site <- NULL

# Assume merged is your full joined dataframe, with 'X' as substrate-site IDs
mat <- merged

# ✅ Set rownames from 'X' column (substrate ID)
rownames(mat) <- merged$X

# 🧹 Remove 'X' column from the data itself
mat$X <- NULL

# ✅ Convert to a pure numeric matrix
mat <- as.matrix(mat)
mode(mat) <- "numeric"
rownames(mat) <- rownames(merged)

#get activities from zscore
act_scores <- run_zscore(mat = data.frame(mat), network = kinase_net)
act_scores
wt_act_scores<-act_scores
