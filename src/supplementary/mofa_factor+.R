library(MOFA2)
library(dplyr)
library(EnsDb.Hsapiens.v86)


setwd("~/Desktop/Melanoma_Resistance//")

MOFAobject.trained<-load_model(file = "./results/mofa/mofa_object.hdf5")

###enrichment of the weights ####

weights <- get_weights(MOFAobject.trained,
                       views = "all",
                       as.data.frame = TRUE
)

sample_metadata <- data.frame(
  sample = samples_names(MOFAobject.trained)[[1]],
  Drug = unlist(map(str_split(samples_names(MOFAobject.trained)[[1]], pattern = "__"), 1)),
  Genetic = unlist(map(str_split(samples_names(MOFAobject.trained)[[1]], pattern = "__"), 2))
)

samples_metadata(MOFAobject.trained) <- sample_metadata


head(MOFAobject.trained@cache$variance_explained$r2_total[[1]]) # group 1
head(MOFAobject.trained@cache$variance_explained$r2_per_factor[[1]]) # group 1

variance_per_view<-plot_variance_explained(MOFAobject.trained, x="view", y="factor")
variance_heat<-plot_variance_explained(MOFAobject.trained, x="group", y="factor", plot_total = T)[[2]]

facotrs_plot<-plot_factor(MOFAobject.trained,
                          factor = 1:3,
                          color_by = "Drug",
                          shape_by = "Genetic",
                          scale = F
)

#WRITE PDFS
data_overview_plot<-plot_data_overview(MOFAobject.trained)

pdf(file = "results/mofa/visualisation/mofa_variance.pdf",width=6,height=10)
data_overview_plot
ggpubr::ggarrange(variance_heat, variance_per_view, nrow = 2,align = "hv") #density plots per batch
dev.off()
pdf(file = "results/mofa/visualisation/mofa_factors.pdf",width=5,height=2)
facotrs_plot
dev.off()
pdf(file = "results/mofa/visualisation/mofa_factors_info.pdf",width=4,height=4)
plot_factor_cor(MOFAobject.trained, method = "spearman")
plot_variance_explained(MOFAobject.trained)
plot_variance_explained(MOFAobject.trained, plot_total = T)[[2]]
dev.off()
