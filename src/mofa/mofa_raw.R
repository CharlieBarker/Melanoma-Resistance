
#find correct environment 
packLib="/usr/lib/R"
if (file.exists(packLib)) {
  reticulate::use_condaenv("py37", required = T)
  setwd(dir = "~/Desktop/Melanoma_Resistance/")
}else {
  setwd(dir = "~/Desktop/Melanoma_Resistance/")
}


library(data.table)
library(MOFA2)
library(dplyr)
library(purrr)
library(stringr)
library(ggplot2)
library(PhosR)
library(EnsDb.Hsapiens.v86)

list_of_inputs<-list(#signalome_view=data.frame(read.csv("input_data/signalome_by_sample.csv")),
  phospho=data.frame(read.csv("./data/input_data/phosphosites.csv"))
  ,protein=data.frame(read.csv("./data/input_data/proteins.csv"))
  ,mRNA=data.frame(read.csv("./data/input_data/rna_expression.csv"))
  )


list_of_inputs<-lapply(list_of_inputs, function(x){rownames(x)<-x$X; return(x)})
list_of_inputs<-lapply(list_of_inputs, function(x){x$X<-NULL; return(x)})
list_of_inputs<-lapply(list_of_inputs, function(x){colnames(x)<-sub("*\\.[0-9]", "", colnames(x)); return(x)})

#anova to find the most important TF features
aov <- matANOVA(mat=list_of_inputs$mRNA, grps=colnames(list_of_inputs$mRNA))
idx <- (aov < 0.001) & (rowSums(abs(list_of_inputs$mRNA) > 1) > 0) #idx is the filter 0.5 and aov of 0.05
list_of_inputs$mRNA <- list_of_inputs$mRNA[idx, ,drop = FALSE]

list_of_inputs$protein<-list_of_inputs$protein[rowSums(is.infinite(data.matrix(list_of_inputs$protein))) == 0,]
aov <- matANOVA(mat=list_of_inputs$protein, grps=colnames(list_of_inputs$protein))
idx <- (aov < 0.001) & (rowSums(abs(list_of_inputs$protein) > 1) > 0) #idx is the filter 0.5 and aov of 0.05
list_of_inputs$protein <- list_of_inputs$protein[idx, ,drop = FALSE]


genename_df <- AnnotationDbi::select(EnsDb.Hsapiens.v86, keys = rownames(list_of_inputs$protein), keytype = "UNIPROTID", columns = "GENENAME")
new_names<-make.unique(genename_df$GENENAME[match(rownames(list_of_inputs$protein), genename_df$UNIPROTID)])
rownames(list_of_inputs$protein)[!is.na(new_names)] <- new_names[!is.na(new_names)]

list_of_inputs_df<-lapply(list_of_inputs, function(dat){reshape2::melt(as.matrix(dat))})

#scale this?
#scale this and views
#list_of_inputs_df$proteome_view<-NULL

inputs_df<-bind_rows(list_of_inputs_df, .id = "view")
colnames(inputs_df)<-c("view", "feature", "sample", "value")
inputs_df$sample<-str_replace(string = inputs_df$sample, pattern = "vemurafenib.trametinib", replacement = "vemurafenib_and_trametinib")


#so we have to aggregate the samples together
#iwould suggest to do an annova
long_mofa_input<-aggregate(.~view+feature+sample, inputs_df, mean)

MOFAobject <- create_mofa(long_mofa_input)
data_overview_plot<-plot_data_overview(MOFAobject)

data_opts <- get_default_data_options(MOFAobject)
head(data_opts)

model_opts <- get_default_model_options(MOFAobject)
model_opts$num_factors <- 4
head(model_opts)

train_opts <- get_default_training_options(MOFAobject)
head(train_opts)
data_opts$scale_views<-T
data_opts$scale_groups<-T
MOFAobject <- prepare_mofa(
  object = MOFAobject,
  data_options = data_opts,
  model_options = model_opts,
  training_options = train_opts
)

#run MOFA
outfile = file.path("./results/mofa/mofa_object.hdf5")
MOFAobject.trained <- run_mofa(MOFAobject, outfile, use_basilisk = T)

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
                          factor = 1:4,
                          color_by = "Drug",
                          shape_by = "Genetic"
)


#WRITE PDFS

pdf(file = "results/mofa/visualisation/mofa_variance.pdf",width=6,height=10)
data_overview_plot
ggpubr::ggarrange(variance_heat, variance_per_view, nrow = 2,align = "hv") #density plots per batch
dev.off()
pdf(file = "results/mofa/visualisation/mofa_factors.pdf",width=7,height=2)
facotrs_plot
dev.off()

