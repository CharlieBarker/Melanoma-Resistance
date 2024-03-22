# Define the drugs and their corresponding colors
drug_colours <- c(
  "Untreated" = "#A2AEBB",
  "Vemurafenib" = "#FFBA08",
  "Trametinib" = "#D00000",
  "Combinations" = "#3F88C5"
)
gene_colours <- c(
  "WT" = "#608569",
  "ARID1A_KO" = "#FE5F55"
)


#replacement vector
replacement_Vec<-c("WT","Untreated","Vermurafenib_1uM","Trametinib_10nM","vemurafenib_and_trametinib", "vemurafenib.trametinib", "vemurafenib+trametinib")
names(replacement_Vec)<- c("Untreated", "Untreated", "Vemurafenib", "Trametinib", "Combinations", "Combinations", "Combinations")
