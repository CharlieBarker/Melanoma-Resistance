library("reshape2")
library("ggplot2")
library("ggrepel")
library('dplyr')
library(reshape2)
library(edgeR)
library(limma)
library('pheatmap')
library(RColorBrewer)


#find correct environment 
packLib="/usr/lib/R"
if (file.exists(packLib)) {
  setwd(dir = "~/Desktop/Melanoma_Resistance/")
}else {
  setwd(dir = "~/Desktop/Melanoma_Resistance/")
}


All_counts<- read.csv("./data/RNAseq/data/geneCounts_fixed.csv", header=T, fill=T)
All_counts_melted<- melt(All_counts)

#basic QC plots:
#pdf("./All_counts_all.pdf", width= 8, height =8)
#p <- ggplot(All_counts_melted, aes(variable,log(value)))
#p1 <- p + geom_boxplot()+
#  xlab("Sample")+ 
#  ylab ('Normalised counts (log10)')+
#  theme_bw() + 
#  theme(axis.text.y = element_text(angle = 0, hjust = 1,size=14 ), 
#        axis.text.x = element_text(angle = 90, size=14, hjust=1,vjust=0.5), 
#        axis.title.x = element_text(size=14), 
#        axis.title.y = element_text(size=14),
#        plot.background = element_blank() ,
#        panel.grid.major = element_blank(),
#        panel.grid.minor = element_blank())+
#  theme(legend.title=element_blank())+
#  theme(plot.title = element_text(size = 14))
#print(p1)
#dev.off()

Sample_info<- read.csv("./data/RNAseq/Study_design.csv", header=T)

#remove MED12
grps <- Sample_info$New_Sample_name[match(colnames(All_counts),Sample_info$Study_ID)]
All_counts<-All_counts[,!grepl(grps, pattern = "MED12")]
Sample_info<-Sample_info[!grepl(Sample_info$New_Sample_name, pattern = "MED12"),]


select <- order(rowMeans(All_counts[,2:dim(All_counts)[2]]), decreasing=TRUE)[1:5000]
highexprgenes_counts <- All_counts[select,]
highexprgenes_counts<- highexprgenes_counts[,2:dim(All_counts)[2]]
data_for_PCA <- t(highexprgenes_counts)
mds <- cmdscale(dist(data_for_PCA), k=2, eig=TRUE) 
Data_to_plot<- data.frame(mds$points)
Data_to_plot$Sample<- rownames(Data_to_plot)
Data_to_plot$group<- Sample_info$Sample_name[match(Data_to_plot$Sample,Sample_info$Study_ID)]

# pdf("./PCA_raw_counts.pdf", width= 8, height =6)
# p <- ggplot(Data_to_plot, aes(X1,X2, color=group))
# p1 <- p + geom_point()+
#   xlab("Principle component 1")+
#   ylab ('Principle component 2')+
#   theme_bw() +
#   theme(axis.text.y = element_text(angle = 0, hjust = 1,size=14 ),
#         axis.text.x = element_text(angle = 90, size=14, hjust=1,vjust=0.5),
#         axis.title.x = element_text(size=14),
#         axis.title.y = element_text(size=14),
#         plot.background = element_blank() ,
#         panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank())+
#   geom_text_repel(data = Data_to_plot, aes(label = Sample), size= 3 , color= 'blue' , fontface='italic', show.legend = FALSE, box.padding = unit(0.2, "line"))+
#   theme(legend.title=element_blank())+
#   theme(plot.title = element_text(size = 14))+
#   theme(legend.position= 'right')
# print(p1)
# dev.off()

#Let's normalise the counts now:
mycounts <- All_counts
rownames(mycounts) <- All_counts$ENSEMBL_ID
mycounts<- mycounts[,2:dim(All_counts)[2]]
isexpr <- rowSums(cpm(mycounts)>2) >= 1
mycounts <- mycounts[isexpr,]
experiment_design<-Sample_info
rownames(experiment_design) <- experiment_design$Study_ID
group<-factor(experiment_design$Type)
design <- model.matrix(~0+group)
colnames(design)<- gsub("group","",colnames(design))
nf <- calcNormFactors(mycounts)
y <- voom(mycounts,design,lib.size=colSums(mycounts)*nf, plot = T)
# voom_norm<-data.frame(y$E)
# voom_norm$gene<-rownames(voom_norm)
# rownames(voom_norm)<-NULL
# write_csv(x =voom_norm, file = "./counts_voom_normalised.csv", )
counts.voom <- y$E


ENSEMBL_convert<- read.table("./data/ENSEMBL_mapping_human.txt", header=T)
rownames(counts.voom)<-ENSEMBL_convert$symbol[match(rownames(counts.voom), ENSEMBL_convert$ENSEMBL_ID)]
colnames(counts.voom)<-Sample_info$New_Sample_name[match(Sample_info$Study_ID, colnames(counts.voom))]
write.csv(x = data.frame(counts.voom), file ="./data/input_data/rna_expression.csv")
