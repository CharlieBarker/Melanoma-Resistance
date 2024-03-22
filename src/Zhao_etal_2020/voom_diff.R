library(tidyverse)
library(ggrepel)
library(limma)
library(edgeR)
setwd("~/MelanomaProject/A375r/data/counts/")
# from https://www.sciencedirect.com/science/article/pii/S2352340920305047
count_files<-list.files()
list_counts<-list()
for (file in count_files) {
  out<-read.delim(file, header = F)
  colnames(out)[2]<-file
  list_counts[[file]]<-data.frame(out)
}
cont_mat<-list_counts %>% reduce(inner_join, by = "V1")     
colnames(cont_mat) <- gsub("\\..*","", colnames(cont_mat))

key<-data.frame(code=c("SRR10982881", "SRR10982882", "SRR10982883", 
                       "SRR10982884", "SRR10982885", "SRR10982886"),
                condition=c("A375R", "A375R", "A375R", 
                            "A375", "A375", "A375"))

select <- order(rowMeans(cont_mat[,2:ncol(cont_mat)]), decreasing=TRUE)[1:5000]
highexprgenes_counts <- cont_mat[select,]
highexprgenes_counts<- highexprgenes_counts[,2:ncol(highexprgenes_counts)]
data_for_PCA <- t(highexprgenes_counts)
mds <- cmdscale(dist(data_for_PCA), k=2, eig=TRUE) 
Data_to_plot<- data.frame(mds$points)
Data_to_plot$Sample<- rownames(Data_to_plot)
Data_to_plot$group<- key$condition[match(Data_to_plot$Sample,key$code)]

p <- ggplot(Data_to_plot, aes(X1,X2, color=group))
p1 <- p + geom_point()+
  xlab("Principle component 1")+
  ylab ('Principle component 2')+
  theme_bw() +
  theme(axis.text.y = element_text(angle = 0, hjust = 1,size=14 ),
        axis.text.x = element_text(angle = 90, size=14, hjust=1,vjust=0.5),
        axis.title.x = element_text(size=14),
        axis.title.y = element_text(size=14),
        plot.background = element_blank() ,
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  geom_text_repel(data = Data_to_plot, aes(label = group), size= 3 , color= 'blue' , fontface='italic', show.legend = FALSE, box.padding = unit(0.2, "line"))+
  theme(legend.title=element_blank())+
  theme(plot.title = element_text(size = 14))+
  theme(legend.position= 'right')
print(p1)

mycounts <- cont_mat
rownames(mycounts) <- cont_mat$V1
mycounts$V1<-NULL
isexpr <- rowSums(cpm(mycounts)>1) >= 1
mycounts <- mycounts[isexpr,]
experiment_design<-key
rownames(experiment_design) <- experiment_design$code
group<-factor(experiment_design$condition)
design <- model.matrix(~0+group)
colnames(design)<- gsub("group","",colnames(design))
nf <- calcNormFactors(mycounts)
y <- voom(mycounts,design,lib.size=colSums(mycounts)*nf, plot = T)

counts.voom <- y$E
# save normalised expression data into output dir
Normalised_data_to_print<- as.data.frame(counts.voom)

#Fit a linear model:
fit <- lmFit(y,design)
cont.matrix <- makeContrasts(A375R - A375,levels=design)
fit <- contrasts.fit(fit, cont.matrix)
fit <- eBayes(fit)
mycoef="A375R - A375"
limma.res <- topTable(fit,coef=mycoef,n=dim(fit)[1])
limma.res.pval <- topTable(fit,coef=mycoef,n=dim(fit)[1],p.val=1)
ENSEMBL_convert<- read.table("~/MelanomaProject/Sumana_RNAseq_analysis/ENSEMBL_mapping_human.txt", header=T)

#Set the fold change level:
FCL_level<- 2.5
padj_level <- 0.01
diff_df<-limma.res.pval
diff_df$genes<- ENSEMBL_convert$symbol[match(rownames(diff_df), ENSEMBL_convert$ENSEMBL_ID)]
diff_df[which(diff_df['adj.P.Val'] > padj_level & abs(diff_df['logFC']) < FCL_level ),"group"] <- "Not_significant"
diff_df[which(diff_df['adj.P.Val'] < padj_level & abs(diff_df['logFC']) < FCL_level ),"group"] <- "Significant"
diff_df[which(diff_df['adj.P.Val'] > padj_level & abs(diff_df['logFC']) > FCL_level ),"group"] <- "FoldChange"
diff_df[which(diff_df['adj.P.Val'] < padj_level & abs(diff_df['logFC']) > FCL_level ),"group"] <- "Significant&FoldChange"
diff_df$Pos_neg<- NA
diff_df[which(diff_df['adj.P.Val'] < padj_level & diff_df['logFC'] > FCL_level & diff_df['group']=='Significant&FoldChange'),"Pos_neg"] <- "Differentially_high"
diff_df[which(diff_df['adj.P.Val'] < padj_level & diff_df['logFC'] < -FCL_level& diff_df['group']=='Significant&FoldChange'),"Pos_neg"] <- "Differentially_low"

p2 <- ggplot(diff_df,aes(logFC,-log10(adj.P.Val), color=group))+
  theme_bw()+
  geom_point()+
  theme(axis.text.y = element_text(angle = 0, hjust = 1,size=14), 
        axis.text.x = element_text(angle = 0, size=14, hjust=0.5,vjust=1), 
        axis.title.x = element_text(size=14), 
        axis.title.y = element_text(size=14),
        plot.background = element_blank() ,
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  theme(axis.line = element_line(color = 'black'))+
  scale_color_manual(values=c("gray", "gray", "gray", "#ca0020"))+
  ylab('-log(p-value)')+
  xlab('Fold change')+
  theme(legend.position="none")
print(p2)

#write.csv(diff_df, file = "../lfc.csv")

#GSEA 

library(fgsea)
pathways.hallmark <- gmtPathways("~/phenotype_networks/data/gene_set_files/h.all.v7.1.symbols.gmt") #Human_GO_AllPathways_with_GO_iea_April_01_2020_symbol.gmt

ranks<-diff_df$logFC
names(ranks)<-diff_df$genes
# Load the pathways into a named list

fgseaRes <- fgsea(pathways=pathways.hallmark, stats=ranks, nperm=1000)
fgseaRes[fgseaRes$padj<0.1,]

#DOROTHEA
library(dorothea)
library(ggplot2)
library(dplyr)
library(viper)

net <- dorothea::dorothea_hs

###subset to the threshold - keep only the most confident TFs
Regulon_file<- net[net$confidence=='A'
                            | net$confidence=='B'
                            | net$confidence=='C'
                            | net$confidence=='D'
                            ,]

# Estimatez-score values for the GES. Check VIPER manual for details
myStatistics = matrix(diff_df$logFC, dimnames = list(diff_df$genes, 'log2FC') )
myPvalue = matrix(diff_df$adj.P.Val, dimnames = list(diff_df$genes, 'adj.P.Val') )
mySignature = (qnorm(myPvalue/2, lower.tail = FALSE) * sign(myStatistics))[, 1]
mySignature = mySignature[order(mySignature, decreasing = T)]
# Estimate TF activities
mrs = msviper(ges = mySignature, regulon = df2regulon(Regulon_file), ges.filter = F, minsize = 4)


TF_activities = data.frame(Regulon = names(mrs$es$nes),
                           Size = mrs$es$size[ names(mrs$es$nes) ], 
                           NES = mrs$es$nes, 
                           p.value = mrs$es$p.value, 
                           FDR = p.adjust(mrs$es$p.value, method = 'fdr'))
TF_activities = TF_activities[ order(TF_activities$p.value), ]
TF_activities[TF_activities$FDR<0.1,]

library(progeny)
library(pheatmap)

progeny_in<-data.frame(counts.voom)
progeny_in$gene<-rownames(progeny_in)
Normalised_counts_matrix <- as.tibble(progeny_in) %>% 
  dplyr::mutate_if(~ any(is.na(.x)),~ if_else(is.na(.x),0,.x)) %>% 
  tibble::column_to_rownames(var = "gene") %>% 
  as.matrix()
colnames(Normalised_counts_matrix) <- key$condition[match(colnames(Normalised_counts_matrix), key$code)]
rownames(Normalised_counts_matrix) <- ENSEMBL_convert$symbol[match(rownames(Normalised_counts_matrix), ENSEMBL_convert$ENSEMBL_ID)]

PathwayActivity_counts <- progeny(Normalised_counts_matrix, scale=TRUE, 
                                  organism="Human", top = 100)
Activity_counts <- as.vector(PathwayActivity_counts)
paletteLength <- 100
myColor <- 
  colorRampPalette(c("darkblue", "whitesmoke","indianred"))(paletteLength)

progenyBreaks <- c(seq(min(Activity_counts), 0, 
                       length.out=ceiling(paletteLength/2) + 1),
                   seq(max(Activity_counts)/paletteLength, 
                       max(Activity_counts), 
                       length.out=floor(paletteLength/2)))

progeny_hmap <- pheatmap(t(PathwayActivity_counts),fontsize=14, 
                         fontsize_row = 10, fontsize_col = 10, 
                         color=myColor, breaks = progenyBreaks, 
                         main = "PROGENy (100)", angle_col = 45,
                         treeheight_col = 0,  border_color = NA)