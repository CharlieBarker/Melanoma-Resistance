
fcVSfc<-function(results, var1_name, var2_name, gene_format="uniprot", data_point="protein"){
  #get the relevant variables 

  results_all<-results
  var1_df<-results_all[results_all$experiment==var1_name,]
  var2_df<-results_all[results_all$experiment==var2_name,]
  if (data_point=="protein") {
    rownames(var1_df)<-var1_df$uniprot_accession
    rownames(var2_df)<-var2_df$uniprot_accession
  }
  if (data_point=="psite") {
    rownames(var1_df)<-paste0(var1_df$uniprot_accession, "_",
                              var1_df$residue, "_",
                              var1_df$site)
    rownames(var2_df)<-paste0(var2_df$uniprot_accession, "_",
                              var2_df$residue, "_",
                              var2_df$site)
    }
  
  var_1<-unlist(strsplit(x = unique(var1_df$experiment),split = "_"))
  var_2<-unlist(strsplit(x = unique(var2_df$experiment),split = "_"))
  
  
  diff_only_TFs<- dplyr:::select(var1_df, adj.P.Val, delabel, logFC, diffexpressed, uniprot_accession)
  diff_only_TFs$lfc_FC_WT_WT_combo<- var2_df$logFC[match(rownames(var1_df), rownames(var2_df))]
  diff_only_TFs$Pos_neg_Var2<- var2_df$diffexpressed[match(rownames(var1_df), rownames(var2_df))]
  diff_only_TFs$adj.P.Val_Var2<- var2_df$adj.P.Val[match(rownames(var1_df), rownames(var2_df))]
  diff_only_TFs$X<-rownames(diff_only_TFs)
  colnames(diff_only_TFs)<- c('adj.P.Val_Var1', 'genes', 'LFC_Var1', 'Pos_neg_Var1','uniprot_accession' , 
                              'LFC_Var2', 'Pos_neg_Var2', "adj.P.Val_Var2", "X")
  #diff_only_TFs<- diff_only_TFs[complete.cases(diff_only_TFs), ]
  diff_only_TFs$Pos_neg_Var1[is.na(diff_only_TFs$Pos_neg_Var1)]<-"NO"
  diff_only_TFs$Pos_neg_Var2[is.na(diff_only_TFs$Pos_neg_Var2)]<-"NO"
  
  diff_only_TFs$match<-paste0(diff_only_TFs$Pos_neg_Var1,"__", diff_only_TFs$Pos_neg_Var2)
  diff_only_TFs$label<- 'NO__NO'
  diff_only_TFs[diff_only_TFs$match!="NO__NO",'label']= 'Interested_genes'
  
  diff_only_TFs$match<- as.factor(diff_only_TFs$match)
  diff_only_TFs$label<- as.factor(diff_only_TFs$label)
  
  diff_only_TFs$alpha_value<- 1
  diff_only_TFs[diff_only_TFs$label!='NO__NO'| diff_only_TFs$match!='NO__NO', 'alpha_value']<- 1
  
  Fill <- c("UP__UP" = "#ca0020", 
            "UP__NO"= "deepskyblue2", 
            "NO__UP"= 'green', 
            "NO__DOWN" = "yellow", 
            "DOWN__NO" = "black", 
            "NO__NO"= 'grey')
  
  cols <- c("NO__NO" = "#FFFFFF00",  "Interested_genes"= "black")
  
  #pdf(paste0(mypath, "LFC-LFC_plot.pdf"), width= 5, height =4.75)
  #png(paste0(mypath, "LFC-LFC_plot.png"), width = 5, height = 4.75, units = 'in', res = 600)
  diff_only_TFs[diff_only_TFs == "not"]<-NA
  plusone<-diff_only_TFs
  p2 <- ggplot(diff_only_TFs,aes(x = LFC_Var1, y = LFC_Var2, fill= match, color=label, alpha= alpha_value))+
    geom_hline(yintercept=0, color='grey', linetype="dashed")+
    geom_vline(xintercept=0, color='grey', linetype="dashed")+
    geom_point(shape=21)+ geom_point(data = subset(diff_only_TFs, label!= 'No'), aes(LFC_Var1,  LFC_Var2 , fill= match, color=label, alpha= alpha_value), shape=21)+
    theme_bw()+
    theme(axis.text.y = element_text(angle = 0, hjust = 1,size=14),
          axis.text.x = element_text(angle = 0, size=14, hjust=0.5,vjust=1),
          axis.title.x = element_text(size=14),
          axis.title.y = element_text(size=14),
          plot.background = element_blank() ,
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()) +
    theme(axis.line = element_line(color = 'black'))+
    scale_fill_manual(values=Fill)+
    scale_color_manual(values=cols)+
    geom_text_repel(data = subset(diff_only_TFs, label!='NO__NO'), mapping=aes(label = genes), size= 4 , alpha=1, color='black', fontface='italic', show.legend = FALSE, box.padding = unit(0.25, "line"))+
    xlab(paste(var_1, collapse = " "))+
    ylab(paste(var_2, collapse = " "))+
    theme(legend.position="none") + stat_cor(method = "pearson", aes(x = LFC_Var1, y = LFC_Var2, fill= NULL, color=NULL))
  return(p2)
}