#plot a 3d pca
pca_3D<-function(data, colour, out_file){
  dat_pc<-t(data.frame(na.omit(data)))
  pca <- prcomp(dat_pc)
  scores = as.data.frame(pca$x) 
  scores$colour <- as.factor(design[match(rownames(scores), paste0("X", design$code)),colour])
  palette(rainbow(length(unique(scores$colour))))
  plot(PC1~PC3, scores, col = scores$colour)
  library(rgl)
  par3d(cex=.4)
  with(scores, plot3d(PC1,PC2,PC3, col = as.numeric(scores$colour), size = 5))
  legend3d("topright", legend = levels(scores$colour), col = as.factor(scores$colour), pch=19)
  rgl.postscript(out_file,fmt="pdf")
}

