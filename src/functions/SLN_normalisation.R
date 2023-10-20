
#prep phosR data
normSLN<-function(x){
  # first basic normalization is to adjust each TMT experiment to equal signal per channel
  # figure out the global scaling value
  target <- mean(colSums(x, na.rm = T))
  # do the sample loading normalization before the IRS normalization
  # there is a different correction factor for each column
  norm_facs <- target / colSums(x, na.rm = T)
  data_sl <- sweep(x, 2, norm_facs, FUN = "*")
  return(data_sl)
}