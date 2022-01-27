quality = function(ADJtrue,ADJ){
  dims <- dim(ADJtrue)
  J <- dims[2]  
  tn <-0
  fn <-0
  tp <-0
  fp <-0
  for (j in 2:J) {
    for (j2 in 1:(j-1)) {
      if ((ADJtrue[j,j2] == 0) & (ADJ[j,j2]== 0)) tn<-tn+1
      if ((ADJtrue[j,j2] == 1) & (ADJ[j,j2]== 1)) tp<-tp+1
      if ((ADJtrue[j,j2] == 0) & (ADJ[j,j2]== 1)) fp<-fp+1
      if ((ADJtrue[j,j2] == 1) & (ADJ[j,j2]== 0)) fn<-fn+1
      }
  }
  tpr <- tp/(tp+fn)
  tnr <- tn/(tn+fp)
  g   <- sqrt(tpr*tnr)
  return(c("tp"=tp, "tn"=tn, "fp"=fp, "fn"=fn, "tpr"=tpr,"tnr"=tnr,"g"=g))
}