RMSECV = function(X,LV,split){
  X = as.matrix(X)
  size = dim(X)
  I = size[1]
  J = size[2]
  NTest = ceiling(I/split); # maximum number of samples in each split
  TESTSET = 1:(NTest*split)*NA 
  TESTSET[1:I] = (1:I);
  TEST = matrix(TESTSET,nrow = split) 
  # Matrix TEST shows how the sample numbers are distributed over the splits. 
  # Some splits might have NAs if the number of samples do not perfectly fit in the splits.
  PRESS = 1:LV * 0;
  RMSECV = PRESS
  Xhat = X*0
  XX = X;
  for(lv in 1:LV){
    for(s in 1:split){
      Train = 1:I;
      Test = TEST[s,]; # Select the samples that go into Test set
      Test = Test[!is.na(Test)]
      Train = Train [! Train %in% Test] # Samples that go into Train set
      Xtrain = XX[Train,] # 
      Xtest = XX[Test,]
      
      # PCA on Train
      svd = svd(Xtrain);
      v = svd$v[,1]
      
      for (j in 1:J){
        # for each variable j; predict score of testset without using this variable
        # and then predict Xhat for test set for this variable
        # This guarantees independence between Xhat and X.
        vj = v[j]
        v_j = v[-j]
        t_test_j = as.matrix(Xtest[,-j]) %*% as.vector(v_j) / sum(v_j^2)
        Xhat[Test,j] = Xhat[Test,j] + t_test_j * vj
      } # end J
      
    } # end split
    
    E = X - Xhat;
    PRESS[lv] = sum(as.vector(E^2))
    RMSECV[lv] = sqrt(PRESS[lv]/(I*J));
    
    # Deflate XX with 1 component
    svd_all = svd(XX);
    XX = XX - svd_all$d[1] * (svd_all$u[,1]%*%t(svd_all$v[,1]))
  } # end LV
 
  return(RMSECV)
} # end function
