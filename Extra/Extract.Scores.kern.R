Extract.Scores.kern <- function(obj.data,Z1,kernel = "linear.weighted", Is.Common=FALSE, weights.beta=c(1,25), weights = NULL
                                 , impute.method = "fixed",r.corr=0,   is_check_genotype=TRUE, is_dosage = FALSE, missing_cutoff=0.15, estimate_MAF=1,max_maf=1,verbose = TRUE){
  
  Z1 <- Genotype.Kernels(Z1,obj.data,kernel = kernel,Is.Common = Is.Common, weights.beta = weights.beta,weights = weights, impute.method = impute.method,r.corr=r.corr,is_check_genotype = is_check_genotype,
                         is_dosage = is_dosage, missing_cutoff = missing_cutoff, estimate_MAF = estimate_MAF,max_maf = max_maf,verbose = verbose)
  
  n.pheno <- obj.data$n.pheno
  R = diag(n.pheno)
  R1 = mat.sqrt(R)
  res.mat1 = matrix(obj.data$res.V, byrow=TRUE, ncol=n.pheno)
  q1 = obj.data$q1
  N = obj.data$n.all
  V.item.inv = obj.data$V.item.inv
  
  Q.temp1<-t(res.mat1)%*% Z1;
  m = ncol(Z1)
  
  colnames(Q.temp1) <- colnames(Z1);
  if (Is.Common) {
    Q <- sum(colSums(Q.temp1)^2)/2
    m1 = m
  }
  if (!Is.Common) {
    tr <- matrix(0, ncol = m, nrow = n.pheno)
    for (i in 1:m) tr[, i] <- R1 %*% Q.temp1[, i]
    Q <- sum(tr^2)
    m1 <- m * n.pheno
  }
  resid = NULL
  re <- list(Score.Matrix = Q.temp1, Score.Matrix.kern = tr, n.pheno = n.pheno, y.cov = Get_GenInverse(obj.data$V.item.inv),weights.beta = weights.beta)
  class(re) <- "MetaMultiSKAT.Score.Object";
  return(re)
}
