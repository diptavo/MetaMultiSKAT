Extract.Scores <-  function(obj.data,Z1){
    n.pheno <- obj.data$n.pheno
    res.mat1 = matrix(obj.data$res.V, byrow=TRUE, ncol=n.pheno)
    Q.temp1<-t(res.mat1)%*% Z1;
    re <- list(Score.Matrix = Q.temp1, n.pheno = n.pheno, y.cov = solve(obj.data$V.item.inv))
    colnames(Q.temp1) <- colnames(Z1);
    class(re) <- "MetaMultiSKAT.Score.Object";
    return(re)
}

Extract.kins_Scores.kern <-  function(obj.data,Z1,kernel = "linear.weighted", Is.Common=FALSE, weights.beta=c(1,25), weights = NULL,
                                      impute.method = "fixed",r.corr=0,   is_check_genotype=TRUE, is_dosage = FALSE, missing_cutoff=0.15, estimate_MAF=1,max_maf=1,verbose = FALSE){
  
  Z1 <- Genotype.Kernels(Z1,obj.data,kernel = kernel,Is.Common = Is.Common, weights.beta = weights.beta,weights = weights, impute.method = impute.method,r.corr=r.corr,is_check_genotype = is_check_genotype,
                         is_dosage = is_dosage, missing_cutoff = missing_cutoff, estimate_MAF = estimate_MAF,max_maf = max_maf,verbose = verbose)
  
    Sigma_p = diag(obj.res$n.pheno)    
    R = Sigma_p
    R1 = mat.sqrt(Sigma_p)
    m = ncol(Z1)
    n = nrow(Z1)
    n.pheno = obj.res$n.pheno
    q1 = obj.res$q1
    N = obj.res$n.all
    V.item.inv = obj.res$V.item.inv
    
    res.mat = matrix(obj.res$kins.adj.res.V, byrow=TRUE, ncol=n.pheno)
    Q.temp1<-t(res.mat) %*% Z1;
    
    if(Is.Common){
      Q<-sum(colSums(Q.temp1)^2)/2;
      m1=m;
    }
    if(!Is.Common){
      tr <- matrix(0,ncol = m, nrow = n.pheno);
      #R <- Sigma_P;
      #       R = cov(matrix(obj.res$res.V,ncol = obj.res$n.pheno,byrow=T))
      for(i in 1:m)
        tr[,i] <- R1%*%Q.temp1[,i];
      Q <- sum(tr^2);
      m1 <- m*n.pheno
    }
    colnames(tr) <- colnames(Z1)
    resid = NULL;
    re <- list(Score.Matrix = tr, n.pheno = n.pheno, y.cov = Get_GenInverse(obj.res$V.item.inv))  
    class(re) <- "MetaMultiSKAT.Score.Object";
    return(re)
}

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


Extract.Test.Info <- function(obj.data,Z1,kernel = "linear.weighted", Is.Common=FALSE, weights.beta=c(1,25), weights = NULL
                              ,impute.method = "fixed",r.corr=0,   is_check_genotype=TRUE, is_dosage = FALSE, missing_cutoff=0.15, estimate_MAF=1,max_maf=1,verbose = FALSE){
  
  
  
  if(class(obj.data) == "MultiSKAT_NULL_unrelated"){
    Sc = Extract.Scores.kern(obj.data,Z1,kernel = kernel,Is.Common = Is.Common, weights.beta = weights.beta,weights = weights, impute.method = impute.method,r.corr=r.corr,
                             is_check_genotype = is_check_genotype, is_dosage = is_dosage, missing_cutoff = missing_cutoff, estimate_MAF = estimate_MAF,max_maf = max_maf,verbose = verbose)
    }else if(class(obj.data) == "MultiSKAT_NULL_related"){
    Sc = Extract.kins_Scores.kern(obj.data,Z1,kernel = kernel,Is.Common = Is.Common, weights.beta = weights.beta,weights = weights, impute.method = impute.method,r.corr=r.corr,
                            is_check_genotype = is_check_genotype, is_dosage = is_dosage, missing_cutoff = missing_cutoff, estimate_MAF = estimate_MAF,max_maf = max_maf,verbose = verbose)}
  
  n.pheno <- obj.data$n.pheno; Sigma_p <- diag(n.pheno);

  maf <- MultiSKAT:::MAF(Z1); wts <- dbeta(maf,weights.beta[1],weights.beta[2])
 
  Sigma_g <-  diag(wts)%*%((1-r.corr)*diag(dim(Z1)[2]) + r.corr*matrix(1,ncol = dim(Z1)[2],nrow = dim(Z1)[2]))%*%diag(wts)

  Z1 <- Genotype.Kernels(Z1,obj.data,kernel = kernel,Is.Common = Is.Common, weights.beta = weights.beta,weights = weights, impute.method = impute.method,r.corr=r.corr,is_check_genotype = is_check_genotype,
                         is_dosage = is_dosage, missing_cutoff = missing_cutoff, estimate_MAF = estimate_MAF,max_maf = max_maf,verbose = verbose)
  
  pv <- MultiSKAT:::MultiSKAT_base(obj.data,Z1,Sigma_p = Sigma_p)
  re <- list(Score.Object = Sc,Test.Stat = pv$Q,p.value = pv$p.value,Regional.Info.Pheno.Adj = pv$W,Method.Sigma_P = Sigma_p,Method.Sigma_g = Sigma_g,maf = maf)
  class(re) <- "MultiSKAT.Test.Item"
  return(re)
}
