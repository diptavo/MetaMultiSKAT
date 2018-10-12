Extract.Scores <-
function(obj.data,Z1){
  n.pheno <- obj.data$n.pheno
  res.mat1 = matrix(obj.data$res.V, byrow=TRUE, ncol=n.pheno)
  Q.temp1<-t(res.mat1)%*% Z1;
  re <- list(Score.Matrix = Q.temp1, n.pheno = n.pheno, y.cov = solve(obj.data$V.item.inv))
  colnames(Q.temp1) <- colnames(Z1);
  class(re) <- "MetaMultiSKAT.Score.Object";
  return(re)
}

Extract.Test.Info <- function(obj.data,Z1,kernel = "linear.weighted", Is.Common=FALSE, weights.beta=c(1,25), weights = NULL
                              ,impute.method = "fixed",r.corr=0,   is_check_genotype=TRUE, is_dosage = FALSE, missing_cutoff=0.15, estimate_MAF=1,max_maf=1,verbose = TRUE){
  
 
  if(class(obj.data) == "MultiSKAT_NULL_unrelated"){
    Sc = Extract.Scores.kern(obj.data,Z1,kernel = kernel,Is.Common = Is.Common, weights.beta = weights.beta,weights = weights, impute.method = impute.method,r.corr=r.corr,is_check_genotype = is_check_genotype,
                             is_dosage = is_dosage, missing_cutoff = missing_cutoff, estimate_MAF = estimate_MAF,max_maf = max_maf,verbose = verbose)}else{
      Sc = Extract.kins_Scores.kern(obj.data,Z1)}
  
  n.pheno <- obj.data$n.pheno; Sigma_p <- diag(n.pheno);
  
  Z1 <- Genotype.Kernels(Z1,obj.data,kernel = kernel,Is.Common = Is.Common, weights.beta = weights.beta,weights = weights, impute.method = impute.method,r.corr=r.corr,is_check_genotype = is_check_genotype,
                         is_dosage = is_dosage, missing_cutoff = missing_cutoff, estimate_MAF = estimate_MAF,max_maf = max_maf,verbose = verbose)
  
  pv <- MultiSKAT:::MultiSKAT_base(obj.data,Z1,Sigma_p = Sigma_p)
  re <- list(Score.Object = Sc,Test.Stat = pv$Q,P.value = pv$p.value,Regional.Info.Pheno.Adj = pv$W,Method.Sigma_P = Sigma_p)
  class(re) <- "Study.Info.Item"
  return(re)
}
