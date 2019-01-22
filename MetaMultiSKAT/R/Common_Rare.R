MMS.MAF.Objs <- function(study.list,maf.cutoff = 0.05,weights.beta.common = c(0.5,0.5),weights.beta.rare = c(1,25)){
  n.studies <- length(study.list)
  list1 <- list2 <- list()
  for(i in 1:n.studies){
    if(class(study.list[[i]])!="MultiSKAT.Test.Item")
      stop("Please provide a MultiSKAT.Test.Item")
    A <- Extract.MAF.Objs(study.list[[i]],maf.cutoff = maf.cutoff,weights.beta.common = weights.beta.common,weights.beta.rare = weights.beta.rare)
    list1[[i]] <- A$Obj.Common; list2[[i]] <- A$Obj.rare
  }
  return(list(common = list1,rare = list2))
}

Extract.MAF.Objs <- function(Score.object,maf.cutoff = 0.05,weights.beta.common = c(0.5,0.5),weights.beta.rare = c(1,25)){
  
  wts <- dbeta(Score.object$maf,Score.object$Score.Object$weights.beta[1], Score.object$Score.Object$weights.beta[2])
  out1 <- Score.object;
  out1$Score.Object$Score.Matrix <- Score.object$Score.Object$Score.Matrix%*%diag(1/wts); out1$Score.Object$Score.Matrix.kern <- Score.object$Score.Object$Score.Matrix.kern%*%diag(1/wts)
  out1$Regional.Info.Pheno.Adj <- diag(1/rep(wts,each =Score.object$Score.Object$n.pheno))%*%Score.object$Regional.Info.Pheno.Adj%*%diag(1/rep(wts,each = Score.object$Score.Object$n.pheno))
  out1$Method.Sigma_g <- diag(1/wts)%*%out1$Method.Sigma_g%*%diag(1/wts)
  
  ind.CR <- as.numeric(out1$maf > maf.cutoff)
  wts.common <- dbeta(out1$maf,weights.beta.common[1],weights.beta.common[2]); Sigma_g.common <- diag(wts.common)%*%out1$Method.Sigma_g%*%diag(wts.common)
  Sigma_g.common <- diag(ind.CR)%*%Sigma_g.common%*%diag(ind.CR)
  
  out.common <- Transform.MultiSKAT(out1,Sigma_g = Sigma_g.common,Sigma_p = out1$Method.Sigma_P); colnames(out.common$Score.Object$Score.Matrix) <- colnames(Score.object$Score.Object$Score.Matrix)
  out.common$weights.beta <- weights.beta.common
  
  ind.CR <- as.numeric(out1$maf <= maf.cutoff)
  wts.rare <- dbeta(out1$maf,weights.beta.rare[1],weights.beta.rare[2]); Sigma_g.rare <- diag(wts.rare)%*%out1$Method.Sigma_g%*%diag(wts.rare)
  Sigma_g.rare <- diag(ind.CR)%*%Sigma_g.rare%*%diag(ind.CR)
  
  out.rare <- Transform.MultiSKAT(out1,Sigma_g = Sigma_g.rare,Sigma_p = out1$Method.Sigma_P);colnames(out.rare$Score.Object$Score.Matrix) <- colnames(Score.object$Score.Object$Score.Matrix)
  out.rare$weights.beta <- weights.beta.rare
  
  return(list(Obj.Common = out.common,Obj.rare = out.rare))
}

Meta.MultiSKAT.Common_Rare <- function(study.list,Sigma_s=NULL,method.s = c("Hom"),maf.cutoff = 0.05,weights.beta.common = c(0.5,0.5),weights.beta.rare = c(1,25),Sigma_g = NULL,Sigma_p = NULL){
  
  CR.obj <- MMS.MAF.Objs(study.list,maf.cutoff=maf.cutoff,weights.beta.common = weights.beta.common,weights.beta.rare = weights.beta.rare)
  H2 <-  Harmonize.Test.Info(CR.obj$common); H3 <- Harmonize.Test.Info(CR.obj$rare)
  
  common <- Meta.MultiSKAT.base(H2,Sigma_s=NULL,method.s = method.s,Sigma_g = Sigma_g,Sigma_p = NULL)
  rare <- Meta.MultiSKAT.base(H3,Sigma_s=NULL,method.s = method.s,Sigma_g = Sigma_g,Sigma_p = NULL)
  wt.common <- sum(common$lambda); wt.rare <- sum(rare$lambda);
  wt.common <- wt.common/(wt.common+wt.rare); wt.rare = 1-wt.common
  Q.CR <- wt.common*common$Q + wt.rare*rare$Q
  W.CR <- wt.common*common$W + wt.rare*rare$W
  pv <- SKAT:::Get_Davies_PVal(Q.CR/2,W.CR)$p.value
  re <- list(Q = Q.CR, p.value = pv, Sigma_s = common$Sigma_s,method = method.s,W = W.CR)
  
}
