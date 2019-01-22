Meta.MultiSKAT.base <- function(study.list,Sigma_s=NULL,method.s = c("Hom"),Sigma_g = NULL,Sigma_p = NULL){
  
  H1 <- Harmonize.Test.Info(study.list)
  n.studies <- length(H1)
  
  score.vector <- NULL
  l1 <- list()
  for(i in 1:n.studies){
    score.vector <- c(score.vector,H1[[i]]$Score.Object$Score.Matrix.kern)
    l1[[i]] <- H1[[i]]$Regional.Info.Pheno.Adj
  }
  
  p1 <- dim(H1[[1]]$Score.Object$Score.Matrix)[2]
  n.pheno <-  dim(H1[[1]]$Score.Object$Score.Matrix)[1]
  
  n.studies <- length(H1)
  if(method.s == "Hom"){
    Sigma_s = matrix(1,ncol = n.studies,nrow = n.studies)
  }else if(method.s == "Het"){
    Sigma_s = diag(n.studies)
  }
  
  
  phi <- blockMatrixDiagonal(l1)
  if(is.null(Sigma_g))
    Sigma_g = diag(p1); 
  
  if(is.null(Sigma_p))
    Sigma_p  = diag(n.pheno)
  
  kernel.score = Sigma_s%x%Sigma_g%x%Sigma_p
  kernel.var = Sigma_s%x%diag(ncol(Sigma_g))%x%diag(ncol(Sigma_p))
  Test.stat = t(score.vector)%*%kernel.score%*%score.vector
  L.mat <- MultiSKAT:::mat.sqrt(kernel.var)
  Vmat <- t(L.mat)%*%phi%*%L.mat
  eigs <- eigen(Vmat)
  pv <- SKAT:::Get_PValue.Lambda(Re(eigs$values),Test.stat)$p.value
  re <- list(Q = Test.stat, p.value = pv, Sigma_s = Sigma_s,method = method.s,lambda = Re(eigs$values),W = Vmat)
  class(re) <- "Meta.MultiSKAT_Object"
  return(re)
}



minP.Meta <- function(list.MM.objects){
  
  n.test <- length(list.MM.objects)
  
  if(n.test == 1){
    re <- list(Meta.MultiSKAT.Objects = list.MM.objects,p.value = list.MM.objects[[1]]$p.value,n.test = 1)
}else{
  for(i in 1:n.test)
    if(class(list.MM.objects[[i]]) != "Meta_MultiSKAT.wResample")
      stop("Wrong Class of objects")
  
  pv <- NULL; pv.null <- NULL
  for(i in 1:n.test){
    pv <- c(pv,list.MM.objects[[i]]$p.value)
    pv.null <- cbind(pv.null,list.MM.objects[[i]]$null.p.values)
  }
  
  minP <- min(pv,na.rm = T)
  R2<- cor(pv.null,method = "kendall");
  obj.t.Copula = tCopula(R2[lower.tri(R2)],dim = nrow(R2),dispstr = "un");
  minp.meta <- 1-pCopula(1-minP,obj.t.Copula)
  if(minp.meta < minP){
    minp.meta <- n.test*minP
  }else if(minp.meta > n.test*minP){
    minp.meta <- n.test*minP
  }
  re <- list(Meta.MultiSKAT.Objects = list.MM.objects,p.value = minp.meta,n.test = n.test)}
  class(re) <- "minP.Meta.Object"
  return(re)
}



Transform.MultiSKAT.fixedP.fixedG <- function(list.score,Sigma_p=NULL,r.corr=0){
  
  ds <- dim(list.score[[1]]$Score.Object$Score.Matrix)[1]
  if(is.null(Sigma_p)){
    Sigma_p <- diag(ds)
  }
  
  l1 <- length(list.score)
  for(i in 1:l1){
    d1 <- dim(list.score[[i]]$Score.Object$Score.Matrix)[2]
    Sigma_g <- (1-r.corr)*diag(d1) + (r.corr)*matrix(1,ncol = d1,nrow = d1)
    list.score[[i]] <- Transform.MultiSKAT(list.score[[i]],Sigma_g,Sigma_p)
  }
  return(list.score)
}

Transform.MultiSKAT.fixedG <- function(list.score,Sigma_p=NULL,r.corr=0){
  
  n.studies <- length(list.score)
  
  ds <- dim(list.score[[1]]$Score.Object$Score.Matrix)[1]
  if(is.null(Sigma_p)){
    Sigma_p <- list()
    for(i in 1:n.studies){
      ds <- dim(list.score[[i]]$Score.Object$Score.Matrix)[1]
      Sigma_p <- diag(ds)
    }
  }else if(length(Sigma_p) != n.studies){
    stop("Provide a list of Sigma_P's as long as list.score")
  }
  
  l1 <- length(list.score)
  for(i in 1:l1){
    d1 <- dim(list.score[[i]]$Score.Object$Score.Matrix)[2]
    Sigma_g <- (1-r.corr)*diag(d1) + (r.corr)*matrix(1,ncol = d1,nrow = d1)
    list.score[[i]] <- Transform.MultiSKAT(list.score[[i]],Sigma_g,Sigma_p[[i]])
  }
  return(list.score)
}





Meta.default <- function(list.score,Sigma_s=NULL,method.s = c("Hom"),method.p = c("PhC"),r.corr = 0,resample = 500,verbose = FALSE){
  
  n.studies <- length(list.score)
  
  l.p <- length(method.p)
  lst.vec <- NULL; lst <- list()
  for(i in 1:l.p){
    sgp <- as.character(method.p[i])
    lst[[i]] <- Meta.MultiSKAT.wResample.VarP.VarG(list.score,r.corr = r.corr,resample = resample,verbose=verbose,method.p = sgp,Sigma_s=Sigma_s,method.s = method.s)
    lst.vec <- c(lst.vec,lst[[i]])
  }
  
  out <- minP.Meta(lst.vec)
  return(out)
}


Meta.Hom <- function(list.score,resample = 500,verbose = FALSE){
  
  A1 <- MetaMultiSKAT:::Meta.default(list.score,method.s = c("Hom"),method.p = c("PhC","Het"),r.corr = c(0,1),verbose = verbose,resample = resample)
  
  return(A1)
}

Meta.Het <- function(list.score,resample = 500,verbose = FALSE){
  
  A1 <- MetaMultiSKAT:::Meta.default(list.score,method.s = c("Het"),method.p = c("PhC","Het"),r.corr = c(0,1),verbose = verbose,resample = resample)
  
  return(A1)
}

Meta.Com <- function(list.score,resample = 500,verbose = FALSE){
  A1 <- MetaMultiSKAT:::Meta.default(list.score,method.s = c("Het"),method.p = c("PhC","Het"),r.corr = c(0,1),verbose = verbose,resample = resample)
  A2 <- MetaMultiSKAT:::Meta.default(list.score,method.s = c("Hom"),method.p = c("PhC","Het"),r.corr = c(0,1),verbose = verbose,resample = resample)
  
  out <- minP.Meta(c(A1$Meta.MultiSKAT.Objects,A2$Meta.MultiSKAT.Objects))
  return(c(out,"p.Meta.Hom" = A2$p.value,"p.Meta.Het" = A1$p.value))
}
