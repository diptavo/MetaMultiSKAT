Meta.MultiSKAT.base <- function(H1,Sigma_s=NULL,method.s = c("Hom")){
  
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
  Sigma_g = diag(p1); Sigma_p  = diag(n.pheno)
  
  kernel.score = Sigma_s%x%Sigma_g%x%Sigma_p
  kernel.var = Sigma_s%x%diag(ncol(Sigma_g))%x%diag(ncol(Sigma_p))
  Test.stat = t(score.vector)%*%kernel.score%*%score.vector
  L.mat <- MultiSKAT:::mat.sqrt(kernel.var)
  eigs <- eigen(t(L.mat)%*%phi%*%L.mat)
  pv <- SKAT:::Get_PValue.Lambda(Re(eigs$values),Test.stat)$p.value
  re <- list(Q = Test.stat, p.value = pv, Sigma_s = Sigma_s,method = method.s)
  class(re) <- "Meta.MultiSKAT_Object"
  return(re)
}


Meta.MultiSKAT.wResample <- function(H1, method.s = "Hom",resample = 500,with.resample = TRUE){
  
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
  Sigma_g = diag(p1); Sigma_p  = diag(n.pheno)
  
  kernel.score = Sigma_s%x%Sigma_g%x%Sigma_p
  kernel.var = Sigma_s%x%diag(ncol(Sigma_g))%x%diag(ncol(Sigma_p))
  Test.stat = t(score.vector)%*%kernel.score%*%score.vector
  L.mat <- MultiSKAT:::mat.sqrt(kernel.var)
  eigs <- eigen(t(L.mat)%*%phi%*%L.mat)
  pv <- SKAT:::Get_PValue.Lambda(Re(eigs$values),Test.stat)$p.value
  
  kp <- length(H1[[1]]$Score.Object$Score.Matrix)
  
  
  chols <- list()
  for(i in 1:n.studies)
    chols[[i]] <- mat.sqrt(H1[[i]]$Regional.Info.Pheno.Adj)
  
  p1.norm <- array();
  s.norm <- list()
  for(i in 1:n.studies){
    s.norm[[i]] <- matrix(rnorm(resample*kp),nrow = kp,ncol = resample)
  }
  
  system.time(for(id in 1:resample){
    for(js in 1:n.studies){
      s.norm[[js]][,id] <- t(chols[[js]])%*%s.norm[[js]][,id];
    }
    
    score.vector.norm <- NULL
    for(js in 1:n.studies){
      score.vector.norm <- c(score.vector.norm,s.norm[[js]][,id])
    }
    
    Test.stat.het = t(score.vector.norm)%*%kernel.score%*%score.vector.norm
    
    p1.norm[id] <- SKAT:::Get_PValue.Lambda(Re(eigs$values),Re(Test.stat.het))$p.value  
    
  })
  
  re <- list(p.value = pv,Q = as.numeric(Re(Test.stat)), n.studies = n.studies, null.p.values = p1.norm)
  class(re) <- "Meta_MultiSKAT.wResample"
  
  
  return(re)
}

minP.Meta <- function(list.MM.objects){
  
  n.test <- length(list.MM.objects)
  
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
  re <- list(Meta.MultiSKAT.Objects = list.MM.objects,p.value = minp.meta,n.test = n.test)
  class(re) <- "minP.Meta.Object"
  return(re)
}





