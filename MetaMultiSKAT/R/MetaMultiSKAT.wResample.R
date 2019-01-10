Meta.MultiSKAT.wResample.work <- function(study.list, method.s = "Hom",resample = 500,Sigma_g = NULL,Sigma_p = NULL,Sigma_s = NULL){
  
  with.resample = TRUE
  H1 <- Harmonize.Test.Info(study.list)
  n.studies <- length(H1)
  
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


Meta.MultiSKAT.wResample <- function(study.list, method.s = c("Hom"),method.p = c("PhC"),resample = 500,r.corr = 0,Sigma_p = NULL,Sigma_s = NULL,verbose = FALSE){
  
  #study.list; method.s = c("Hom");method.p = c("PhC");resample = 500;r.corr = 0;Sigma_p = NULL;Sigma_s = NULL;verbose = FALSE;
  n.studies <- length(study.list)
  
  if(is.null(Sigma_s)){
    if(method.s == "Hom"){
      Sigma_s = matrix(1,ncol = n.studies,nrow = n.studies)
    }else if(method.s == "Het"){
      Sigma_s = diag(n.studies)
    }else{
      stop("method.s can only be one of Hom or Het")
    }
  }
  
  Sigma_p <- assign.Sigma_p(study.list,method.p=method.p,Sigma_p = Sigma_p,verbose=verbose)
  lg <- length(r.corr); lp <- length(Sigma_p);
  
  if(lg == 1 && lp == 1){
    if(verbose)
      print(paste("Note: Using 1 Sigma_g and 1 Sigma_p to create MultiSKAT tests in each study"))
    out1 <- Meta.MultiSKAT.wResample.fixedP.VarG(study.list,Sigma_s=Sigma_s,method.s = method.s,Sigma_p = Sigma_p[[1]],r.corr = r.corr[1],resample = resample,verbose = verbose)
  }else if(lg > 1 && lp == 1){
    if(verbose)
      print(paste("Note: Using multiple Sigma_g and 1 Sigma_p to create MultiSKAT tests in each study"))
    out1 <- MetaMultiSKAT:::Meta.MultiSKAT.wResample.fixedP.VarG(study.list,Sigma_s=Sigma_s,method.s = method.s,Sigma_p = Sigma_p[[1]],r.corr = r.corr,resample = resample,verbose = verbose)
  }else if(lg > 1 && lp > 1){
    if(verbose)
      print(paste("Note: Using multiple Sigma_g and study-specific Sigma_p's to create MultiSKAT tests in each study"))
    out1 <- Meta.MultiSKAT.wResample.VarP.VarG(study.list,Sigma_s=Sigma_s,method.s =method.s,method.p =method.p,Sigma_p =Sigma_p,r.corr = r.corr,resample = resample, verbose = verbose)
  }else if(lg == 1 && lp > 1){
    if(verbose)
      print(paste("Note: Using 1 Sigma_g and study-specific Sigma_p's to create MultiSKAT tests in each study"))
    out1 <- Meta.MultiSKAT.wResample.VarP.VarG(study.list,Sigma_s=Sigma_s,method.s =method.s,method.p =method.p,Sigma_p =Sigma_p,r.corr = r.corr,resample = resample, verbose = verbose) 
  }else{
    stop("Provide correct r.corr and Sigma_p/method.p")
  }
  
  return(out1)
  
}

assign.Sigma_p <- function(list.score,method.p=c("PhC"),Sigma_p = NULL,verbose=TRUE){
  
  n.studies <- length(list.score)
  if((! class(Sigma_p)%in% c("list","matrix")) && !is.null(Sigma_p)){
    stop("Provide Sigma_p as a list or matrix")
  }
  
  if(is.null(Sigma_p)){
    
    Sigma_p <- list()
    
    if(method.p == "Het"){
      for(i in 1:n.studies){
        Sigma_p[[i]] <- diag(dim(list.score[[i]]$Score.Object$Score.Matrix)[1])
      }
    }else if(method.p == "Hom"){
      for(i in 1:n.studies){
        Sigma_p[[i]] <- matrix(1,ncol = dim(list.score[[i]]$Score.Object$Score.Matrix)[1],nrow = dim(list.score[[i]]$Score.Object$Score.Matrix)[1])
      }
    }else if(method.p == "PhC"){
      if(verbose)
        print("study.specific MultiSKAT tests are constructed")
      for(i in 1:n.studies){
        Sigma_p[[i]] <- list.score[[i]]$Score.Object$y.cov
      }
    }else{
      stop("method.p can only be one of Hom, Het, PhC")
    }
  }
  
  if(!is.null(Sigma_p)){
    if(class(Sigma_p) == "matrix"){
      Sigma_p <- list(Sigma_p)
    }else if(class(Sigma_p) == "list"){
      if(length(Sigma_p) != n.studies){
        stop("Provide a list of Sigma_P's as long as list.score")
      }
    }
  }
  
  return(Sigma_p)
}

Meta.MultiSKAT.wResample.fixedP.VarG <- function(list.score,Sigma_s=NULL,method.s = c("Hom"),Sigma_p = NULL,r.corr = 0,resample = 500,verbose = FALSE){
  
  n.studies <- length(list.score)
  
  if(is.null(Sigma_s)){
    if(method.s == "Hom"){
      Sigma_s = matrix(1,ncol = n.studies,nrow = n.studies)
    }else if(method.s == "Het"){
      Sigma_s = diag(n.studies)
    }else{
      stop("method.s can only be one of Hom or Het")
    }
  }
  
  if(verbose){
    print("Running Meta.MultiSKAT.base tests")
  }
  
  ds <- dim(list.score[[1]]$Score.Object$Score.Matrix)[1]
  if(is.null(Sigma_p)){
    Sigma_p <- diag(ds)
  }
  
  if(! is.matrix(Sigma_p))
    stop("Please provide Sigma_P in a matrix form")
  
  l.g <- length(r.corr)
  MM.obj <- list()
  for(i in 1:l.g){
    if(verbose)
      print(paste("Running",i,"th test"))
    lst <-  Transform.MultiSKAT.fixedP.fixedG(list.score,Sigma_p = Sigma_p,r.corr = r.corr[i])
    MM.obj[[i]] <- Meta.MultiSKAT.wResample.work(lst,method.s = method.s,resample = resample)
  }
  
  re <- MM.obj
  return(re)
  
}

Meta.MultiSKAT.wResample.VarP.VarG <- function(list.score,Sigma_s=NULL,method.s = c("Hom"),method.p = c("PhC"),Sigma_p = NULL,r.corr = 0,resample = 500,verbose = FALSE){
  
  n.studies <- length(list.score)
  
  ds <- dim(list.score[[1]]$Score.Object$Score.Matrix)[1]
  
  if(is.null(Sigma_p)){
    Sigma_p <- list()
    if(method.p == "Het"){
      for(i in 1:n.studies){
        Sigma_p[[i]] <- diag(dim(list.score[[i]]$Score.Object$Score.Matrix)[1])
      }
    }else if(method.p == "Hom"){
      for(i in 1:n.studies){
        Sigma_p[[i]] <- matrix(1,ncol = dim(list.score[[i]]$Score.Object$Score.Matrix)[1],nrow = dim(list.score[[i]]$Score.Object$Score.Matrix)[1])
      }
    }else if(method.p == "PhC"){
      for(i in 1:n.studies){
        Sigma_p[[i]] <- list.score[[i]]$Score.Object$y.cov
      }
    }else{
      stop("method.p can only be one of Hom, Het, PhC")
    }
    
  }
  
  if(length(Sigma_p) != n.studies){
    stop("Provide a list of Sigma_P's as long as list.score")
  }
  
  if(is.null(Sigma_s)){
    if(method.s == "Hom"){
      Sigma_s = matrix(1,ncol = n.studies,nrow = n.studies)
    }else if(method.s == "Het"){
      Sigma_s = diag(n.studies)
    }else{
      stop("method.s can only be one of Hom or Het")
    }
  }
  
  if(verbose){
    print("Running Meta.MultiSKAT.base tests")
  }
  
  l.g <- length(r.corr)
  
  MM.obj <- list()
  for(i in 1:l.g){
    if(verbose)
      print(paste("Running",i,"th test"))
    lst <-  Transform.MultiSKAT.fixedG(list.score,Sigma_p = Sigma_p,r.corr = r.corr[i])
    MM.obj[[i]] <- Meta.MultiSKAT.wResample.work(lst,method.s = method.s,resample = resample)
  }
  
  re <- MM.obj
  return(re)
  
}
