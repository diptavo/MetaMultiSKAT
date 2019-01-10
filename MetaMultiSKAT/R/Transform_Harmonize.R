Transform.MultiSKAT <- function(S1,Sigma_g, Sigma_p){
  
  m <- ncol(S1$Score.Object$Score.Matrix); n.pheno <- nrow(S1$Score.Object$Score.Matrix)
  
  Sc <- mat.sqrt(Sigma_p)%*%S1$Score.Object$Score.Matrix%*%mat.sqrt(Sigma_g)
  Ls <- mat.sqrt(Sigma_g%x%Sigma_p)
  S1$Regional.Info.Pheno.Adj <- t(Ls)%*%S1$Regional.Info.Pheno.Adj%*%(Ls)
  S1$Score.Object$Score.Matrix.kern <- Sc; S1$Method.Sigma_P = Sigma_p; S1$Method.Sigma_g = Sigma_g;
  
  Q <- sum(Sc^2)
  m1 <- m*n.pheno
  S1$Test.Stat <- Q; 
  S1$p.value <- SKAT:::Get_Davies_PVal(Q/2, S1$Regional.Info.Pheno.Adj, NULL)$p.value
  return(S1)
}





Harmonize.Test.Info <- function(Info.list){
  
  n.studies <- length(Info.list)
  n.pheno <- nrow(Info.list[[1]]$Score.Object$Score.Matrix)
  snp.list <- vector("list",n.studies)
  for(i in 1:n.studies){
    snp.list[[i]] <- colnames(Info.list[[i]]$Score.Object$Score.Matrix);
  }
  all.snps <- Reduce(union,snp.list)
  for(i in 1:n.studies){
    m1 <- match(all.snps,snp.list[[i]])
    Info.list[[i]]$Score.Object$Score.Matrix <- Info.list[[i]]$Score.Object$Score.Matrix[,m1];
    Info.list[[i]]$Score.Object$Score.Matrix.kern <- Info.list[[i]]$Score.Object$Score.Matrix.kern[,m1];
    
    Info.list[[i]]$Score.Object$Score.Matrix[which(is.na(Info.list[[i]]$Score.Object$Score.Matrix))] = 0
    Info.list[[i]]$Score.Object$Score.Matrix.kern[which(is.na(Info.list[[i]]$Score.Object$Score.Matrix.kern))] = 0
    
    colnames(Info.list[[i]]$Score.Object$Score.Matrix) <- all.snps
    colnames(Info.list[[i]]$Score.Object$Score.Matrix.kern) <- all.snps
    
  }
  m = length(all.snps);
  
  for(i in 1:n.studies){
    r1 = match(all.snps,snp.list[[i]])
    r2 <- r1
    for(idx in 2:n.pheno){
      r2 = c(r2,(r1 + (idx-1)*length(snp.list[[i]])))
    }
    Info.list[[i]]$Regional.Info.Pheno.Adj <- Info.list[[i]]$Regional.Info.Pheno.Adj[r2,r2]
    Info.list[[i]]$Regional.Info.Pheno.Adj[which(is.na(Info.list[[i]]$Regional.Info.Pheno.Adj))] = 0
    colnames(Info.list[[i]]$Regional.Info.Pheno.Adj) = rownames(Info.list[[i]]$Regional.Info.Pheno.Adj) = rep(all.snps,n.pheno);
    
  }
  
  return(Info.list)
  
}
