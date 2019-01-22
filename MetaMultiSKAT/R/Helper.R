mat.sqrt <- function(A)
{
  ei<-eigen(A)
  d<-ei$values
  d<-(d+abs(d))/2
  d2<-sqrt(d)
  ans<-ei$vectors %*% diag(d2) %*% t(ei$vectors)
  return(ans)
}

MAF <- function(G){
  mf <- array()
  for(i in 1:ncol(G)){
    mf[i] <- sqrt(sum(G[,i] == 0,na.rm = T)/(nrow(G)-sum(is.na(G[,i]))));
    if(mf[i] > 0.5) 
      mf[i] = 1- mf[i];
  }
  return(mf)
}

MAC <- function(G){
  mc <- array()
  for(i in 1:ncol(G)){
    mc[i] <- sum(G[,i] == 1,na.rm = T) + 2*min(sum(G[,i] == 2,na.rm = T),sum(G[,i] == 0,na.rm = T));
  }
  return(mc);
}

Get_GenInverse<-function(V.item){
	
	
	#V.item <-cov(y.mat)[1:10,1:10]
	temp<-eigen(V.item)
	idx<-which(abs(temp$values) > 10^-5)
	V.item.inv<- temp$vectors[,idx] %*% diag( 1/temp$values[idx]) %*%  t(temp$vectors[,idx])
	
	return(V.item.inv)
	
}

invnorm <- function(x){
  res = rank(x)
  res = qnorm(res/(length(res)+0.5))
  return(res)
}

Get_EffectSize<-function(MAF, c=0.4, MAF.cutoff=1, p.neg=0){
  
  beta<-abs(log10(MAF)) *c
  n.neg<-round(length(beta) * p.neg)
  
  if(n.neg > 0){
    idx<-sample(1:length(beta), n.neg)
    beta[idx]<- -beta[idx]
  }
  
  beta[MAF > MAF.cutoff] = 0
  return(beta)
}
blockMatrixDiagonal <-
function(...){  
  matrixList<-list(...)
  if(is.list(matrixList[[1]])) matrixList<-matrixList[[1]]
  
  dimensions<-sapply(matrixList,FUN=function(x) dim(x)[1])
  finalDimension<-sum(dimensions)
  finalMatrix<-matrix(0,nrow=finalDimension,ncol=finalDimension)
  index<-1
  for(k in 1:length(dimensions)){
    finalMatrix[index:(index+dimensions[k]-1),index:(index+dimensions[k]-1)]<-matrixList[[k]]
    index<-index+dimensions[k]
  }
  finalMatrix
}
