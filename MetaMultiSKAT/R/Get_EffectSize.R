Get_EffectSize <-
function(MAF, c=0.4, MAF.cutoff=1, p.neg=0){
  
  beta<-abs(log10(MAF)) *c
  n.neg<-round(length(beta) * p.neg)
  
  if(n.neg > 0){
    idx<-sample(1:length(beta), n.neg)
    beta[idx]<- -beta[idx]
  }
  
  beta[MAF > MAF.cutoff] = 0
  return(beta)
}
