MultiSKAT_NULL_allNA <-
function (y.mat, X, is.fast = TRUE, method = "nlminb") 
{
  
  n <- nrow(y.mat)
  n.pheno <- ncol(y.mat)
  
  mis.pheno <- which(colSums(is.na(y.mat)) == dim(y.mat)[1])
  non.mis.pheno <- setdiff(1:n.pheno,mis.pheno)
  
  Corr <- cor(y.mat,use = "pairwise.complete.obs")
  Cov <- cov(y.mat,use = "pairwise.complete.obs")
  require(nlme)
  X <- cbind(X)
  
  q <- ncol(X)
  if (sum(scale(X[, 1], scale = FALSE)^2) != 0) {
    stop("The first column of X should be an intercept!")
  }
  y <- as.vector(t(y.mat))
  X1 <- X %x% diag(n.pheno)
  subject.id <- rep(1:n, each = n.pheno)
  if (is.fast == FALSE) 
    g1 = glsControl(maxIter = 50, msMaxIter = 200, tolerance = 1e-04, 
                    msTol = 1e-04, msVerbose = T, singular.ok = T, returnObject = T, 
                    opt = method)
  if (sum(c(Corr[lower.tri(cor(y.mat))]) > 1e-07,na.rm = T) == 
      0) {
    V.item = cov(y.mat)
    V.item.inv <- solve(V.item)
    gls1 <- lm(y ~ X1 - 1,na.action = na.exclude)
    res <- residuals(gls1)
    res[which(is.na(res))] <- 0
  }
  else if (is.fast == TRUE) {
    gls1 <- lm(y ~ X1 - 1,na.action = na.exclude)
    res <- residuals(gls1)
    V.item <- cov(matrix(res, ncol = n.pheno, byrow = TRUE),use = "pairwise.complete.obs")
  }
  else {
    gls1 <- gls(y ~ X1 - 1, correlation = corSymm(form = ~1 | 
                                                    subject.id), control = g1,na.action = na.exclude)
    V.item <- getVarCov(gls1)
    
    res <- residuals(gls1)
  }
  res[which(is.na(res))] <- 0
  V.item[mis.pheno,] <- 0; V.item[,mis.pheno] <- 0
  V.item.inv <- Get_GenInverse(V.item)
  q1 <- ncol(X1)
  XVX1 <- matrix(rep(0, q1 * q1), ncol = q1)
  res.V <- rep(0, n * n.pheno)
  for (i in 1:n) {
    id <- ((i - 1) * n.pheno + 1):(i * n.pheno)
    XVX1 <- XVX1 + t(X1[id, ]) %*% V.item.inv %*% X1[id,]
    res.V[id] <- V.item.inv %*% res[id]
  }
  
  XVX1_inv <- Get_GenInverse(XVX1)
  L <- princomp(y.mat,na.action = na.exclude,covmat = V.item)$loadings
  res.mat <- matrix(res,byrow = T,nrow = n,ncol = n.pheno)
  res.mat[,mis.pheno] <- 0
  pc.mat <- res.mat %*% L
  V_p = cov(pc.mat,use = "pairwise.complete.obs")
  re <- list(res = res, res.V = res.V, V.item.inv = V.item.inv, 
             cov = X, X1 = X1, XVX1_inv = XVX1_inv, n.pheno = n.pheno, 
             q1 = q1, n.all = n, id_include = 1:n, y.Corr = Corr, 
             y.Cov = Cov, pc.Cov = V_p, L = L, V.item = V.item)
  class(re) <- "MultiSKAT_NULL_unrelated"
  return(re)
}
