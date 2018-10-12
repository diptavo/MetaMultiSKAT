\name{Transform.MultiSKAT}
\alias{Transform.MultiSKAT}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Extract kernelized scores for a Multi-SKAT test
}
\description{
Transform extracted statistics to MultiSKAT test objects with intended phenotype and genotype kernels.
}
\usage{
Transform.MultiSKAT(S1,Sigma_g, Sigma_p)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{S1}{
	Extracted score objects from individual studies (built using Extract.Test.Info function)
}
  \item{Sigma_g}{
Genotype kernel
}
  \item{Sigma_p}{
phenotype kernel
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (obj.data, Z1, Sigma_p) 
{
    R = Sigma_p
    R1 = mat.sqrt(Sigma_p)
    n.pheno <- obj.data$n.pheno
    res.mat1 = matrix(obj.data$res.V, byrow = TRUE, ncol = n.pheno)
    q1 = obj.data$q1
    N = obj.data$n.all
    V.item.inv = obj.data$V.item.inv
    Q.temp1 <- t(res.mat1) \%*\% Z1
    m = ncol(Z1)
    colnames(Q.temp1) <- colnames(Z1)
    if (Is.Common) {
        Q <- sum(colSums(Q.temp1)^2)/2
        m1 = m
    }
    if (!Is.Common) {
        tr <- matrix(0, ncol = m, nrow = n.pheno)
        for (i in 1:m) tr[, i] <- R1 \%*\% Q.temp1[, i]
        Q <- sum(tr^2)
        m1 <- m * n.pheno
    }
    resid = NULL
    re <- list(Score.Matrix = tr, n.pheno = n.pheno, y.cov = Get_GenInverse(obj.data$V.item.inv))
    class(re) <- "MetaMultiSKAT.Score.Object.kern"
    return(re)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
