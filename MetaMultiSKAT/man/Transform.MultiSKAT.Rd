\name{Transform.MultiSKAT}
\alias{Transform.MultiSKAT}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Transform score objects to MultiSKAT test objects
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
	Extracted score objects from individual studies (built using Extract.Test.Info function); class should be "MultiSKAT.Test.Item"
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
List containing the following
  \item{Score.Object }{List containing scores extracted from the MultiSKAT test}
  \item{Test.Stat}{Test statistic of the MultiSKAT test}
  \item{p.value}{p-value of the MultiSKAT test}
  \item{Regional.Info.Pheno.Adj}{Phenotype adjusted region information matrix}
  \item{Method.Sigma_P}{Sigma_P kernel used}
  \item{Method.Sigma_g}{Sigma_G kernel used}
  \item{maf}{MAF(Minor allele frequency)'s of the variants in the region}
}
