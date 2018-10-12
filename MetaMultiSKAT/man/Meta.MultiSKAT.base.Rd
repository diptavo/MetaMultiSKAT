\name{Meta.MultiSKAT.base}
\alias{Meta.MultiSKAT.base}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Perform Meta MultiSKAT test with a given kernel.
}
\description{
}
\usage{
Meta.MultiSKAT.base(H1,Sigma_s=NULL,method.s = c("Hom"))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{H1}{
	Output from Harmonize.Test.Info function. Contains a harmonized list of summary statistics for each study as a list.}
  \item{Sigma_s}{
	Kernel representing the heterogeneity of effects between studies}
  \item{method.s}{
	"Hom": Homogeneous (Sigma_s) or "Het": Heterogeneous (Sigma_s)}

}
 
