\name{Meta.MultiSKAT.wResample}
\alias{Meta.MultiSKAT.wResample}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Perform Meta MultiSKAT test with a given kernel with resample p-values for future omnibus test
}
\description{
}
\usage{
Meta.MultiSKAT.wResample(H1,Sigma_s=NULL,method.s = c("Hom"),resample = 500)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{H1}{
	Output from Harmonize.Test.Info function. Contains a harmonized list of summary statistics for each study as a list.}
  \item{Sigma_s}{
	Kernel representing the heterogeneity of effects between studies}
  \item{method.s}{
	"Hom": Homogeneous (Sigma_s) or "Het": Heterogeneous (Sigma_s)}
\item{resample}{
	Number of resampling iterations.}
}
 
\value{
A list containing
\item{p.value}{p.value of the Meta.MultiSKAT test}
\item{Q}{Test statistic of the test}
\item{n.studies}{number of studies over which the meta analysis was performed}
\item{null.p.values}{null p-values generated through resampling iterations. To be used in minP.Meta tests}
}