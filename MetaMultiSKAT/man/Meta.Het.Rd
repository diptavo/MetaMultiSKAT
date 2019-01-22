\name{Meta.Het}
\alias{Meta.Het}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Meta-MultiSKAT for heterogeneous effects
}
\description{
}
\usage{
Meta.Het(list.score,resample = 500,verbose = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{list.score}{
	list of MultiSKAT.Test.Item's corresponding to each study}
  \item{resample}{
	Number of resampling iterations; default is 500}
  \item{verbose}{
	If TRUE print messages for each step; default is FALSE}
}

\details{Performs minimum p-values based omnibus test for 4 Meta-MultiSKAT tests. Provides an aggregate p-value for the hypothesis (alternative) that the effects across multiple studies are heterogeneous. 
In particular aggregates p-values for the following choices of \eqn{latex}{ascii}{Sigma_P} and \eqn{latex}{ascii}{Sigma_G} : (PhC,SKAT), (PhC,Burden), (Het,SKAT), (Het,Burden)}

\value{
A list containing
\item{Meta.MultiSKAT.Objects}{A list of the individual Meta.MultiSKAT.wResample objects aggregated}
\item{p.value}{p-value for the Meta.Hom test}
\item{n.test}{number of tests performed. Same as the number of objects provided in the input list; by default 4}
}