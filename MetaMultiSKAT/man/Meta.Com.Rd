\name{Meta.Com}
\alias{Meta.Com}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Meta-MultiSKAT for homogeneous effects
}
\description{
}
\usage{
Meta.Com(list.score,resample = 500,verbose = FALSE)
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

\details{Performs minimum p-values based omnibus test for 8 Meta-MultiSKAT tests. Provides an aggregate p-value for the hypothesis (alternative) that the effects across multiple studies are homogeneous or heterogeneous. 
In particular aggregates p-values for the following choices of \eqn{latex}{ascii}{Sigma_P} and \eqn{latex}{ascii}{Sigma_G}$ : (PhC,SKAT), (PhC,Burden), (Het,SKAT), (Het,Burden) with \eqn{latex}{ascii}{Sigma_s} being homogeneous or heterogeneous}

\value{
A list containing
\item{p.Meta.Hom}{Meta.Hom p-value for the input objects}
\item{p.Meta.Het}{Meta.Het p-value for the input objects}
\item{p.value}{Meta.Com p-value for the input objects}
\item{Meta.MultiSKAT.Objects}{Meta-MultiSKAT objects that have been aggregated}
\item{n.test}{Number of tests performed; by default 8}
}