\name{minP.Meta}
\alias{minP.Meta}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
minimum p-value based omnibus test
}
\description{
}
\usage{
minP.Meta(list.MM.objects)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{list.MM.objects}{
	list of Meta.MultiSKAT.wResample objects to be aggregated.}
}

\details{Performs minimum p-values based omnibus test on the list of Meta.MultiSKAT.wResample provided.}

\value{
A list containing
\item{Meta.MultiSKAT.Objects}{A list of the individual Meta.MultiSKAT.wResample objects aggregated}
\item{p.value}{p-value for the minP.Meta test}
\item{n.test}{number of tests performed. Same as the number of objects provided in the input list}
}