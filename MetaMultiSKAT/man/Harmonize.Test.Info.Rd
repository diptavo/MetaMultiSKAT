\name{Harmonize.Test.Info}
\alias{Harmonize.Test.Info}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Harmonized Test information
}
\description{
Harmonize the summary statistic (Extract.Test.Info outputs) from individual studies.
}
\usage{
Harmonize.Test.Info(Info.list)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Info.list}{
	list of summary statistic objects (Extract.Test.Info outputs or transformed) from individual studies.
}
 
}
\value{
A list of harmonized objects as long as the input list containing the objects which have been harmonized by the variants.
}
\details{
The function harmonizes the objects for variants. The desired place for variant names is the column names of the score matrices in the individual score objects. 
Please make sure the genotype matrix had variant-ids/names as column names when running/extracting the objects.}