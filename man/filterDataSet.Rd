\name{filterDataSet}
\alias{filterDataSet}
\title{
  	This function filters our lowly expressed genes in RNAseq data. 
}
\description{
		This function filters our lowly expressed genes in RNAseq data.
}
\usage{
filterDataSet(data,filterPerc=0.75)
}
\arguments{
  \item{data}{
		A dataset with genes as rows and samples as columns.
}
  \item{filterPerc}{
		a number specifying the percent of expression values that are not equal to 0 for a gene.
}
}
\details{
	This function removes any genes in a dataset that have an expression value of 0 for a specified percentage of samles. 
}
\value{
	A data frame is returned.
}
\author{
	Jessica Mar
}
\examples{
data(exprs.dat)
exprs.filtered.dat <- filterDataSet(exprs.dat)
}
\keyword{methods}
