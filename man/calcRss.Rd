\name{calcRss}
\alias{calcRss}
\title{
		Function calculates the average RSS for a set of cluster assignments. 
}
\description{
		Function calculates the average RSS for a set of cluster assignments. 
}
\usage{
calcRss(exprs.dat, cl, class.vector)
}
\arguments{
  \item{exprs.dat}{
		a \code{matrix} of gene expression values. 
 }
  \item{cl}{
		a \code{vector} of cluster assignments.
}
  \item{class.vector}{
		a \code{vector} specifying the group membership of the samples.
}
}
\details{
	This function is called internally by \code{findSynexprs}. 
	For an informative cluster, the RSS values should be very small relative to those produced by the informativeness metric (the MSS values). 
}
\value{
		A numeric value representing the average RSS value for this set of cluster assignments.
}
\author{
	Jessica Mar
}
\examples{
\dontrun{
library(cluster)
data(subset.loring.eset)
clustObj <- agnes(as.dist(1-t(cor(exprs(subset.loring.eset)))))
crss.vals <- NULL 
for( i in 1:10 ){
	crss.vals <- c(crss.vals, calcRss(exprs(subset.loring.eset), cutree(clustObj,i), pData(subset.loring.eset)$celltype))
}
# The RSS values are expected to be smaller than the informativeness metric values in the presence of genuine cluster structure.
}
}
\keyword{methods}
