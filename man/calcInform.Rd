\name{calcInform}
\alias{calcInform}

\title{
	Function calculates the informativeness metric (average MSS) for a set of cluster assignments.
}
\description{
	Function calculates the informativeness metric (average MSS) for a set of cluster assignments.
}
\usage{
calcInform(exprs.dat, cl, class.vector)
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
	This function is also called internally by \code{findSynexprs}. 
}
\value{
		A numeric value representing the average MSS value (informativeness metric) for a set of cluster assignments. 
		For an informative cluster, the RSS values should be very small relative to those produced by the informativeness metric (the MSS values). 
}
\author{
		Jessica Mar
}
\references{
Mar, J., C. Wells, and J. Quackenbush, Defining an Informativeness Metric for Clustering Gene Expression Data. to appear, 2010.
}
\examples{
\dontrun{
library(cluster)
data(subset.loring.eset)
clustObj <- agnes(as.dist(1-t(cor(exprs(subset.loring.eset)))))
cinform.vals <- NULL 
for( i in 1:10 ){
	cinform.vals <- c(cinform.vals, calcInform(exprs(subset.loring.eset), cutree(clustObj,i), pData(subset.loring.eset)$celltype))
}
k <- (1:10)[cinform.vals==max(cinform.vals)]			# gives the optimal number of clusters 
}
}
\keyword{methods}
