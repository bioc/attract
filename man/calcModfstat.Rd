\name{calcModfstat}
\alias{calcModfstat}

\title{
		Function calculates a modified F-statistic for a set of cluster assignments.
}
\description{
		Function calculates a modified F-statistic for a set of cluster assignments.
}
\usage{
calcModfstat(exprs.dat, cl, class.vector)
}

\arguments{
  \item{exprs.dat}{
		a \code{matrix} of gene expression values.
}
  \item{cl}{
		a \code{vector} of cluster assignments. 
}
  \item{class.vector}{
		a \code{vector} specifying group membership of the samples.
}
}
\details{
		This function is called internally by \code{findSynexprs}. 
}
\value{
		a modified F-statistic (average MSS/average RSS) value for a set of cluster assignments.
}
\author{
		Jessica Mar
}
\examples{
\dontrun{
library(cluster)
data(subset.loring.eset)
clustObj <- agnes(as.dist(1-t(cor(exprs(subset.loring.eset)))))
cfmod.vals <- NULL 
for( i in 1:10 ){
	cfmod.vals <- c(cfmod.vals, calcModfstat(exprs(subset.loring.eset), cutree(clustObj,i), pData(subset.loring.eset)$celltype))
}
k <- (1:10)[cfmod.vals==max(cfmod.vals)]			
}
}
\keyword{methods}
