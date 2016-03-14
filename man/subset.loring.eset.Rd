\name{subset.loring.eset}
\alias{subset.loring.eset}
\docType{data}
\encoding{latin1}
\title{ An ExpressionSet Object containing published data from M?ller et al. }
\description{
  This is an \code{ExpressionSet} object containing a subset of the published data from M?ller et al. (NCBI GEO 
  accession id GSE11508). The expression data set contains 5522 probes for 68 samples. This \code{ExpressionSet} object 
  was created specifically to demonstrate the functions in this package. If you're looking for the full M?ller data set, 
  see \code{\link{loring.eset}}.

}
\usage{data(subset.loring.eset)}
\format{
  An \code{ExpressionSet} object.
}
\value{
	An \code{ExpressionSet} object containing a subset of the published data from M?ller et al. (NCBI GEO 
  	accession id GSE11508). The expression data set contains 5522 probes for 68 samples.
}
\seealso{
	\code{\link{exprs.dat}}, \code{\link{samp.info}}, \code{\link{loring.eset}}
}
\references{
M\"{u}ller, F, et al., Regulatory networks define phenotypic classes of human stem cell lines. Nature, 2008. 455(7211): p. 401-405.
}
\examples{
data(subset.loring.eset)
subset.exprs.dat <- exprs(subset.loring.eset)		# gene expression matrix
}
\keyword{datasets}
