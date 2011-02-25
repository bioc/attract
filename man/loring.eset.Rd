\name{loring.eset}
\alias{loring.eset}
\docType{data}
\encoding{latin1}
\title{ An ExpressionSet Object containing published data from Müller et al. }
\description{
  This is an \code{ExpressionSet} object containing the published data from Müller et al. (NCBI GEO 
  accession id GSE11508). The expression data set contains 11044 probes for 68 samples. 

}
\usage{data(loring.eset)}
\format{
  An \code{ExpressionSet} object.
}
\seealso{
	\code{\link{exprs.dat}}, \code{\link{samp.info}}
}
\references{
Müller, F, et al., Regulatory networks define phenotypic classes of human stem cell lines. Nature, 2008. 455(7211): p. 401-405.
}
\examples{
data(loring.eset)
exprs.dat <- exprs(loring.eset)		# gene expression matrix
}
\keyword{datasets}
