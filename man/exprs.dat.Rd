\name{exprs.dat}
\alias{exprs.dat}
\docType{data}
\encoding{latin1}
\title{ Gene Expression Matrix of Published Data }
\description{
  This is a matrix object containing published gene expression data from Mueller et al. (NCBI GEO 
  accession id GSE11508). The data set contains 11044 probes for 68 samples. From the original data set,
  we have selected four cell lines giving a total of 68 samples - embryonic stem cells (12 samples), 
  neural progenitors (31 samples), neural stem cells (8 samples) and teratoma-differentiated cells (17 samples). 
  The lines have also been restricted based on Illumina BeadChip platform, and only those collected using the 
  WG-6 version have been used.
  
  We also applied a quality filter to the original gene expression data where a probe was retained if it passed a 0.99
  detection score in 75% of samples in at least one of the four cell types. 
}
\usage{data(exprs.dat)}
\format{
  A matrix with normalized log2 expression intensities for 11044 probes on 68 samples (representing 4 different cell types).
}
\value{
	A matrix object containing published gene expression data from Mueller et al. (NCBI GEO 
  	accession id GSE11508). The data set contains 11044 probes for 68 samples.
}
\seealso{
\code{\link{samp.info}}, \code{\link{loring.eset}}
}
\references{
M\"{u}ller F, et al., Regulatory networks define phenotypic classes of human stem cell lines. Nature, 2008. 455(7211): p. 401-405.
}
\examples{
data(exprs.dat)
}
\keyword{datasets}
