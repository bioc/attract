\name{samp.info}
\alias{samp.info}
\docType{data}
\encoding{latin1}
\title{ samp.info Contains the Sample Information for the Mueller data set. }
\description{
  This is sample information data frame for the samples in the Mueller data set (NCBI GEO 
  accession id GSE11508). The data frame contains the cell type groups for the 68 samples. 

}
\usage{data(samp.info)}
\format{
  A data.frame object with one column of sample IDs (these are the column IDs of the exprs.dat expression matrix object) 
  and second column indicating which cell type each sample represents. 
  \describe{
    \item{\code{ChipID}}{ A \code{vector} of sample IDs. } 
    \item{\code{celltype}}{ A \code{vector} denoting the cell type a sample represents. }
  }
}
\value{
	A sample data frame for the samples in the Mueller data set (NCBI GEO 
  	accession id GSE11508). The data frame contains the cell type groups for the 68 samples.
}
\seealso{
	\code{\link{exprs.dat}}, \code{\link{loring.eset}}
}
\references{
M\"{u}ller F, et al., Regulatory networks define phenotypic classes of human stem cell lines. Nature, 2008. 455(7211): p. 401-405.
}
\examples{
data(samp.info)
}
\keyword{datasets}
