\name{formatrankpathways}
\alias{formatrankpathways}

\title{
		Internal function - ranks the pathways of the same P-value from GSEAlm according to a secondary filter based on size (number of genes in the pway).
}
\description{
		Internal function - ranks the pathways of the same P-value from GSEAlm according to a secondary filter based on size (number of genes in the pway).}
\usage{
formatrankpathways(gseaResMat, incidenceMat, direction = "upper")
}

\arguments{
  \item{gseaResMat}{
		a \code{matrix} object, output of \code{gsealmPerm}.
}
  \item{incidenceMat}{
		a \code{matrix}, the incidence matrix derived from KEGG used as input for GSEAlm.
}
  \item{direction}{
		a character string indicating the direction of the enrichment, either "upper" or "lower". 
}
}
\details{
	This function is called internally by the \code{\link{findAttractors}} function. Users should not need to call this function directly.
}
\value{
	A \code{matrix} object.
}
\author{
		Jessica Mar
}
\examples{
\dontrun{
	data(subset.loring.eset)
	library(illuminaHumanv2.db)
	kegg.incidence.matrix <- buildKeggIncidenceMatrix(ls(illuminaHumanv2PATH2PROBE), rownames(exprs(subset.loring.eset)), illuminaHumanv2)
	keep.pways <- apply(kegg.incidence.matrix, 1, sum) >= 5
	kegg.incidence.matrix <- kegg.incidence.matrix[keep.pways,]
	p.kegg <- gsealmPerm(subset.loring.eset, as.formula(paste("~ ", cellTypeTag, sep="")), kegg.incidence.matrix, nperm = 100, na.rm=TRUE)
	rankpmat <- formatrankpathways(p.kegg, kegg.incidence.matrix, direction="upper") 
}
}
\keyword{internal}

