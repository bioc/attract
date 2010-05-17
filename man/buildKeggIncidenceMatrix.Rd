\name{buildKeggIncidenceMatrix}
\alias{buildKeggIncidenceMatrix}
\title{
		Internal function - buildings the KEGG incidence matrix requird by GSEAlm.
}
\description{
		Internal function - buildings the KEGG incidence matrix requird by GSEAlm.
}
\usage{
buildKeggIncidenceMatrix(kegg.ids, gene.ids, annotation)
}
\arguments{
  \item{kegg.ids}{
		character vector of KEGG pathway ids. 
 }
  \item{gene.ids}{
		character vector of gene ids. 
}
  \item{annotation}{
		character string that denotes which annotation package to be used, eg. illuminaHumanv1.db. 
}
}
\details{
	This function is called internally by \code{findAttractors}.
}
\value{
	A matrix object with 0 and 1 entries where 1 denotes membership of a gene in a KEGG pathway, 0 denotes non-membership. 
}
\references{
Kanehisa, M. and S. Goto, KEGG: Kyoto Encyclopedia of Genes and Genomes. . Nucleic Acids Res., 2000. 28: p. 27-30.
}
\author{
		Jessica Mar
}
\examples{
\dontrun{
# this takes a long time!
library(illuminaHumanv2.db) 
kegg.ids <- ls(illuminaHumanv2PATH2PROBE) 
gene.ids <- ls(illuminaHumanv2PATH)
imat <- buildKeggIncidence(kegg.ids, gene.ids, illuminaHumanv2.db) 
}
}
\keyword{internal}
