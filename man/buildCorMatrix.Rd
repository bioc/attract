\name{buildCorMatrix}
\alias{buildCorMatrix}
\title{
Internal function - builds the correlation matrix between an average transcriptional module expression profile and a set of other genes.
}
\description{
Internal function - builds the correlation matrix between an average transcriptional module expression profile and a set of other genes.
}
\usage{
buildCorMatrix(dat.fr, module.genes, cor.cutoff)
}

\arguments{
  \item{dat.fr}{
		a \code{matrix} object of gene expression values. 
 }
  \item{module.genes}{
		character vector specifying genes that belong in this pathway module.
}
  \item{cor.cutoff}{
		numeric value specifying the correlation cut-off.
}
}
\details{
	This function is called internally by \code{findCorrPartners} which is easier for the user to call since \code{findCorrPartners} uses the SynExpressionSet and ExpressionSet class objects directly. 
}
\value{
		A character vector of genes that meet the correlation cut-off.
}
\author{
		Jessica Mar 
}
\examples{
\dontrun{
data(subset.loring.eset)
attractor.states <- findAttractors(subset.loring.eset, "celltype", nperm=10, annotation="illuminaHumanv1.db")
remove.these.genes <- removeFlatGenes(subset.loring.eset, "celltype", contrasts=NULL, limma.cutoff=0.05)
mapk.syn <- findSynexprs("04010", attractor.states, remove.these.genes)
cormat <- buildCorMatrix(exprs(subset.loring.eset), mapk.syn, 0.95)
}
}
\keyword{internal}
