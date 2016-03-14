\name{findCorrPartners}
\alias{findCorrPartners}
\title{
		Determines Genes with Highly Correlated Expression Profiles to a Synexpression Group
}
\description{
		This function finds genes with expression profiles highly correlated to a synexpression group. 
}
\usage{
findCorrPartners(mySynExpressionSet, myEset, removeGenes = NULL, cor.cutoff = 0.85, ...)
}
\arguments{
  \item{mySynExpressionSet}{
		\code{SynExpressionSet} object.
}
  \item{myEset}{
		\code{ExpressionSet} object.
}
  \item{removeGenes}{
		\code{vector} of probes that specify those genes who demonstrate little variability across the different celltypes and thus should be removed from 
		downstream analysis. 
}
  \item{cor.cutoff}{
		numeric value specifying the correlation cut-off.
}
  \item{\dots}{
		additional arguments. 
}
}
\details{
		Genes with highly correlated profiles to the synexpression groups (e.g. R > 0.85) are also likely to be integral in 
		maintaining cell type-specific differences, however due to their lack of inclusion in resources like KEGG, would not 
		have been picked up by the first GSEA step using \code{findAttractors}. 
}
\value{
		A \code{SynExpressionSet} object which stores the genes that are highly correlated with the synexpression group provided, and their average expression profile.
}
\author{
		Jessica Mar
}
\examples{
data(subset.loring.eset)
attractor.states <- findAttractors(subset.loring.eset, "celltype", annotation="illuminaHumanv1.db")
remove.these.genes <- removeFlatGenes(subset.loring.eset, "celltype", contrasts=NULL, limma.cutoff=0.05)
mapk.syn <- findSynexprs("04010", attractor.states, remove.these.genes)
mapk.cor <- findCorrPartners(mapk.syn, subset.loring.eset, remove.these.genes)
}
\keyword{methods}
