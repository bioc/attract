\name{findOnepwaySynexprs}
\alias{findOnepwaySynexprs}
\title{
		Internal function - finds the synexpression groups for a single given pathway.
}
\description{
		Internal function - finds the synexpression groups for a single given pathway.
}
\usage{
findOnepwaySynexprs(pathwayIds, myAttractorModuleSet, min.clustersize = 5, removeGenes = NULL, ...)
}

\arguments{
  \item{pathwayIds}{
		a single character string denoting the KEGG ID of the pathway module to be analyzed. 
}
  \item{myAttractorModuleSet}{
		\code{AttractorModuleSet} object, output of the \code{findAttractors} step.
}
  \item{min.clustersize}{
		integer specifying the minimum number of genes that must be present in clusters that are inferred. 
}
  \item{removeGenes}{
		\code{vector} of probes that specify those genes who demonstrate little variability across the different celltypes and thus should be removed from 
		downstream analysis. 
}
  \item{\dots}{
		additional arguments.
}
}
\details{
	This function is called internally by \code{calcFuncSynexprs}. Users should use \code{calcFuncSynexprs} rather than calling \code{findOnepwaySynexprs} directly. 
}
\value{
	A \code{SynExpressionSet} object is returned. 
}
\author{
		Jessica Mar
}
\examples{
\dontrun{
data(subset.loring.eset)
attractor.states <- findAttractors(subset.loring.eset, "celltype", nperm=10, annotation="illuminaHumanv1.db")
remove.these.genes <- removeFlatGenes(subset.loring.eset, "celltype", contrasts=NULL, limma.cutoff=0.05)
map.syn <- findOnepwaySynexprs("04010", attractor.states, removeGenes=remove.these.genes)
}
}
\keyword{internal}
