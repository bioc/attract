\name{findOnepwaySynexprs}
\alias{findOnepwaySynexprs}
\title{
		Internal function - finds the synexpression groups for a single given pathway.
}
\description{
		Internal function - finds the synexpression groups for a single given pathway.
}
\usage{
findOnepwaySynexprs(myIDs, myDataSet, cellTypeTag, min.clustersize = 5, removeGenes = NULL, ...)
}

\arguments{
  \item{myIDs}{
		a single character string denoting the KEGG or reactome ID of the pathway module to be analyzed or a character code{vector} of gene names of a pathway. 
}
  \item{myDataSet}{
		\code{AttractorModuleSet} object, output of the \code{findAttractors} step. This could also be an \code{ExpressionSet} object.
}
  \item{cellTypeTag}{
  	character string of the variable name which stores the cell-lineages or experimental groups of interest for the samples in the data set 
		(this string should be one of the column names of pData(myEset)).   
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
map.syn <- findOnepwaySynexprs("04010", attractor.states, "celltype", removeGenes=remove.these.genes)
vec.geneid <- c("GI_17999531-S","GI_17978503-A")
custom.syn <- findSynexprs(vec.geneid, subset.loring.eset, "celltype", removeGenes=remove.these.genes)
}
}
\keyword{internal}
