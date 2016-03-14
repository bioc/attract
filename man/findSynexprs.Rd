\name{findSynexprs}
\alias{findSynexprs}
\title{
		This function finds the synexpression groups present within a core attractor pathway module. 
}
\description{
		This function takes the modules that were inferred from the GSEA step using (\code{findAttractors}) and finds a set of transcriptionally  
		coherent set of genes associated with a particular core attractor pathway, i.e. the synexpression groups.
}
\usage{
findSynexprs(myIDs, myDataSet, cellTypeTag, removeGenes = NULL, min.clustersize = 5, ...)
}
\arguments{
  \item{myIDs}{
		either a single character string or \code{vector} of character strings denoting the KEGG or reactome IDs of the pathway modules to be analyzed. 
    It may also be a character code{vector} of gene names of a pathway if defining a custom pathway.
}
  \item{myDataSet}{
		\code{AttractorModuleSet} object, output of the \code{findAttractors} step. This could also be an \code{ExpressionSet} object
    if using a custom pathway.
}
  \item{cellTypeTag}{
  	character string of the variable name which stores the cell-lineages or experimental groups of interest for the samples in the data set 
		(this string should be one of the column names of pData(myEset)).
}
  \item{removeGenes}{
		\code{vector} of gene names that specify those genes who demonstrate little variability across the different celltypes and thus should be removed from 
		downstream analysis. 
}
  \item{min.clustersize}{
		\code{integer} specifying the minimum number of genes that must be present in clusters that are inferred. 
}
  \item{\dots}{
		additional arguments. 
}
}
\details{
	This function performs a hierarichical cluster analysis of the genes in a core attractor pathway module, and uses an informativeness metric 
	to determine the number of optimal clusters (syenxpression groups) that describe the data. 
}
\value{
	If a single KEGG or reactome ID is specified in \code{pwayIds}, then a \code{SynExpressionSet} object is returned. 
	If a multiple KEGG or reactome IDs are specified, then an environment object is returned where the keys are labeled "pwayIDsynexprs" 
	(e.g. for MAPK KEGGID = 04010, the key is pway04010synexprs). The value associated with each key is a \code{\link{SynExpressionSet}} object.  
}
\author{
	Jessica Mar
}
\references{
Mar, J., C. Wells, and J. Quackenbush, Identifying the Gene Expression Modules that Represent the Drivers of Kauffman's Attractor Landscape. to appear, 2010.
}
\examples{
data(subset.loring.eset)
attractor.states <- findAttractors(subset.loring.eset, "celltype", annotation="illuminaHumanv1.db")
remove.these.genes <- removeFlatGenes(subset.loring.eset, "celltype", contrasts=NULL, limma.cutoff=0.05)
mapk.syn <- findSynexprs("04010", attractor.states, "celltype", remove.these.genes)
top5.syn <- findSynexprs(attractor.states@rankedPathways[1:5,1], attractor.states, "celltype", removeGenes=remove.these.genes)
vec.geneid <- c("GI_17999531-S","GI_17978503-A")
custom.syn <- findSynexprs(vec.geneid, subset.loring.eset, "celltype", removeGenes=remove.these.genes)
}
\keyword{methods}
