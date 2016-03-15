\name{getPwayGenes}
\alias{getPwayGenes}
\title{
  Function removes genes from the chosen pathway that demonstrate little variation across the cell types.
}
\description{
  Function removes genes from the chosen pathway that demonstrate little variation across the cell types.
}
\usage{
getPwayGenes(pathwayIds, myAttractorModuleSet, removeGenes = NULL)
}

\arguments{
  \item{pathwayIds}{
    a single character string denoting the KEGG or reactome ID of the pathway module to be analyzed.
  }
  \item{myAttractorModuleSet}{
    \code{AttractorModuleSet} object, output of the \code{findAttractors} step.
  }
  \item{removeGenes}{
  	\code{vector} of probes or gene IDs (RNAseq) that specify those genes who demonstrate little variability across the different celltypes and thus should be removed from 
		downstream analysis. 
}
}
\details{
  This function is also called internally by \code{findSynexprs}. 
}
\value{
  	A \code{vector} of gene names that have variable expression across the different cell types in a pathway.
}
\author{
  	Jessica Mar
}
\references{
Mar, J., C. Wells, and J. Quackenbush, Defining an Informativeness Metric for Clustering Gene Expression Data. to appear, 2010.
}
\examples{
\dontrun{
data(subset.loring.eset)
attractor.states <- findAttractors(subset.loring.eset, "celltype", nperm=10, annotation="illuminaHumanv1.db")
remove.these.genes <- removeFlatGenes(subset.loring.eset, "celltype", contrasts=NULL, limma.cutoff=0.05)
map.syn <- getPwayGenes("04010", attractor.states, removeGenes=remove.these.genes)}
}
\keyword{internal}