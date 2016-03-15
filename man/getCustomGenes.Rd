\name{getCustomGenes}
\alias{getCustomGenes}
\title{
  Function removes genes from the custom pathway that demonstrate little variation across the cell types.
}
\description{
  Function removes genes from the custom pathway that demonstrate little variation across the cell types.
}
\usage{
getCustomGenes(vec.geneid, removeGenes = NULL)
}

\arguments{
  \item{vec.geneid}{
  	a \code{vector} of character strings denoting a custom gene set.
  }
  \item{removeGenes}{
    \code{vector} of gene names that specify those genes who demonstrate little variability across the different cellt ypes and thus should be removed from
    downstream analysis.
  }
}
\details{
  This function is also called internally by \code{findSynexprs}. 
}
\value{
  	A \code{vector} of gene names that have variable expression across the different cell types.
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
remove.these.genes <- removeFlatGenes(subset.loring.eset, "celltype", contrasts=NULL, limma.cutoff=0.05)
vec.geneid <- c("GI_17999531-S","GI_17978503-A")
customGenes <- getCustomGenes(vec.geneid, removeGenes=NULL)
}
}
\keyword{internal}