\name{removeFlatGenes}
\alias{removeFlatGenes}
\title{
		Flags a set of genes which demonstrates little variation across the celltypes or experimental groups of interest.
}
\description{
		This function uses a linear model set up in \code{limma} to assess the degree of association between celltype and a gene's expression 
		profile. In this way, we can flag those genes whose profiles show very little change across different celltype groups, or in other words 
		are "flat". 
}
\usage{
removeFlatGenes(eSet, cellTypeTag, contrasts = NULL, limma.cutoff = 0.05, ...)
}
\arguments{
  \item{eSet}{
		\code{ExpressionSet} object.
}
  \item{cellTypeTag}{
		character string of the variable name which stores the cell-lineages or experimental groups of interest for the samples in the data set 
		(this string should be one of the column names of pData(myEset)).   
}
  \item{contrasts}{
		optional vector of contrasts that specify the comparisons of interest. By default, all comparisons between the differnt groups are generated. 
}
  \item{limma.cutoff}{
		numeric specifying the P-value cutoff. Genes with P-values greater than this value are considered "flat" and will be included in the 
		set of flat genes. 
}
  \item{\dots}{
		additional arguments.
}
}
\details{
		Flat genes are removed from the analysis after the core attractor pathway modules are first inferred (i.e. the \code{findAttractors} step).
}
\value{
		A \code{vector} with gene names (as defined in the eset) of those genes with expression profiles that hardly vary across 
		different celltype or experimental groups. 
}
\references{
		\code{limma} package.
		
		Smyth, G. K. (2004). Linear models and empirical Bayes methods for assessing differential expression in microarray experiments. Statistical Applications in Genetics and Molecular Biology 3, No. 1, Article 3.
}
\author{
		Jessica Mar
}
\examples{
data(subset.loring.eset)
remove.these.genes <- removeFlatGenes(subset.loring.eset, "celltype", contrasts=NULL, limma.cutoff=0.05)
}
\keyword{methods}
