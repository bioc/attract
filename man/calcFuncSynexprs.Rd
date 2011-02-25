\name{calcFuncSynexprs}
\alias{calcFuncSynexprs}
\title{
		Functional enrichmental analysis for a set of synexpression groups.
}
\description{
		This function performs functional enrichment for a given set of synexpression groups.
}
\usage{
calcFuncSynexprs(mySynExpressionSet, myAttractorModuleSet, ontology = "BP", min.pvalue = 0.05, min.pwaysize = 5, annotation = "illuminaHumanv2BeadID.db", ...)
}
\arguments{
  \item{mySynExpressionSet}{
		\code{SynExpressionSet} object.
}
  \item{myAttractorModuleSet}{
		\code{AttractorModuleSet} object.
}
  \item{ontology}{
		character string specifying which GO ontology to use, either "MF", "BP", or "CC"; defaults to "BP".
}
  \item{min.pvalue}{
		numeric value specifying adjusted P-value cut-off to use, categories with P-values <= min.pvalue will be reported.
}
  \item{min.pwaysize}{
		\code{integer} specifying minimum size of the pathway or category to consider for enrichment analysis.
}
  \item{annotation}{
		character string specifying the annotation package that corresponds to the chip platform the data was generated from. 
}
  \item{\dots}{
		additional arguments.
}
}
\details{
		This function performs a functional enrichment analysis on each synexpression group using the \code{hyperGTest} from 
		the \code{GOstats} package. P-values are adjusted using the Benjamini-Hochberg correction method. Results are returned
		only if they satisfy the minimum P-value level, as specified by the \code{min.pvalue} argument. 
}
\value{
	A \code{list} object. 
}
\author{
		Jessica Mar
}
\references{
Falcon, S. and R. Gentleman, Using GOstats to test gene lists for GO term association. Bioinformatics, 2007. 23(2): p. 257-8.
}
\examples{
\dontrun{
data(subset.loring.eset)
attractor.states <- findAttractors(subset.loring.eset, "celltype", nperm=10, annotation="illuminaHumanv1.db")
remove.these.genes <- removeFlatGenes(subset.loring.eset, "celltype", contrasts=NULL, limma.cutoff=0.05)
mapk.syn <- findSynexprs("04010", attractor.states, remove.these.genes)
mapk.func <- calcFuncSynexprs(mapk.syn, attractor.states, "CC", annotation="illuminaHumanv1.db") 
}
}
\keyword{methods}
