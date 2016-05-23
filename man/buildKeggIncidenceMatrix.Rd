\name{buildKeggIncidenceMatrix}
\alias{buildKeggIncidenceMatrix}
\title{
		Internal function - buildings the KEGG incidence matrix requird by GSEAlm.
}
\description{
		Internal function - buildings the KEGG incidence matrix requird by GSEAlm.
}
\usage{
buildKeggIncidenceMatrix(kegg.ids, gene.ids, annotation, database, analysis, envPos, expressionSetGeneFormat)
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
  \item{database}{
    a character string specifiying what pathway database you would like to use.  
}
  \item{analysis}{
    a character string specifying what type of experiment you performed, microarray or RNAseq.
}
  \item{envPos}{
  	the position of the annotation package in the R search path. 
}
  \item{expressionSetGeneFormat}{
  a character string specifying the type of identifier for a gene in your expression data set.
  The default value is NULL. (ex. SYMBOL, ENTREZID, REFSEQ, ENSEMBL)
}
}
\details{
	This function is called internally by \code{findAttractors}.
}
\value{
	A matrix object with 0 and 1 entries where 1 denotes membership of a gene in a KEGG or reactome pathway, 0 denotes non-membership. 
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
require("illuminaHumanv2.db", character.only=TRUE)
loadNamespace("illuminaHumanv2.db")
envPos <- match(paste("package:", "illuminaHumanv2.db", sep=""), search())
kegg.ids <- ls(illuminaHumanv2PATH2PROBE) 
gene.ids <- ls(illuminaHumanv2PATH)
database <-"KEGG"
analysis <- "microarray"
imat <- buildKeggIncidence(kegg.ids, gene.ids, illuminaHumanv2.db,database, analysis, envPos) 
}
}
\keyword{internal}
