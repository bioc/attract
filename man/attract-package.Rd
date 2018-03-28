\name{attract-package}
\alias{attract-package}
\alias{attract}
\docType{package}
\encoding{latin1}
\title{
Methods to find the Gene Expression Modules that Represent the Drivers of Kauffman's Attractor Landscape
}
\description{
This package contains functions used to determine the gene expression modules that represent the drivers of Kauffman's attractor landscape.
}
\details{
\tabular{ll}{
Package: \tab attract\cr
Type: \tab Package\cr
Version: \tab 1.31.1\cr
Date: \tab 2018-03-28\cr
License: \tab \cr
LazyLoad: \tab yes\cr
}
The method can be summarized in the following key steps: 
(1) Determine core KEGG or reactome pathways that discriminate the most strongly between celltypes or experimental groups of interest (see \code{findAttractors)}). 
(2) Find the different synexpression groups that are present within a core attractor pathway (see \code{findSynexprs}).
(3) Find sets of genes that show highly similar profiles to the synexpression groups within an attractor pathway module (see \code{findCorrPartners}). 
(4) Test for functional enrichment for each of the synexpression groups to detect any potentially shared biological themes (see \code{calcFuncSynexprs}). 
}
\author{
		Jessica Mar   <jess@jimmy.harvard.edu>
}
\references{
	Kauffman S. 2004. A proposal for using the ensemble approach to understand genetic regulatory networks. J Theor Biol. 230:581. 
	Mar JC, Wells CA, Quackenbush J. 2010. Identifying Gene Expression Modules that Represent the Drivers of Kauffman's 
Attractor Landscape. To Appear. 
M\"{u}ller F et al. 2008. Regulatory networks define phenotypic classes of human stem cell lines. Nature. 455(7211): 401. 
Mar JC, Wells CA, Quackenbush J. 2010. Defining an Informativeness Metric for Clustering Gene Expression Data. To Appear. 
}
\examples{
\dontrun{
data(subset.loring.eset)
attractor.states <- findAttractors(subset.loring.eset, "celltype", nperm=10, annotation="illuminaHumanv1.db")
remove.these.genes <- removeFlatGenes(subset.loring.eset, "celltype", contrasts=NULL, limma.cutoff=0.05)
mapk.syn <- findSynexprs("04010", attractor.states, remove.these.genes)
mapk.cor <- findCorrPartners(mapk.syn, subset.loring.eset, remove.these.genes)
mapk.func <- calcFuncSynexprs(mapk.syn, attractor.states, "CC", annotation="illuminaHumanv1.db") 
}
}
\keyword{ package }
