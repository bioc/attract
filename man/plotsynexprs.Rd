\name{plotsynexprs}
\alias{plotsynexprs}
\title{
		Visualizing the Average Expression Profile of a Synexpression Group.
}
\description{
		This function plots the average expression profile for a specific synexpression group. 
}
\usage{
		plotsynexprs(mySynExpressionSet, tickMarks, tickLabels, vertLines, index=1, ...)
}

\arguments{
  \item{mySynExpressionSet}{
		\code{\link{SynExpressionSet}} object.
}
  \item{tickMarks}{
		numeric vector of specifying the location of the tick marks along the x-axis. There should be one tick for each cell type or group.
}
  \item{tickLabels}{
		character vector specifying the labels to be appear underneath the tick marks on the x-axis. These should correspond to the cell type or group names.
}
  \item{vertLines}{
		numeric vector specifying the location of the vertical lines that indicate the cell type or group-specific regions along the x-axis. 
}
  \item{index}{
		numeric value specifying which synexpression group should be plotted.
}
  \item{\dots}{
		additional arguments.
}
}
\details{
		Generic plotting parameters can be passed to this function to create a more sophisticated plot, e.g \code{col="blue"}, \code{main="Synexpression Group 1"}.
}
\value{
		A plot showing the average expression profile for the synexpression group specified. 
}

\author{
		Jessica Mar
}
\examples{
data(subset.loring.eset)
attractor.states <- findAttractors(subset.loring.eset, "celltype", nperm=10, annotation="illuminaHumanv1.db")
remove.these.genes <- removeFlatGenes(subset.loring.eset, "celltype", contrasts=NULL, limma.cutoff=0.05)
mapk.syn <- findSynexprs("04010", attractor.states, remove.these.genes)
par(mfrow=c(2,2)) 
pretty.col <- rainbow(3) 
for( i in 1:3 ){
	plotsynexprs(mapk.syn, tickMarks=c(6, 28, 47, 60), tickLabels=c("ESC", "PRO", "NSC", "TER"), vertLines=c(12.5, 43.5, 51.5), index=i, 
			main=paste("Synexpression Group ", i, sep=""), col=pretty.col[i])
 }
}
\keyword{aplot}
