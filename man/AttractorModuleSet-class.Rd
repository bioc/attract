\name{AttractorModuleSet-class}
\docType{class}
\alias{AttractorModuleSet}
\alias{AttractorModuleSet-class}

\title{Class AttractorModuleSet }
\description{ This is a class representation for storing the output of the \code{findAttractors} function.	}
\section{Objects from the Class}{
	Objects are output by the function \code{\link{findAttractors}}. 
	Objects can also be created by using \code{new("AttractorModuleSet", ...)}.
}
\section{Slots}{
	 \describe{
    \item{\code{eSet}:}{ ExpressionSet which primarily stores the expression data and the phenotype/sample data sets. }
    \item{\code{cellTypeTag}:}{ character string of the tag which stores the group membership information for the samples. Must be a column name of the data frame pData(eset). }
	\item{\code{incidenceMatrix}:}{ incidence matrix used as input to GSEAlm. }
	\item{\code{rankedPathways}:}{ Data frame of significantly enriched pathways, ranked first by significance and then by size. }
  }
}
\section{Methods}{
No methods have yet been defined with class "AttractorModuleSet" in the signature.
}
\author{ Jessica Mar \email{jess@jimmy.harvard.edu} }
\examples{
\dontrun{
new.attractmodule <- new("AttractorModuleSet", eSet=new("ExpressionSet"), cellTypeTag=character(1), incidenceMatrix=matrix(0), rankedPathways=data.frame())
}
}
\note{ This class is better describe in the vignette. }
\keyword{classes}
