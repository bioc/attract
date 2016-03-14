\name{SynExpressionSet-class}
\docType{class}
\alias{SynExpressionSet}
\alias{SynExpressionSet-class}

\title{Class SynExpressionSet }
\description{ This is a class representation for storing synexpression group information.	}
\section{Objects from the Class}{
	Objects are output by the function \code{\link{findSynexprs}}. 
	Objects can also be created by using \code{new("SynExpressionSet", ...)}.
}
\section{Slots}{
	 \describe{
    \item{\code{groups}:}{A \code{list} object denoting the probes or gene IDs (rnaseq) belonging to each synexpression group. }
    \item{\code{profiles}:}{ A \code{matrix} of average expresson profiles for each synexpression group, each profile is stored as a row. }
  }
}
\section{Methods}{
No methods have yet been defined with class "SynExpressionSet" in the signature.
}
\author{ Jessica Mar \email{jess@jimmy.harvard.edu} }
\examples{
new.synexpressionset <- new("SynExpressionSet", groups=list(), profiles=matrix(0))
}
\note{ This class is described in more detail in the vignette. }
\keyword{classes}
