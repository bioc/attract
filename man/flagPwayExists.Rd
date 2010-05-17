\name{flagPwayExists}
\alias{flagPwayExists}
\title{
		Internal function - flags a gene if it exists in a pathway. 
}
\description{
		Internal function - flags a gene if it exists in a pathway. 
}
\usage{
flagPwayExists(x)
}

\arguments{
  \item{x}{
		A \code{vector} of pathway ids.
}
}
\details{
		This function is called internally by the \code{findAttractors} function. Function returns TRUE if the probe exists in at least one pways, FALSE if there are no pway annotations.
}
\value{
		A \code{logical} value.
}
\author{
		Jessica Mar
}
\examples{
\dontrun{
library(illuminaHumanv2.db)
flag.check <- flagPwayExists(ls(illuminaHumanv2PATH2PROBE))
}
}
\keyword{internal}

