\name{O2DIAparms}
\alias{O2DIAparms}
\alias{O2DIApor}
\alias{O2DIAdx}
\alias{O2DIAdepth}
\alias{O2DIAgrid}
\title{
  Functions to retrieve parameters, porosity, depth and sediment grid for the O2DIA model.
}
\description{

  \code{O2DIAparms, O2DIAdepth, O2DIAdx} retrieve the parameters, sediment depths and layer thicknesses of O2DIA model solutions.
  
  \code{O2DIApor} retrieves the porosity of O2DIA model solutions.
}
\usage{
  O2DIAparms(out = NULL, as.vector = FALSE, which = NULL)
  O2DIAdepth(out = NULL)
  O2DIAgrid (out = NULL)
  O2DIAdx   (out = NULL)
  O2DIApor  (out = NULL)
}
\arguments{

  \item{out }{an output object returned by \link{O2DIAsteady} or \link{O2DIAdyna}. 
   If \code{NULL}, \code{O2DIAparms} will return the default (parameter) values.
  }
  \item{as.vector }{if \code{TRUE} will return the parameter \code{vector}, else a \code{data.frame} that also contains the units.
  }
  \item{which }{if not \code{NULL}, a vector with names of the variables/parameters to return.
  }
}

\author{
  Karline Soetaert
}

\examples{

  defparms  <- O2DIAparms(as.vector = TRUE)
  
  std1 <- O2DIAsteady()
  std2 <- O2DIAsteady(parms = list(por = 0.5))
  
  O2DIAparms(std1, as.vector = FALSE, which = NULL)
  O2DIAdepth(out = std1)
  O2DIAdepth(out = std2)

}

\references{
 to be added
}

\keyword{ utilities }
