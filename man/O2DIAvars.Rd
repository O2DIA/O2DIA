\name{O2DIAvars}
\alias{O2DIA0D}
\alias{O2DIA1D}
\alias{O2DIAsvar}
\title{
  Functions to retrieve variables for the O2DIA model.
}
\description{

  \code{O2DIA0D, O2DIA1D, O2DIAsvar} retrieve the (0-dimensional or 1-dimensional) output variables or the state variables of O2DIA model solutions.  When called with \code{out} not specified, will return the names of these variables and their units.
}
\usage{
  O2DIA0D  (out, as.vector = FALSE, which = NULL) 
  O2DIA1D  (out, which = NULL) 
  O2DIAsvar(out, which = NULL) 
}
\arguments{

  \item{out }{an output object returned by \link{O2DIAsteady}, or \link{O2DIAdyna}. 
   If not specified, the names, units and description will be returned, not the (mean) values.
  }
  \item{as.vector }{if \code{TRUE} will return the parameter vector, else a data.frame that also contains the units.
  }
  \item{which }{if not \code{NULL}, a vector with names of the variables/parameters to return.
  }
}
\value{
 \code{O2DIA0D} and \code{O2DIA1D} return the output variables of the solution as a vector or data.frame. 
 For dynamic runs, the output is averaged over the mean of the run.
 
 \code{O2DIA1D} always returns the sediment depth and the porosity as the first two columns.
}

\author{
  Karline Soetaert
}


\examples{

# defaults
  O2DIAsvar()
  head(O2DIA0D())
  O2DIA1D()
  
# some runs to work with  
  defsteady <- O2DIAsteady()
  defdyn    <- O2DIAdyna()

# altered steady-state run
  out <- O2DIAsteady(Temperature = 20)

# 0-D outputs
  cbind(steady = O2DIA0D(defsteady, as.vector = TRUE), 
        dyna   = O2DIA0D(defdyn   , as.vector = TRUE), 
        ref    = O2DIA0D(out))

# 1-D outputs
  head(O2DIA1D(out))
  head(O2DIA1D(defdyn, which = c("O2", "ODU")))

}

\references{
  Soetaert K, PMJ Herman and JJ Middelburg, 1996a.
  A model of early diagenetic processes from the shelf to abyssal depths.
  Geochimica Cosmochimica Acta, 60(6):1019-1040.

  Soetaert K, PMJ Herman and JJ Middelburg, 1996b.
  Dynamic response of deep-sea sediments to seasonal variation: a model.
  Limnol. Oceanogr. 41(8): 1651-1668.
}

\keyword{ utilities }
