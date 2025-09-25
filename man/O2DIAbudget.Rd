\name{O2DIA budget}
\alias{O2DIAbudgetO2}
\title{
  Budget for the O2DIA model.
}
\description{
 \code{O2DIAbudgetO2} estimates the oxygen mass 
 budget from O2DIA model solutions.
}
\usage{
  O2DIAbudgetO2(out, ..., which = c("All", "Rates", "Fluxes")) 
}
\arguments{

  \item{out }{an output object returned by \link{O2DIAsteady} or \link{O2DIAdyna}. 
  }

  \item{which }{if not \code{NULL}, a vector with names of the items to return.
  }
  \item{... }{unused.
  }
}
\value{
 \code{O2DIAbudgetO2} returns the oxygen budget of the solution as a \code{list}, with the following items. 
 
 \itemize{
  \item \code{Fluxes}, the boundary fluxes at the surface and bottom of the sediment, 
  and the net input, in mol/m2/d. Positive fluxes are directed into the sediment for the surface, and out of the sediment at the bottom. For dynamic runs, fluxes are averaged over the simulation period.
  \item \code{Rates}, the integrated process rates, in mol/m2/d
  \item \code{dC}, the rate of change of integrated values of state vairables, in mmol/m2/d, defined as ([concentration at the end] - [concentration at beginning])/[length of simulation]
 }

If more than one O2DIA object is passed to these functions, a matrix is returned, one column for each object. 

For dynamic runs, the budget is taken over the mean of the run.

}

\author{
  Karline Soetaert
}


\examples{

# steady-state runs
  std1  <- O2DIAsteady()
  std2  <- O2DIAsteady(Temperature = 20)
  O2DIAbudgetO2(std1, std2)         # temperature effect
 
# dynamic
  times <- seq(0, 1, length.out = 100)
  out1  <- O2DIAdyna(times = times)
  
  out2  <- O2DIAdyna(times = times,
                     Temperature = 20)
  out3  <- O2DIAdyna(times = times,
                     Temperature = function(t) 20 - 2*sin(2*pi*t))

  plot(out2, out3, 
       which=c("Temperature", "O2SWIflux"), 
       lty = 1)
       
# budgets
  O2DIAbudgetO2(std2, out2)         # the same
  
  O2DIAbudgetO2(out2, out3)$Rates   # very similar
  
  O2DIAbudgetO2(out1, out2, out3, which = "Fluxes")
}

\references{

 to be added
}

\keyword{ utilities }
