\name{O2DIA-package}
\alias{O2DIA-package}
\alias{O2DIA}
\docType{package}
\title{
  An early diagenetic model for the O2 cycle in marine sediments.
}
\description{

  O2DIA is a 1-D model of oxygen diagenesis in a (marine) sediment.

  The model describes two state variables, in \bold{300} layers:
  \itemize{
    \item Oxygen (O2), dissolved substance.
    \item Oxygen demand unit (ODU), dissolved substance, as a
      lump sum of all reduced substances that can be reoxidised consuming oxygen, including ammonium.
  }
  Concentrations of are expressed in [mol O2/m3 liquid].

  The model is implemented in fortran and linked to R.

}
\author{
  Karline Soetaert, Qi Liu
}
\references{
  Soetaert K, PMJ Herman and JJ Middelburg, 1996a.
  A model of early diagenetic processes from the shelf to abyssal depths.
  Geochimica Cosmochimica Acta, 60(6):1019-1040.

  Soetaert K, PMJ Herman and JJ Middelburg, 1996b.
  Dynamic response of deep-sea sediments to seasonal variation: a model.
  Limnol. Oceanogr. 41(8): 1651-1668.
}

\keyword{ package }
