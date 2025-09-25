
## ====================================================================
## Helper functions
## ====================================================================

CreateO2DIAgrid <- function(parms){
  N        <- .O2DIA$N
  DBL      <- parms["x_DBL"]
  x.down   <- parms["x_down"]
  x.up     <- -DBL 
  x.SWI    <- 0
  
  # grid with equally sized boxes (not used here, see next R-command)
  dblx1    <- min(DBL/20, 1e-6)
  Grid     <- setup.grid.1D(x.up = x.up,
                            L    = c(DBL, x.down), 
                            dx.1 = c(dblx1, 1e-6), 
                            N    = c(20,    N-20)) #grid is in m
  Grid
  
}

## ====================================================================
## Steady-state function
## ====================================================================

O2DIAsteady <- function(parms = list(), yini = NULL, Grid = NULL, times = 0, 
                        Temperature = 10, PAR0 = 1000, Waterheight = 0.01,
                        Tfunction = "Q10"){
  
  # parameter vector
  P   <- .O2DIA$Parms
  np  <- names(parms)
  if (any (!np %in% names(P))) 
    stop("parameter ", np[(!np %in% names(P))], "not known")
  P[np] <- unlist(parms)
  
  if (Tfunction != "Q10")
    P["T0"] <- 99999.   # above 50
  
  N        <- .O2DIA$N

  # spatial domain
  if (is.null(Grid)) Grid <- CreateO2DIAgrid(P)
  
  # parameters in dll - IN THAT ORDER
  dllparms <- as.double (c(P[-(1:2)], 
                           Grid$x.mid, Grid$x.int,
                           Grid$dx,    Grid$dx.aux ))
  
  # forcings: temp, PAR, waterlevel   - KEEP ORDER
  forcs <- vector()
  
  if (is.vector(Temperature))
    forcs[1] <- Temperature[1]
  else if (is.function(Temperature))
    forcs[1] <- Temperature(times[1])
  else if (is.matrix(Temperature) | is.data.frame(Temperature))
    forcs[1] <- approx(Temperature[,1:2], xout = times[1])$y
  
  if (is.vector(PAR0))
    forcs[2] <- PAR0[1]
  else if (is.function(PAR0))
    forcs[2] <- PAR0(times[1])
  else if (is.matrix(PAR0) | is.data.frame(PAR0))
    forcs[2] <- approx(PAR0[,1:2], xout = times[1])$y
  
  if (is.vector(Waterheight))
    forcs[3] <- Waterheight[1]
  else if (is.function(Waterheight))
    forcs[3] <- Waterheight(times[1])
  else if (is.matrix(Waterheight) | is.data.frame(Waterheight))
    forcs[3] <- approx(Waterheight[,1:2], xout = times[1])$y
  
  if (any(is.na(forcs))) 
    stop ("NAs in Temperature, PAR0 (light) or Waterheight forcing function")
  
  # initial condition
  if (is.null(yini)) yini <- runif(2*N)
  
  # names of variables
  names    <- .O2DIA$ynames     # state variable names
  outnames <- .O2DIA$outnames   # output variable names
  
  # solve 
  std <- steady.1D(
    y = yini, names = names, 
    parms = dllparms, forcings = forcs,
    func = "o2mod", initfunc = "inito2mod",  initforc="inito2forc", 
    dllname = "O2DIA",
    times = times[1], atol = 1e-8, rtol=1e-8, pos=TRUE,
    nspec = length(names), dimens = N, 
    outnames=outnames, nout=length(outnames))
  
  std$Grid  <- Grid  
  std$Depth <- Grid$x.mid   
  std$porosity <-  rep(P["por"], times = N)
  
  std$Parms <- P
  class(std) <- c("O2DIAstd", class(std))
  
  std  
}

## ====================================================================
## Dynamic function
## ====================================================================

O2DIAdyna <- function(parms = list(),  times = seq(0, 1, length.out=10),
                      yini = NULL, Grid = NULL, 
                      Temperature = 10, PAR0 = 1000, Waterheight = 0.01,
                      Tfunction = "Q10", verbose = FALSE){
  
  # parameter vector
  P   <- .O2DIA$Parms
  np <- names(parms)
  if (any (!np %in% names(P))) 
    stop("parameter ", np[(!np %in% names(P))], "not known")
  P[np] <- unlist(parms)
  
  if (Tfunction != "Q10")
    P["T0"] <- 99999.   # above 50
  
  # Spatial domain
  N    <- .O2DIA$N
  
  if (is.null(Grid)) Grid <- CreateO2DIAgrid(P)
  
  # parameters in dll - IN THAT ORDER
  dllparms <- as.double (c(P[-(1:2)], # remove first 2 parameters (DBL, x)
                           Grid$x.mid, Grid$x.int,
                           Grid$dx,    Grid$dx.aux))
  
  # forcings: temp, PAR, waterlevel   - KEEP ORDER
  
  forcings <- list()
  if (is.function(Temperature))
    forcings[[1]] <- cbind(times, Temperature(times))
  else if (is.vector(Temperature) & length(Temperature) == 1)
    forcings[[1]] <- data.frame(times, Temperature)
  else if (is.matrix(Temperature) | is.data.frame(Temperature))
    forcings[[1]] <- as.data.frame(approx(Temperature[,1:2], xout = times))
  if (any(is.na(forcings[[1]]))) 
    stop ("NAs in temperature forcing function")
  
  if (is.function(PAR0))
    forcings[[2]] <- cbind(times, PAR0(times))
  else if (is.vector(PAR0) & length(PAR0) == 1)
    forcings[[2]] <- data.frame(times, PAR0)
  else if (is.matrix(PAR0) | is.data.frame(PAR0))
    forcings[[2]] <- as.data.frame(approx(PAR0[,1:2], xout = times))
  if (any(is.na(forcings[[2]]))) stop ("NAs in PAR0 (light) forcing function")
  
  if (is.function(Waterheight))
    forcings[[3]] <- cbind(times, Waterheight(times))
  else if (is.vector(Waterheight) & length(Waterheight) == 1)
    forcings[[3]] <- data.frame(times, Waterheight)
  else if (is.matrix(Waterheight) | is.data.frame(Waterheight))
    forcings[[3]] <- as.data.frame(approx(Waterheight[,1:2], xout = times))
  if (any(is.na(forcings[[3]]))) stop ("NAs in Waterheight forcing function")
  
  if (is.null(yini))
    yini <- O2DIAsteady(times = times[1], parms = parms, 
                        Temperature = Temperature, PAR0 = PAR0, 
                        Waterheight = Waterheight, 
                        Grid = Grid)
  
  # names of variables
  names    <- .O2DIA$ynames     # state variable names
  outnames <- .O2DIA$outnames   # output variable names
  
  # Solve it
  ZZ <- NULL
  if (! verbose)
    ZZ <- capture.output(suppressWarnings(DYN <- ode.1D(
      y = yini$y, names = names, 
      parms = dllparms, forcings = forcings,
      func = "o2mod", initfunc = "inito2mod",  initforc="inito2forc", 
      dllname = "O2DIA",
      times = times, atol = 1e-8, rtol=1e-8, 
      nspec = length(names), dimens = N, 
      outnames = outnames, nout = length(outnames))))
  else
    DYN <- ode.1D(
      y = yini$y, names = names, 
      parms = dllparms, forcings = forcings,
      func = "o2mod", initfunc = "inito2mod",  initforc="inito2forc", 
      dllname = "O2DIA",
      times = times, atol = 1e-8, rtol=1e-8, 
      nspec = length(names), dimens = N, 
      outnames = outnames, nout = length(outnames))
  
  colnames(DYN)[2:(2*N+1)] <- as.vector(
    sapply(names, FUN = function(x) rep(x, times = N)))
  
  # add attributes
  attributes(DYN)$Grid     <- Grid
  attributes(DYN)$Depth    <- Grid$x.mid
  attributes(DYN)$porosity <- rep(P["por"], times = N)
  attributes(DYN)$Parms <- P
  class(DYN) <- c("O2DIAdyn", class(DYN))
  DYN
}
