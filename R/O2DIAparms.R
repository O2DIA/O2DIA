##===========================================
## Interrogation Functions for O2DIA models
##===========================================

##------------------------------------
## Get parameters and values
##------------------------------------

O2DIAparms <- function(out = NULL, as.vector = FALSE, which = NULL) {
  if (is.null(out))
    Parms <- .O2DIA$Parms
  else if (inherits(out, "steady1D"))
    Parms <- out$Parms[1:length(.O2DIA$Parms)]
  else if (inherits(out, "deSolve"))
    Parms <- attr(out, "Parms")[1:length(.O2DIA$Parms)]
  else stop("object 'out' not supported")
  
  if (as.vector) {
    if (! is.null(which))
      Parms <- Parms[which]
    return(Parms)
  } else {
    Units <- .O2DIA$Parunit
    if  (Parms["BCupLiq"] == 1) {
      i1 <- which(names(Parms) == "O2bw")
      Units[i1:(i1+6)] <- "mmol/m2/d"
    }  
    if  (Parms["BCdownLiq"] == 1) {
      i <- which(names(Parms) == "O2dw")
      Units[i1:(i1+6)] <- "mmol/m2/d"
    }  
    Parms <- data.frame(parms = Parms, units = Units, description = .O2DIA$Pardesc)
    if (! is.null(which))
      Parms <- Parms[which, ]
    return(Parms)
  }
}

##------------------------------------
## Grid
##------------------------------------

O2DIAgrid <- function(out = NULL) {
  if (is.null(out))
    D <- .O2DIA$Grid 
  else if (inherits(out, "steady1D"))
    D <- out$Grid
  else if (inherits(out, "deSolve"))
    D <- attributes(out)$Grid
  else stop("object 'out' not supported for grid calculation - try O2DIAdepth intstead")
  D
} 

O2DIAdepth <- function(out = NULL) {
 if (is.null(out))
   D <- .O2DIA$Grid$x.mid
 else if (inherits(out, "steady1D"))
   D <- out$Depth
 else if (inherits(out, "deSolve"))
   D <- attr(out, "Depth")
 else stop("object 'out' not supported")
 D
} 

O2DIAdx <- function(out = NULL) {
  if (is.null(out))
    D <- .O2DIA$Grid$dx
  else if (inherits(out, "steady1D"))
    D <- out$Grid$dx
  else if (inherits(out, "deSolve"))
    D <- attr(out, "Grid")$dx
  else stop("object 'out' not supported")
  D
} 

O2DIApor <- function(out = NULL) {
  if (is.null(out))
    p0 <- .O2DIA$Parms$por
  if (inherits(out, "steady1D"))
    p0 <- out$porosity
  else if (inherits(out, "deSolve"))
    p0 <- attr(out, "porosity")
  else stop("object 'out' not supported")
  x <- O2DIAdepth((out))
  pmin(1.0, p0 + (1.0-p0)*exp(-x/0.01e-3)) 
} 


##------------------------------------
## Get variables 
##------------------------------------

MeanVal <- function(out)  # takes into account unequal timing 
  (colSums(diff(out[,1])*(out[-1,]+out[-nrow(out),])*0.5)/(out[nrow(out),1]-out[1,1]))[-1]

O2DIA0D <- function(out, as.vector = FALSE, which = NULL) {
  if (missing(out)) {
    Dnames <- c(.O2DIA$var0D,.O2DIA$varforc)
    D <- rep(NA, times = length(Dnames))
    names(D) <- Dnames
  }  else if (inherits(out, "steady1D"))
   D <- unlist(out[c(.O2DIA$var0D,.O2DIA$varforc)])
  else if (inherits(out, "deSolve"))
    D <- MeanVal(out[, c("time",.O2DIA$var0D,.O2DIA$varforc)])
  else stop("object 'out' not supported")
 
 if (! as.vector)
   D <- data.frame(names = names(D), values = D, 
      units = c(.O2DIA$unit0D,.O2DIA$unitforc),
      description = c(.O2DIA$descrip0D,.O2DIA$descripforc))
 
 if (! is.null(which)){
   if (is.vector(D))
     D <- D[which]
   else D <- D[which,]  
 } 
 row.names(D) <- NULL
 D
} 

O2DIA1D <- function(out, which = NULL) {
  
  if (missing(out)) {
     D <- data.frame(names = c(.O2DIA$svar, .O2DIA$var1D), 
            units = c(.O2DIA$yunits, .O2DIA$unit1D), 
            description = c(.O2DIA$ydescrip, .O2DIA$descrip1D))
     if (! is.null(which))
       D <- D[(D$names %in% which), ]
     return(D)
  }
  
  if (inherits(out, "steady1D"))
    D <- cbind(out$y, as.data.frame(out[.O2DIA$var1D]))
  else if (inherits(out, "deSolve")) {
    D <- NULL
    for (cc in c(.O2DIA$svar, .O2DIA$var1D))
      D <- cbind(D,MeanVal(cbind(out[,1],subset(out, which = cc))))
    rownames(D) <- NULL
    colnames(D) <- c(.O2DIA$svar,.O2DIA$var1D)
    D <- as.data.frame(D)  
  }
  else stop("object 'out' not supported")
 
  D <- cbind(x = O2DIAdepth(out), por = O2DIApor(out), D)
  if (! is.null(which))
    D <- D[ ,c( "x", "por", which)]
  D  
} 
  
O2DIAsvar <- function(out, which = NULL) {
  if (missing(out)) 
    return(data.frame(names = .O2DIA$svar, 
                      units = .O2DIA$yunits, 
                      description = .O2DIA$ydescrip))
  
  if (inherits(out, "steady1D"))
    D <- out$y
  else if (inherits(out, "deSolve")) {
    D <- NULL
    for (cc in .O2DIA$svar)
      D <- cbind(D,MeanVal(cbind(out[,1],subset(out, which = cc))))
    rownames(D) <- NULL
    colnames(D) <- .O2DIA$svar
    D <- as.data.frame(D)  
  }
  else stop("object 'out' not supported")
  
  D <- cbind(x = O2DIAdepth(out), por = O2DIApor(out), D)
  if (! is.null(which))
    D <- D[ ,c( "x", "por", which)]
  D  
} 


## ============================================================================
## ============================================================================
##   Functions to extract parameters and variables from MPBDIA models
## ============================================================================
## ============================================================================

O2DIAparms <- function(out = NULL, as.vector = FALSE, which = NULL) {
  if (is.null(out))
    Parms <- c(.O2DIA$Parms)
  else if (inherits(out, "steady1D"))
    Parms <- out$Parms
  else if (inherits(out, "deSolve"))
    Parms <- attr(out, "Parms")
  else stop("object 'out' not supported")
  
  Parms <- Parms[!(names(Parms) %in% c("" , "model"))]
  if (as.vector) {
    if (! is.null(which))
      Parms <- Parms[which]
    return(Parms)
  } else {
    Pn <- names(Parms)
    Units <- c(.O2DIA$Parunit)
    Parms <- data.frame(parms = Parms, units = Units, 
       description = c(.O2DIA$Pardesc))
    row.names(Parms) <- Pn   
    if (! is.null(which))
      Parms <- Parms[which, ]
    return(Parms)
  }
}

