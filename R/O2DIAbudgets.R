#####################################################################################################
######                           O2DIA: C, N, P, O2 diagenesis                                ######
######                                     BUDGETTING                                          ######
#####################################################################################################

# Utility functions

getMean0D <- function(out){
  OUT <- O2DIA0D(out)
  
  on  <- OUT$name
  OUT <- as.list(OUT$value)
  names(OUT) <- on
  return(OUT)
}

#====================#
# budget wrapper     #
#====================#

## -----------------------------------------------------------------------------

O2DIAbudget_all <- function(out, ..., 
   which = c("All", "Rates", "Fluxes"), 
   func = O2DIAbudgetO2_one, args) 
  
  {
     which <- match.arg(which, 
                     choices = c("All", "Rates", "Fluxes"))
     if ("All" %in% which)
       which <- c("Rates", "Fluxes", "dC")
     
     ALL  <- list(out, ...)
     
     NM <- unlist(lapply(args[-1], as.character))
     if (! is.null(names(NM)))
       NM <- NM[!names(NM) == "which"]
     if (length(NM) != length(ALL)) NM <- paste("out", 1:length(ALL), sep = "_")
     
     budg <- func(out)  
     budgFlux    <- budg$Fluxes
     budgRates   <- budg$Rates
     budgdC      <- budg$dC

     if (length(ALL) > 1) {

       budgFlux    <- unlist(budgFlux)
       budgRates   <- unlist(budgRates)
       for ( i in 2:length(ALL)) {
        budg <- func(ALL[[i]])
    
        budgFlux    <- cbind(budgFlux,    unlist(budg$Fluxes))
        budgRates   <- cbind(budgRates,   unlist(budg$Rates))
        budgdC      <- cbind(budgdC,          budg$dC)
    } 
  
    cn <- rep(names(budg$Fluxes), each = 3)
    nc <- nchar(cn)
    cn <- paste (cn, c("surf", "deep",  "net"), sep="")

    rownames(budgFlux) <- cn
     }
     
     budg <- list()
     if ("Rates" %in% which){
      if(length(NM) > 1) colnames(budgRates) <- NM
      budg$Rates <- budgRates
    }
    if ("Fluxes" %in% which){
      if(length(NM) > 1) colnames(budgFlux) <- NM
      budg$Fluxes <- budgFlux
    }
    if ("dC" %in% which){
      if(length(NM) > 1) colnames(budgdC) <- NM
      budg$dC <- budgdC
    }
    return(budg)
}    
  
O2DIAbudgetO2 <- function(out, ..., 
                           which = c("All", "Rates", "Fluxes")) 
  O2DIAbudget_all(out, ..., which = which, 
                   func = O2DIAbudgetO2_one, args = sys.call())  

#====================#
#  O2 budgets        #
#====================#

O2DIAbudgetO2_one <- function(out) {

  dC <- NULL
  if (inherits(out, "deSolve")){
    out <- getMean0D(out)
  }
  
  Cflux <- out$O2SWIflux - out$O2deepflux

  Fluxes <- data.frame(O2  = c(out$O2SWIflux, out$O2deepflux),
                       ODU = c(out$ODUSWIflux, out$ODUdeepflux))
  
  Fluxes           <- rbind(Fluxes, Fluxes[1,] - Fluxes[2,])                     
  rownames(Fluxes) <- c("surface", "bottom", "netinflux")

  Rates <- data.frame(
    O2Production       = out$Totalprod,
    Growthresp_min     = out$Totalmin,
    Basalresp          = out$Totalbasalresp,
    Photoresp          = out$Totalphotoresp,
    ODUreoxidation     = out$Totalodureox,
    TotalConsumption   = out$Totalcons
  )

  Rates$NetConsumption = Rates$TotalConsumption - Rates$O2Production

  rownames(Rates) <- "molO2/m2/d" 
  
  dC <- c(O2  = out$O2SWIflux  - out$O2deepflux  - Rates$NetConsumption,
          ODU = out$ODUSWIflux - out$ODUdeepflux - Rates$ODUreoxidation)
  
  return(list(
     Fluxes = Fluxes, Rates = Rates, dC = dC))                 
 
}

## -----------------------------------------------------------------------------
