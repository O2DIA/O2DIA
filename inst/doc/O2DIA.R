## ----setup, include=FALSE---------------------------------------------------------------------------------------------
options(width = 120)

## ----message=FALSE----------------------------------------------------------------------------------------------------
require(O2DIA)

## ---------------------------------------------------------------------------------------------------------------------
ls("package:O2DIA")

## ----comment=NA-------------------------------------------------------------------------------------------------------
args(O2DIAsteady)

## ----comment=NA-------------------------------------------------------------------------------------------------------
args(O2DIAdyna)

## ---------------------------------------------------------------------------------------------------------------------
P  <- O2DIAparms()
head(P)

## ---------------------------------------------------------------------------------------------------------------------
DYN <- O2DIAdyna()
head(O2DIAparms(DYN))

## ---------------------------------------------------------------------------------------------------------------------
O2DIAparms(DYN, which = "Q10_min")

## ---------------------------------------------------------------------------------------------------------------------
std <- O2DIAsteady()

print(O2DIAbudgetO2(std))

## ---------------------------------------------------------------------------------------------------------------------
range(O2DIAdx(std))
range(O2DIAgrid()$x.int)
head(O2DIAdepth(std))

## ---------------------------------------------------------------------------------------------------------------------
head(O2DIA1D(std), n = 3)

## ----fig.width=8, fig.height=8----------------------------------------------------------------------------------------
STD1 <- O2DIAsteady (Temperature =  0)   
STD2 <- O2DIAsteady (Temperature = 10)   
STD3 <- O2DIAsteady (Temperature = 20)
plot(STD1, STD2, STD3, 
     which = c("O2", "ODU"), 
     ylim  = c(0.01, 0), lty = 1, lwd = 2)
legend("bottom", legend = c(0, 10, 20), title = "dgC", 
       lty = 1, col = 1:3)

## ----fig.width=8, fig.height=6----------------------------------------------------------------------------------------
O2DIAparms(which = "prodfacdry")
out    <- O2DIAsteady()
outdry <- O2DIAsteady(Waterheight = 0)

plot(out, outdry, 
     which = c("O2", "ODU"), 
     lty = 1, lwd = 2)
legend("topright", title = "exchange", legend = c("water","dry"), 
       col = 1:2, lty = 1) 

## ---------------------------------------------------------------------------------------------------------------------
print(O2DIAbudgetO2(out, outdry))

## ----fig.width=8, fig.height=8----------------------------------------------------------------------------------------
times    <- seq (0, 1, length.out = 100)
PAR0_fun <- function(t) pmax(0, 1000*sin(2*pi*t - pi/2))

DIA <- O2DIAdyna (times = times,
                  PAR0 = PAR0_fun)

image2D(DIA, 
        mfrow = c(2, 3), las = 1)

plot(DIA, which = c("PAR0", "Waterheight"), mfrow=NULL)
matplot.0D(DIA, 
           which = c("O2SWIflux", "ODUSWIflux"), 
           main = "sediment-water exchange fluxes", 
           ylab = "mol/m2/d",
           mfrow = NULL, lty = 1, lwd = 2)
matplot1D(DIA, which = c("O2"), 
     mfrow = NULL, lwd = 2)

## ----fig.width=8, fig.height=8----------------------------------------------------------------------------------------
DIA2 <- O2DIAdyna (times = times,
                  PAR0 = PAR0_fun,
                  Waterheight = 0)

image2D(DIA2, 
        mfrow = c(2, 3), las = 1)

plot(DIA2, which = c("PAR0", "Waterheight"), mfrow=NULL)
matplot.0D(DIA2, 
           which = c("O2SWIflux", "ODUSWIflux"), 
           main = "sediment-water exchange fluxes", 
           ylab = "mol/m2/d",
           mfrow = NULL, lty = 1, lwd = 2)
matplot1D(DIA2, which = c("O2"), 
     mfrow = NULL, lwd = 2)

## ---------------------------------------------------------------------------------------------------------------------
O2DIAbudgetO2(DIA, DIA2)

## ----fig.width=8, fig.height=8----------------------------------------------------------------------------------------
times    <- seq (0, 1, length.out = 100)
WH_fun <- function(t) pmax(0, 2*sin(2*pi*t/0.517))  # tidal period = 12.4/24 days
PAR0_fun <- function(t) pmax(0, 1000*sin(2*pi*t - pi/2))*exp(-1*WH_fun(t))

DIA3 <- O2DIAdyna (times = times,
                  Waterheight = WH_fun,
                  PAR0 = PAR0_fun)

image2D(DIA3, 
        mfrow = c(2, 3), las = 1)

plot(DIA3, which = c("PAR0", "Waterheight"), mfrow=NULL)
matplot.0D(DIA3, 
           which = c("O2SWIflux", "ODUSWIflux"), 
           main = "sediment-water exchange fluxes", 
           ylab = "mol/m2/d",
           mfrow = NULL, lty = 1, lwd = 2)
matplot1D(DIA3, which = c("O2"), 
     mfrow = NULL, lwd = 2)

## ---------------------------------------------------------------------------------------------------------------------
knitr:::kable(O2DIAparms())

## ---------------------------------------------------------------------------------------------------------------------
knitr:::kable(O2DIAsvar())

## ---------------------------------------------------------------------------------------------------------------------
knitr:::kable(O2DIA0D())

## ---------------------------------------------------------------------------------------------------------------------
knitr:::kable(O2DIA1D())

