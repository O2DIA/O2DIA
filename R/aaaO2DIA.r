## ====================================================================
## A local environment for non user-visible data,
## ====================================================================
.O2DIA <- new.env()

.O2DIA$N     <- 300

# length units are m
.O2DIA$Grid  <- setup.grid.1D(x.up=0, dx.1=0.0001, N = .O2DIA$N, L = 1)

##------------------------------------
## Parameters
##------------------------------------

.O2DIA$Parms <- c(
    x_DBL    = 0.2e-3,       # [m] thickness of diffusive boundary layer 
    x_down   = 1e-2,         # [m] thickness of sediment   
    por      = 0.9,          # [-] volumetric porosity 
    salinity = 31,           # [-] salinity  
    v_adv    = 0,            # [m/d] sediment accretion (=advection) velocity
    ks_min   = 0.00001,      # [mol/m3] O2 affinity in mineralisation
    ks_odureox   = 0.00001,  # [mol/m3] O2 limitation in ODU oxidation
    ks_basalresp = 0.00001,  # [mol/m3] O2 limitation in dark respiration
    ks_WL    = 0.00001,      # [m] limitation in GPP-water level relation
    biomass  = 1,            # [molC/m2] MPB biomass
    c_min0       = 0.01,     # [-] proportion of MPB production exudate EPS at T0
    r_odureox0   = 100,      # [/d] 1-order O2 consumption rate of odu reoxidation at T0
    r_mpbprod0   = 50,       # [/d] 1-order O2 production rate at T0 (photosynthesis)
    r_basalresp0 = 0.01,     # [/d] 1-order O2 dark respiration rate at T0
    c_photoresp0 = 0.001,    # [-] fraction O2 photorespiration 
    kext_w   = 0.0001e3,     # [/m] light extinction coefficient in water
    kext_s   = 20e3,         # [/m] light extinction coefficient in sediment
    ki_PI    = 400,          # [umol/m2/s] optimal light intensity in photosynthesis (was kI)
    cc_PI          = 1,      # [-] parameter to determine the shape of P-I curve (was cc)
    mu_bio        = 0.3e-3,  # [m] the depth where the maximum biomass shows (was u)
    sd_bio        = 0.15e-3, # [m] the standard deviation of mpb biomass distribution (was sd)
    deepODU       = 50,      # [mol/m3] lower boundary concentration of ODU (was bsODU)
    m_photoresp   = 4,       # [-] order of photorespiration vs. O2  (was m)
    prodfacdry    = 0.4,     # [-] mpb production reduction fraction during dry conditions (was aa)
    T0            = 10,      # [dgC] base temperature
    Q10_min       = 2.4,     # [-]   Q10 coefficient mineralisation
    Q10_odureox   = 1,       # [-]   Q10 coefficient reoxidation 
    Q10_mpbprod   = 1.47,    # [-]   Q10 coefficient photosynthesis 
    Q10_basalresp = 1,       # [-]   Q10 coefficient dark respiration
    Q10_photoresp = 1,       # [-]   Q10 coefficient photorespiration
    Tmax_MPBprod  = 42,      # [dgC] Max temperature where photosynthesis = 0
    Tmin_MPBprod  = -15,     # [dgC] Min temperature where photosynthesis = 0
    Topt_MPBprod  = 37,      # [dgC] Optimal temperature for photosynthesis
    Tmax_min      = 44,      # [dgC] Max temperature where mineralisation = 0
    Tmin_min      = 10,      # [dgC] Min temperature where mineralisation = 0
    Topt_min      = 35,      # [dgC] Optimal temperature for mineralisation
    Tmax_odureox  = 55,      # [dgC] Max temperature where reoxidation = 0
    Tmin_odureox  = 0,       # [dgC] Min temperature where reoxidation = 0
    Topt_odureox  = 38,      # [dgC] Optimal temperature for reoxidation
    Tmax_basalresp= 45,      # [dgC] Max temperature where respiration = 0
    Tmin_basalresp= 10,      # [dgC] Min temperature where respiration = 0
    Topt_basalresp= 35,      # [dgC] Optimal temperature for respiration
    Tmax_photoresp= 41,      # [dgC] Max temperature where photorespiration = 0
    Tmin_photoresp= 12,      # [dgC] Min temperature where photorespiration = 0
    Topt_photoresp= 36,      # [dgC] Optimal temperature for photorespiration
    psat_O2   = 1        # [-] saturation fraction of bottom water oxygen (was p_saturation)
  )
  
  .O2DIA$Parunit <- c("m", "m", "-", "-", "m/d", "mol/m3", "mol/m3", 
                     "mol/m3", "m", "molC/m2", "-", "-", "/d", "/d", "/d", 
                     "/m", "/m", "umol/m2/s", "-", "m", "m", "mol/m3", "-", 
                     "-", "dgC", "-", "-", "-", "-", "-", 
                     "dgC", "dgC", "dgC", "dgC", "dgC", "dgC", "dgC", "dgC", 
                     "dgC", "dgC", "dgC", "dgC", "dgC", "dgC", "dgC",  "-") 
  
  
  .O2DIA$Pardesc <- c(
    "thickness of diffusive boundary layer ", 
    "thickness of sediment" ,
    "volumetric porosity" ,
    "salinity",
    "sediment accretion (=advection) velocity",
    "O2 affinity",
    "O2 limitation in ODU oxydation",
    "O2 limitation in dark respiration" ,
    "limitation of GPP by water level in GPP-water level relation" ,
    "MPB biomass",
    "growth respiration fraction (function of photosynthesis)",
    "1-order ODU reoxidation rate at T0",
    "1-order O2 production rate at T0 (photosynthesis)",
    "1-order O2 dark respiration rate at T0",
    "O2 photorespiration fraction (function of photosynthesis)",
    "light extinction coefficient in water",
    "light extinction coefficient in sediment",
    "optimal light intensity in photosynthesis",
    "parameter to determine the shape of P-I curve",
    "the depth of maximum microphytobenthos biomass",
    "the standard deviation of mpb biomass distribution",
    "lower boundary concentration of ODU" ,
    "order of photorespiration vs. O2",
    "fraction of photosynthesis during dry conditions" ,
    "Q10 base temperature",
    "Q10 coefficient growth respiration+mineralisation",
    "Q10 coefficient reoxidation",
    "Q10 coefficient photosynthesis",
    "Q10 coefficient dark respiration",
    "Q10 coefficient photorespiration",
    "Max temperature where photosynthesis = 0",
    "Min temperature where photosynthesis = 0",
    "Optimal temperature for photosynthesis",
    "Max temperature where mineralisation = 0",
    "Min temperature where mineralisation = 0",
    "Optimal temperature for mineralisation",
    "Max temperature where odu reoxidation = 0",
    "Min temperature where odu reoxidation = 0",
    "Optimal temperature for odu reoxidation rate",
    "Max temperature where basal respiration = 0",
    "Min temperature where basal respiration = 0",
    "Optimal temperature for basal respiration",
    "Max temperature where photorespiration = 0",
    "Min temperature where photorespiration = 0",
    "Optimal temperature for photorespiration" ,
    "saturation fraction of bottom water oxygen")


##------------------------------------
## State variables
##------------------------------------

.O2DIA$ynames <- c("O2", "ODU")
.O2DIA$svar <- .O2DIA$ynames

.O2DIA$yunits <- c("molO2/m3 liquid", "molO2/m3 liquid")

.O2DIA$ydescrip <- c("Oxygen (liquid)", 
                     "Oxygen Demand Units (liquid)")

.O2DIA$ynamesall <- as.vector(sapply(.O2DIA$ynames, 
        FUN = function(x) rep(x, times = .O2DIA$N)))

##------------------------------------
## 0-D Variables
##------------------------------------

.O2DIA$var0D <- c(
  "Totalcons",  "Totalmin",  "Totalbasalresp", 
  "Totalphotoresp", "Totalodureox", "Totalprod", 
  "rODUreox", "rMPBprod", "rbasalresp", 
  "O2SWIflux",  "O2deepflux", 
  "ODUSWIflux", "ODUdeepflux")

.O2DIA$unit0D <- rep("molO2/m2/d", times = 13)
      
.O2DIA$descrip0D <- c(
      "Total oxygen consumption", 
      "Total oxygen consumption by growth respiration+mineralisation",
      "Total oxygen consumption by dark (basal) respiration",
      "Total oxygen consumption by photo respiration",
      "Total oxygen consumption by ODU reoxidation",
      "Total oxygen production by photosynthesis",
      "actual rate of ODU reoxidation", 
      "actual rate of MPB production",
      "actual rate of basal respiration",      
      "O2 influx sediment-water", "O2 efflux lower boundary", 
      "ODU influx sediment-water", "ODU efflux lower boundary")
      
##------------------------------------
## forcing functions 
##------------------------------------

.O2DIA$varforc <- c("Temperature", "PAR0", "Waterheight")

.O2DIA$unitforc <- c("dgC", "uEinst/m2/s", "m")

.O2DIA$descripforc <- c( "Sediment temperature", 
      "PAR (Light) on top of sediment", 
      "Bottom water height"
      )
      
##------------------------------------
## 1D variables
##------------------------------------
      
.O2DIA$var1D <- c("O2consumption",  
                   "O2prod",        
                   "O2min",     
                   "O2reox",     
                   "O2basalresp",        
                   "O2photoresp",       
                   "ODUconsumption", 
                   "BioDist", 
                   "Light",
                   "kext", 
                   "PI")

.O2DIA$unit1D <- c( rep("molO2/m3 liquid/d", times = 7), 
                    "molC/m3solid",
                    "uEinst/m2/s", "/m", "-")
      
.O2DIA$descrip1D <- c(      
  "profile of Total oxygen consumption", 
  "profile of oxygen production by photosynthesis",
  "profile of oxygen consumption by mineralisation",
  "profile of oxygen consumption by ODU reoxidation",
  "profile of oxygen consumption by dark (basal) respiration",
  "profile of oxygen consumption by photo respiration",
  "profile of ODU reoxidation",
  "profile of microphytobenthos biomass", 
  "light profile", 
  "extinction coefficient profile",
  "profile of relative photosynthesis rate (light limitation)")
.O2DIA$var1Dall <- as.vector(
  sapply(.O2DIA$var1D, FUN = function(x) rep(x, times = .O2DIA$N)))

.O2DIA$outnames <- c(.O2DIA$var1Dall, .O2DIA$var0D, .O2DIA$varforc)

.O2DIA$nout     <- length(.O2DIA$outnames)

# text used for labeling plots
.O2DIA$labels <- data.frame(
    Units = c("molO/m3 liquid", 
        "molC/m3 solid", 
        "molO2/m2/d", "molO2/m3/d", "uEinst/m2/s"),
    
    Labels = c("mol/m3.l", "mol/m3.s", 
        "mol/m2/d", "mol/m3.l/d", "uE/m2/s")    )

row.names(.O2DIA$labels) <- .O2DIA$labels$Units

.O2DIA$getplot1D <- 
    rbind(data.frame(names = .O2DIA$svar, 
                     units = .O2DIA$labels[.O2DIA$yunits,2]),
          data.frame(names = .O2DIA$var1D, 
                     units = .O2DIA$labels[.O2DIA$unit1D,2]))
