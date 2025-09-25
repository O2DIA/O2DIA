!==========================================================================
! O2DIA, the Diagenetic model of oxygen dynamics 
!
! implemented in FORTRAN
!
! module declarations and common blocks
!
! Karline Soetaert, nioz-yerseke
!==========================================================================
! ==============================================================================
! ==============================================================================
! Module with dimensions
! ==============================================================================
! ==============================================================================

MODULE O2dim

  INTEGER, PARAMETER :: N = 300                       ! number of grid points
  INTEGER, PARAMETER :: nparms = 2*N + 2*(N+1) + 44   ! number of parameters
  INTEGER, PARAMETER :: nforc = 3                     ! number of forcings
  INTEGER, PARAMETER :: nout = 11*N + 13              ! number of output vars

END MODULE O2dim

! ==============================================================================
! Module with common blocks
! ==============================================================================

MODULE O2module
USE O2dim

IMPLICIT NONE

! ------------------------------------------------------------------------------
! 36 parameter values
! ------------------------------------------------------------------------------

  DOUBLE PRECISION ::      &
   por0,         & ! 0.9,      [-] volumetric porosity 
   S,            & ! 31,       [-] salinity  
   v_adv,        & ! 0,        [m/d]   sediment accretion rate 
   ks_min,       & ! 0.00001,  [mol/m3] O2 affinity
   ks_odureox,   & ! 0.00001,  [mol/m3] O2 limitation in ODU oxydation
   ks_basalresp, & ! 0.00001,  [mol/m3]   O2 limitation in dark respiration
   ks_WL,       & ! 0.00001,   [m] ks for waterlevel impact on production
   biomass,     & ! 1,        [molC/m2] MPB integrated biomass
   cmin0,       & ! 0.16,     [-]  proportion of MPB exudate EPS at T0
   rodureox0,   & ! 0.36,     [/d]  1-order O2 consumption rate at T0
   rMPBprod0,   & ! 8.55,     [/d]  1-order O2 production rate at T0 (photosynthesis)
   rbasalresp0, & ! 0.01,     [/d] 1-order O2 dark respiration rate at T0
   cphotoresp0, & ! 0.001,    [-]  proportion O2 photorespiration rate at T0
   kext_w,      & ! 0.1,      [/m] light extinction coefficient in water
   kext_s,      & ! 20e3,     [/m] light extinction coefficient in sediment
   kI_PI,       & ! 100,      [umol/m2/s] optimal light intensity (P-I) photosynthesis
   cc_PI,       & ! 1,        [-]  par to determine the shape of the P-I curve
   mu_bio,      & ! 0.35e-3,  [m]  the depth where the maximum biomass shows
   sd_bio,      & ! 0.23e-3,  [m]  the standard deviation of biomass distribution
   bsODU,       & ! 50,       [mol/m3] lower boundary concentration of ODU
   m_photoresp, & ! 4,        [-]  order of photorespiration vs. O2
   prodfacdry,  & ! 0.4,      [-]  reduction of photsynthesis when dry
   
   T0,          & ! 10,       [dgC] base temperature if 100: use arrhenius
   Q10_min,        & ! 2.4,      [-]   Q10 coefficient mineralisation
   Q10_odureox,    & ! 1,        [-]   Q10 coefficient reoxidation 
   Q10_mpbprod,    & ! 1.47,     [-]   Q10 coefficient photosynthesis 
   Q10_basalresp,  & ! 1,        [-]   Q10 coefficient dark respiration
   Q10_photoresp,  & ! 1,        [-]   Q10 coefficient photorespiration

   Tmax_MPBprod,    & ! 42,            temperature dependency constants
   Tmin_MPBprod,    & ! -15,           ... for microphytobenthos production
   Topt_MPBprod,    & ! 37,
   Tmax_min,        & ! 44,            ... for mineralisation
   Tmin_min,        & ! 10,
   Topt_min,        & ! 35,
   Tmax_odureox,    & ! 55,            ... forodu reoxidation
   Tmin_odureox,    & ! 0,
   Topt_odureox,    & ! 38,
   Tmax_basalresp,  & ! 45,            ... for dark respiration
   Tmin_basalresp,  & ! 10,
   Topt_basalresp,  & ! 35,
   Tmax_photoresp,  & ! 41,            ... for photo respiration
   Tmin_photoresp,  & ! 12,
   Topt_photoresp,  & ! 36,
   p_saturation   ! 1         [-] saturation percentage of bottom water oxygen

! parameter vectors

  DOUBLE PRECISION ::   &
    x(N),               &
    xint(N+1),          &
    dx(N),              &
    dxint(N+1) 
  
  common/xcbpar/por0, S, v_adv, ks_min, ks_odureox, ks_basalresp,     &
      ks_WL, biomass, cmin0, rodureox0, rMPBprod0, rbasalresp0,       &
      cphotoresp0, kext_w, kext_s, kI_PI, cc_PI, mu_bio, sd_bio,      &
      bsODU, m_photoresp, prodfacdry,   & 
      T0, Q10_min, Q10_odureox, Q10_mpbprod, Q10_basalresp,           &
      Q10_photoresp,                                                  &
      Tmax_MPBprod, Tmin_MPBprod, Topt_MPBprod, Tmax_min,             &
      Tmin_min, Topt_min, Tmax_odureox, Tmin_odureox,                 &
      Topt_odureox, Tmax_basalresp, Tmin_basalresp, Topt_basalresp,   &
      Tmax_photoresp,                                                 &
      Tmin_photoresp, Topt_photoresp, p_saturation,                   &
      x, xint, dx, dxint         
 
! ------------------------------------------------------------------------------
! 2 state variables and derivatives
! ------------------------------------------------------------------------------

  DOUBLE PRECISION ::   &
   O2(N),   & !  [molO2/m3] Oxygen concentration
   ODU(N)     !  [molO2/m3] Oxygen demand units concentration

  DOUBLE PRECISION ::   &
   dO2(N),  & !  [molO2/m3/d] ! derivative of O2
   dODU(N)    !  [molO2/m3/d] 

! ------------------------------------------------------------------------------
! 3 forcing functions
! ------------------------------------------------------------------------------

  DOUBLE PRECISION ::   &
   temp,     & ! [dgC]          temperature
   PAR,      & ! [uEinst/m2/s]  photosynthetically active radiation ON TOP OF SEDIMENT
   waterlevel  ! [m]            height of overlying water

  common/xcbforc/      temp, PAR, waterlevel

! ------------------------------------------------------------------------------
! output variables
! ------------------------------------------------------------------------------

  DOUBLE PRECISION ::   &
     O2consumption(N),  & !     # Total O2 consumption rate, a vector
     O2prod(N),         & !           # O2 production rates, a vector
     O2min(N),          & !
     O2odureox(N),      & !
     O2basalresp(N),    & !
     O2photoresp(N),    & !
     ODUconsumption(N), & !
     BioDist(N),        & !       = BioDist,  # distribution of mpb biomass
     Light(N),          & !         
     kext(N),           & !
     PI(N),             & !    [-] relative photosynthesis rate
! for creating budgets
     TotalCons,         & !    [molO2/m2 BULK/d]
     TotalMin,          & !    [molO2/m2 BULK/d]
     Totalbasalresp,    & !
     Totalphotoresp,    & !
     TotaloduReox,      & !    [molO2/m2 BULK/d]
     TotalProd,         & !
     rodureox,          & !    [/d] actual rates at prevailing temperature
     rMPBprod,          & !
     rbasalresp,        & !
     O2SWIflux,         & !
     O2deepflux,        & ! 
     ODUSWIflux,        & !
     ODUdeepflux 

  common/xcbout/ O2consumption, O2prod, O2min, O2odureox, O2basalresp,    &  
     O2photoresp, ODUconsumption, BioDist, Light, kext, PI, TotalCons,    &  
     TotalMin, Totalbasalresp, Totalphotoresp, TotaloduReox, TotalProd,   & 
     rodureox, rMPBprod, rbasalresp,                                      &
     O2SWIflux, O2deepflux, ODUSWIflux, ODUdeepflux

! ------------------------------------------------------------------------------
! ordinary variables
! ------------------------------------------------------------------------------

 DOUBLE PRECISION :: Flux(N+1), Ds(N+1), porfac(N), Dsporfac(N+1)
 DOUBLE PRECISION :: por(N), porInt(N+1)


END MODULE O2module

