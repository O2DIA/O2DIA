!==========================================================================
! O2DIA, the simple oxygen diagenetic model  
! implemented in FORTRAN
!
! Biogeochemistry
!
! Karline Soetaert, nioz-yerseke
!==========================================================================

!==========================================================================
!==========================================================================
! subroutine calculating the rate of change of
! the oxygen model - here the quantities in the common blocks are named
!==========================================================================
!==========================================================================
      
SUBROUTINE o2mod (neq, t, Conc, dConc, yout, ip)

USE O2module
IMPLICIT NONE
      
!......................... declaration section.............................

INTEGER          :: neq, ip(*), i, BCdown
DOUBLE PRECISION :: t, Conc(12*N), dConc(12*N), yout(*)

DOUBLE PRECISION :: dispO2, dispODU, Ds_a
DOUBLE PRECISION :: bwO2, maxO2, cWH, I0
DOUBLE PRECISION :: cmin, cphotoresp
DOUBLE PRECISION :: sumfac

DOUBLE PRECISION :: Tfac, Q10fac, O2solubility, diffO2, diffH2S

!............................ statements ..................................

!     check memory allocated to output variables
      IF (ip(1) < nout)  CALL rexit("nout not large enough") 
      IF (neq  .NE. 2*N)   CALL rexit("number of state variables not ok") 

! -------------------------------------------------------------------------
! from Conc to o2, odu
! -------------------------------------------------------------------------

      DO I = 1, N
        O2(I)  = Conc(I)
        ODU(I) = Conc(  N +I)
      ENDDO

! -------------------------------------------------------------------------
! -------------------------------------------------------------------------
! Transport
! -------------------------------------------------------------------------
! -------------------------------------------------------------------------
 
! ------------------------
! transport of O2  
! ------------------------
! boundary conditions 
      BCdown = 3  ! zerogradient
      maxO2  = O2solubility(temp, S) * 0.21 / 1000d0
      bwO2   =  maxO2 * p_saturation 
      
! diffusion coefficient 
      ! molecular diffusion in water
      dispO2 = diffO2(temp, S)   
      
      ! effective diffusion coefficient in WATER
      Ds = dispO2 * Dsporfac  
      
      ! diffusive coefficient in air 
      Ds_a  = 2.14e-5*3600*24 
      
      ! a factor that varies between 0 (no water) and 1 (waterlevel = Inf)
      cWH = waterlevel**2.0/(waterlevel**2.0+ks_WL**2) 

      DO I = 1, N+1
        IF (xint(I) > 0.d0) EXIT
        Ds(I) = cWH*Ds(I) + (1.d0-cWH)*Ds_a
      END DO

! transport 
      CALL tran1D (N, O2, bwO2, 0.d0, BCdown, v_adv, Ds,                  &
                   porint, por, dx, dxint, Flux, dO2)
      O2SWIflux  = Flux(1) 
      O2deepflux = Flux(N+1)

! ------------------------
! transport of ODU    
! ------------------------
      BCdown = 2  ! imposed concentration
      dispODU = diffH2S(temp)
      Ds      = dispODU*Dsporfac

      ! diffusive coefficient in air 
      DO I = 1, N+1
        IF (xint(I) > 0.d0) EXIT
        Ds(I) = cWH*Ds(I) + (1.d0-cWH)*Ds_a
      END DO
      
      CALL tran1D (N, ODU, 0.d0, bsODU, BCdown, v_adv, Ds,                &
                   porint, por, dx, dxint, Flux, dODU) 
      ODUSWIflux  = Flux(1) 
      ODUdeepflux = Flux(N+1)

! -------------------------------------------------------------------------
! -------------------------------------------------------------------------
! Biogeochemistry
! -------------------------------------------------------------------------
! -------------------------------------------------------------------------

! ------------------------
! temperature dependency  
! ------------------------
      IF (T0 > 50d0) THEN
       cmin       = cmin0       * Tfac(temp, Tmax_min,                     &
                                       Tmin_min,        Topt_min)
       rodureox   = rodureox0   * Tfac(temp, Tmax_odureox,                 &
                                       Tmin_odureox,    Topt_odureox)
       rMPBprod   = rmpbprod0   * Tfac(temp, Tmax_MPBprod,                 &
                                       Tmin_MPBprod,    Topt_MPBprod)
       rbasalresp = rbasalresp0 * Tfac(temp, Tmax_basalresp,               &
                                       Tmin_basalresp,  Topt_basalresp)
       cphotoresp = cphotoresp0 * Tfac(temp, Tmax_photoresp,               &
                                       Tmin_photoresp,  Topt_photoresp)
      ELSE
       cmin       = cmin0       * Q10fac(temp, T0, Q10_min)
       rodureox   = rodureox0   * Q10fac(temp, T0, Q10_odureox)
       rMPBprod   = rmpbprod0   * Q10fac(temp, T0, Q10_mpbprod)
       rbasalresp = rbasalresp0 * Q10fac(temp, T0, Q10_basalresp)
       cphotoresp = cphotoresp0 * Q10fac(temp, T0, Q10_photoresp)
      END IF
! ------------------------
! distribution of microphytobenthos
! ------------------------
    CALL MPBdist(x, dx, mu_bio, sd_bio, biomass, Biodist)
   
! ------------------------
! MPB photosynthesis
! ------------------------
    
    ! the bulk light extinction coefficient [/m]
    kext  =  kext_w*por + kext_s * (1.d0-por) 
    
    ! the light distribution and photosynthesis
    
    ! photosynthesis reduced under dry conditions with factor prodfacdry
    rMPBprod  = rMPBprod  * (prodfacdry+(1.d0-prodfacdry)*cWH)
    
    ! PAR is the photosynthetically active radiation on top of the sediment
    
    I0  = PAR  ! light on top of sediment
    DO I = 1, N
      
      ! Light in the middle of the box
      IF (x(I) >= 0.d0) THEN  ! Light 
        Light(I) = I0 * exp(-kext(i)* dx(I)/2)   ! light in middle   
        I0       = I0 * exp(-kext(i)* dx(I))     ! light at bottom
        PI(I) = 2d0 * (1.d0+cc_PI)*(Light(I)/kI_PI)/                 &
            ((Light(I)/kI_PI)**2.d0 + 2.d0*cc_PI*(Light(I)/kI_PI)+1.d0)
      ELSE
        Light(I) = I0
        PI(I) = 0.d0
      END IF
      
    END DO
    
    O2prod = rMPBprod * BioDist *  PI * porfac
    
! ------------------------
! basal respiration, 
! ------------------------

    O2basalresp = rbasalresp * BioDist * O2/(O2+ks_basalresp) * porfac

! ------------------------
! mineralization
! ------------------------
!   1st-order to production, limited by O2    
    O2min  = cmin * O2prod * O2/(O2+ks_min) 

! ------------------------
! Light respiration
! ------------------------

! when O2 concentration is very high: proportion of O2 production respired
    O2photoresp = cphotoresp * O2prod * (O2/maxO2)**m_photoresp

! ------------------------
! ODU reoxidation
! ------------------------
   DO I = 1, N
    IF (x(I) <= 0) THEN
      O2odureox(I) = 0.d0
    ELSE  
      O2odureox(I) = rodureox * ODU(I) * O2(I)/(O2(I)+ks_odureox)
    END IF 
   END DO

! total consumption 
    O2consumption  = O2basalresp + O2min + O2photoresp + O2odureox
    ODUconsumption = O2odureox

! -------------------------------------------------------------------------
! -------------------------------------------------------------------------
! mass balances
! -------------------------------------------------------------------------
! -------------------------------------------------------------------------

    dO2  = dO2  + O2prod - O2consumption
    dODU = dODU - ODUconsumption
    
! depth-integrated rates:
    Totalbasalresp = 0.d0
    Totalphotoresp = 0.d0
    TotalMin       = 0.d0
    TotaloduReox   = 0.d0
    TotalCons      = 0.d0
    TotalProd      = 0.d0
    
    DO I = 1, N
      sumfac = dx(i)*por(i)
      Totalbasalresp = Totalbasalresp + O2basalresp(I)   * sumfac
      Totalphotoresp = Totalphotoresp + O2photoresp(I)   * sumfac
      TotalMin       = TotalMin       + O2min(I)         * sumfac
      TotalOduReox   = TotalOduReox   + O2odureox(I)     * sumfac
      TotalCons      = TotalCons      + O2consumption(I) * sumfac
      TotalProd      = TotalProd      + O2prod(I)        * sumfac
    END DO
    
! from do2,... to dconc

    DO I = 1, N
       dConc(I     )  =  dO2(I)
       dConc(N   +I)  =  dODU(I) 
    ENDDO

    CALL getout(yout)  ! output in one long vector

RETURN
END SUBROUTINE o2mod

