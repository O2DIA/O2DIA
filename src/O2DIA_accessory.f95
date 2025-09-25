!==========================================================================
! O2DIA, the simple oxygen diagenetic model  
! implemented in FORTRAN
!
! functions and subroutines 
!
! Karline Soetaert, nioz-yerseke
!==========================================================================

! ==============================================================================
! ==============================================================================
! Accessory functions
! ==============================================================================
! ==============================================================================

! ==============================================================================
! viscosity of seawater (from marelac package)
! ==============================================================================

DOUBLE PRECISION FUNCTION viscosity(T, S)

IMPLICIT NONE

DOUBLE PRECISION S, T, P  ! salinity, temperature, pressure

P = 1.013253

   viscosity = 1.791 - T * (0.06144 - T * (0.001451 - T*1.6826e-05)) +    &       
        -0.0001529 * P + 8.3885e-08 * P * P + 0.0024727 * S +             &
        +(6.0574e-06 * P - 2.676e-09 * P * P) * T +                       &
        (T * (4.8429e-05 - T * (4.7172e-06 - T * 7.5986e-08))) * S              

END FUNCTION viscosity

! ==============================================================================
! diffusion coefficient of oxygen (from marelac package) 
! ==============================================================================

DOUBLE PRECISION FUNCTION diffO2(T, S)

IMPLICIT NONE
DOUBLE PRECISION S, T, TK, mu_0
DOUBLE PRECISION :: viscosity

   mu_0 = viscosity(T, 0.d0)
   TK   = T + 273.15

   diffO2 = (0.2604 + 0.006383 * (TK/mu_0)) * 1e-09* 86400d0  ! m2/d

END FUNCTION diffO2

! ==============================================================================
! diffusion coefficient of sulphide (from marelac package) 
! ==============================================================================

DOUBLE PRECISION FUNCTION diffH2S(T)

IMPLICIT NONE
DOUBLE PRECISION S, T, P, TK, mu_0
DOUBLE PRECISION :: viscosity

   mu_0 = viscosity(T, 0.d0)
   TK   = T + 273.15

   diffH2S = 4.72e-07 * TK/(mu_0 * 35.2**0.6) * 1e-04*86400d0  ! m2/d

END FUNCTION diffH2S

! ==============================================================================
! solubility of oxygen (from marelac package) 
! ==============================================================================

DOUBLE PRECISION FUNCTION O2solubility(T, S)

IMPLICIT NONE
DOUBLE PRECISION :: T, S
DOUBLE PRECISION A1, A2, A3, A4, B1, B2, B3, bet, TK

   TK = T + 273.15

   A1 = -58.3877  
   A2 =  85.8079  
   A3 =  23.8439  
   A4 =  0.0
   B1 = -0.0348920  
   B2 =  0.0155680 
   B3 = -0.0019387

   bet = A1 + A2 * (100/TK) + A3 * log(TK/100) +                   &
         A4 * (TK/100)**2 + S * (B1 + B2 * TK/100 +                &
         B3 * (TK/100)**2)
   O2solubility = exp(bet)/22.4136 * 1e6/1.013253

END FUNCTION O2solubility

! ==============================================================================
! distribution of microphytobenthos in the sediment
! ==============================================================================

SUBROUTINE MPBdist(x, dx, mu, sd, B, MPB)
USE O2dim

IMPLICIT NONE

DOUBLE PRECISION :: mu, sd  ! mean and standard deviation of profile
DOUBLE PRECISION :: B       ! mean concentration of MPB

DOUBLE PRECISION :: x(N), dx(N)  ! position and thickness of layers
DOUBLE PRECISION :: sumB, tmp, dtot

DOUBLE PRECISION :: MPB(N)
DOUBLE PRECISION, PARAMETER :: pi = 3.141593

INTEGER I

! normal distribution

  sumB = 0.d0
  dtot = 0.d0

  DO I = 1, N

    IF (x(I) >= 0.d0) THEN
  
      tmp = 1/(sqrt(2.d0*pi)*sd)*exp(-(x(I)-mu)**2d0/2d0/sd**2)
    ELSE
  
      tmp = 0.d0
    END IF
  
    sumB = sumB + tmp*dx(I)  
    dtot = dtot + dx(I)   ! do not use this if biomass = integrated
  
    MPB(I) = tmp

  END DO

  ! rescale
  
!  Tmp = B*dtot/sumB
  
  Tmp = B/sumB
  DO I = 1, N
    MPB(I) = MPB(I) * Tmp
  END DO

END SUBROUTINE MPBdist

! ==============================================================================
! Temperature scaling factor
! ==============================================================================

DOUBLE PRECISION FUNCTION Tfac (T, Tmax, Tmin, Topt)
IMPLICIT NONE

DOUBLE PRECISION :: T, Tmax, Tmin, Topt
  
  IF (T < Tmin .or. T > Tmax) THEN
     Tfac = 0.d0
  ELSE     
     Tfac = (T-Tmax)*(T-Tmin)**2.d0/                                     &
        ((Topt-Tmin)*((Topt-Tmin)*(T-Topt)-(Topt-Tmax)*(Topt+Tmin-2.d0*T)))
  END IF

END FUNCTION Tfac


! ==============================================================================
! Temperature scaling factor
! ==============================================================================

DOUBLE PRECISION FUNCTION Q10fac (T, T0, Q10)
IMPLICIT NONE

DOUBLE PRECISION :: T, T0, Q10
  
  Q10fac = Q10 ** ((T-T0)/10d0)

END FUNCTION Q10fac

