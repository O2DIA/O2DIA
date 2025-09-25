!==========================================================================
! O2DIA, the oxygen diagenetic model , 
! implemented in FORTRAN
!
! Main functions
!
! Karline Soetaert, nioz-yerseke
!==========================================================================


!==========================================================================
! initialise the common block with parameter values 
!==========================================================================

SUBROUTINE inito2mod (steadyparms)
USE O2dim
       
IMPLICIT NONE
EXTERNAL steadyparms

DOUBLE PRECISION parms(nparms)
COMMON /xcbpar/parms

     CALL steadyparms(nparms, parms)
     CALL initporosity

RETURN
END SUBROUTINE inito2mod

!==========================================================================

SUBROUTINE initporosity
USE O2module
IMPLICIT NONE

INTEGER :: I

! porosity as an exponential function
DO I = 1, N
  por(I) = min(1d0, por0 + (1d0-por0)*exp(-x(I)/0.01e-3)) 
END DO

DO I = 1, N+1
  porInt(I) = min(1d0, por0 + (1d0-por0)*exp(-xInt(I)/0.01e-3)) 
END DO

! factor to go from rate per solid to rate per liquid
DO I = 1, N
  porfac(I) = (1.d0-por(I))/por(I)
END DO

! factor to go from molecular diffusion to effective sediment diffusion
! taking into account the tortuosity
DO I = 1, N+1
  Dsporfac(I) = 1.d0/(1-log(porInt(I)**2d0)) 
END DO

END SUBROUTINE initporosity

!==========================================================================
! Initialise the forcing function common block
!==========================================================================

SUBROUTINE inito2forc (steadyforcs)
USE O2dim
      
IMPLICIT NONE
EXTERNAL steadyforcs

DOUBLE PRECISION forcs(nforc)
COMMON /xcbforc/forcs

    CALL steadyforcs(nforc, forcs)
       
RETURN
END SUBROUTINE inito2forc
