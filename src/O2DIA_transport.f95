!==========================================================================
! O2DIA, the simple oxygen diagenetic model  
! implemented in FORTRAN
!
! Transport
! 
! Karline Soetaert, nioz-yerseke
!==========================================================================
      
!==============================================================================
! Diffusion in a 1-dimensional finite difference grid
! all inputs are vectors
! subroutine from ReacTran in isnt\doc\fortran directory
!==============================================================================

SUBROUTINE tran1d (N, C, Cup, Cdown, BcDown, v, D,                             &
                   VF, VFmid, dx, dxaux, Flux, dC)

IMPLICIT NONE
INTEGER N                  ! length of C

! input
DOUBLE PRECISION C(N)

! Boundary concentrations (used if Bc.. = 2,4), fluxes (used if Bc= 1)
! and convection coeff (used if Bc=4) Cup, Cdown: either conc or flux
DOUBLE PRECISION Cup, Cdown 

! Diffusion,     volume fraction, advection
DOUBLE PRECISION D(N+1), VF(N+1), VFmid(N), v

! grid size, distance from mid to mid
      DOUBLE PRECISION dx(N), dxaux(N+1)

! boundary condition (1 = flux, 2 = conc, 3 = 0-grad)
      INTEGER BcDown

! output: fluxes and rate of change
      DOUBLE PRECISION Flux(N+1), dC(N) 

! locals
      INTEGER I
      DOUBLE PRECISION AVF, Amid, Cbnd

! -------------------------------------------------------------------------------

! Flux - first internal cells

      IF (v >= 0) THEN
        DO I = 2,N
         Flux(I) = -VF(I)*D(I) * (C(I)-C(I-1)) /dxaux(I)                         &
                  + VF(I)*v*C(I-1)
       ENDDO
      ELSE
       DO I = 2,N
        Flux(I) = -VF(I)*D(I) * (C(I)-C(I-1)) /dxaux(I)                         &
                 + VF(I)*v*C(I)
       ENDDO
      ENDIF
      
! Then the outer cells
! upstream boundary
      IF (v >= 0) THEN
        Cbnd = Cup
      ELSE
        Cbnd = C(1)
      ENDIF
      
 ! imposed concentration
        Flux(1) = -VF(1)*D(1) * (C(1)-Cup) /dxaux(1)                            &
     &           + VF(1)*v*Cbnd

! downstream boundary
      IF (v >= 0 .OR. BcDown .eq. 3) THEN
        Cbnd = C(N)
      ELSE
        Cbnd = Cdown
      ENDIF

      IF (BcDown .EQ. 1) THEN
        Flux(N+1) = Cdown

      ELSE IF (BcDown .EQ. 2) THEN
        Flux(N+1) = -VF(N+1)*D(N+1) * (Cdown-C(N)) /dxaux(N+1)                &
     &              + VF(N+1) * v * Cbnd

      ELSE IF (BcDown .EQ. 3) THEN
        Flux(N+1) = VF(N+1) * v * Cbnd

      ENDIF


! Rate of change = negative flux gradient
      DO I = 1,N
        dC(I) = -(Flux(I+1) - Flux(I))/VFmid(I)/dx(I)
      ENDDO

      RETURN
      END SUBROUTINE tran1D

