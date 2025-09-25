!==========================================================================
!==========================================================================
! O2DIA, the simple oxygen diagenetic model  
! implemented in FORTRAN
!
! Output
!
! Karline Soetaert, nioz-yerseke
!==========================================================================
!==========================================================================
! put output variables in one vector
!==========================================================================

SUBROUTINE getout(yout)
USE O2dim
      
IMPLICIT NONE
INTEGER :: i
DOUBLE PRECISION :: yout(*), out(nout), forc(3)

COMMON /xcbout  /out
COMMON /xcbforc /forc
      
  DO i = 1, nout
     yout(i) = out (i)
  ENDDO       

  DO i = 1, nforc
     yout(nout+i) = forc (i)
  ENDDO       
 
END SUBROUTINE getout
       
