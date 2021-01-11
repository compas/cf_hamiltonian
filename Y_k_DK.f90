!***********************************************************************
!                                                                      *
      SUBROUTINE Y_K_DK (K,KM,GTheta,Phi,VAL_REAL,VAL_IM)
!                                                                      *
!     Calculates the spheric function Y^k _q by its algebraic formulae *
!     (see  subsection 5.2  Eq. (1) page 133 in                        *
!     D.A. Varshalovich, A.N. Moskalev and V.K. Khersonskii,           *
!     Quantum Theory of Angular Momentum; Berkeley, CA, 1981).         *
!                                                                      *
!     Written by D. Kato and G. Gaigalas                               *
!                                        NIFS, Japan, September 2013   *
!     The last modification made by G. Gaigalas       June      2018   *
!                                                                      *
!***********************************************************************
!
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE vast_kind_param, ONLY: DOUBLE
      USE FACTS_C
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      USE plgndr_I
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER, INTENT(IN)       :: K, KM
      REAL(DOUBLE), INTENT(IN)  :: GTheta, Phi
      REAL(DOUBLE), INTENT(OUT) :: VAL_REAL, VAL_IM
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
!      INTEGER, PARAMETER :: MFACT = 500
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!----------------------------------------------
      INTEGER      :: KMM, KMP
      REAL(DOUBLE) :: PI, PK, X, FACT0, FACT
!----------------------------------------------
!
      pi=4.d0*datan(1.d0)
      fact0=0.25d0*dble(2*K+1)/pi
      kmm=K-abs(KM)
      kmp=K+abs(KM)
!GG      fact=fact0*factrl(kmm)/factrl(kmp)
      fact=fact0*DEXP(GAM(kmm+1))/DEXP(GAM(kmp+1))
      fact=dsqrt(fact)
!
!      CALL Y_K (K,KM,GTheta,Phi,VAL_REAL,VAL_IM)
!
      x=dcos(GTheta)
      pk=plgndr(K,abs(KM),x)
      VAL_REAL=fact*pk*dcos(KM*Phi)
      VAL_IM=fact*pk*dsin(KM*Phi)
      if (KM < 0) then
        VAL_REAL=VAL_REAL*(1-2*mod(abs(KM),2))
        VAL_IM=VAL_IM*(1-2*mod(abs(KM),2))
      end if
!
      RETURN
      END SUBROUTINE Y_K_DK
