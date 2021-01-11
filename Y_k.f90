!***********************************************************************
!                                                                      *
      SUBROUTINE Y_K(K,KM,GTheta,Phi,VAL_REAL,VAL_IM)
!                                                                      *
!     Calculates the spheric function Y^k _q by its algebraic formulae *
!     (see  subsection 5.2.2  Eq. (9) page 134 in                      *
!     D.A. Varshalovich, A.N. Moskalev and V.K. Khersonskii,           *
!     Quantum Theory of Angular Momentum; Berkeley, CA, 1981).         *
!                                                                      *
!     Written by G. Gaigalas and E. Gaidamauskas                       *
!                                         Last revision: 08 Sep 2008   *
!                                                                      *
!***********************************************************************
!
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE vast_kind_param, ONLY: DOUBLE
      USE CONS_C,          ONLY: ZERO, ONE
      USE FACTS_C
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
      INTEGER :: I, MAX_REZIS, MIN_REZIS
      REAL(DOUBLE) :: PI, SUM, COEF
!----------------------------------------------
!
      IF (K >= IABS(KM)) THEN
         PI = 4.0D00*DATAN (ONE)
         IF(KM < 0) THEN
            MAX_REZIS = K
            MIN_REZIS = -KM
         ELSE
            MAX_REZIS = K
            MIN_REZIS = 0
         END IF
         SUM = ZERO
         DO I = MIN_REZIS, MAX_REZIS
            IF((2*I+KM) == 0) THEN
               SUM = (-1)**I * DEXP(GAM(K+I+1))/DEXP(GAM(K-I+1))       &
                   / (DEXP(GAM(I+1))*DEXP(GAM(I+KM+1)))                &
                   + SUM
            ELSE
               SUM = (-1)**I * DEXP(GAM(K+I+1))/DEXP(GAM(K-I+1))       &
                   * DSIN(GTheta*0.5)**(2*I+KM)                        &
                   / (DEXP(GAM(I+1))*DEXP(GAM(I+KM+1)))                &
                   + SUM
            END IF
         END DO
         COEF = (-1)**KM*DSQRT((2.0*K+1.0)*DEXP(GAM(K+KM+1))/          &
                (4.0*PI*DEXP(GAM(K-KM+1))))
         IF(KM == 0) THEN
            COEF = COEF
         ELSE
            COEF = COEF/(DCOS(GTheta*0.5))**KM
         END IF
         VAL_REAL = DCOS(KM*Phi)*COEF*SUM
         VAL_IM   = DSIN(KM*Phi)*COEF*SUM
      ELSE
         VAL_REAL = ZERO
         VAL_IM   = ZERO
      END IF
      RETURN
      END SUBROUTINE Y_K
