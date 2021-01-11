!***********************************************************************
!                                                                      *
      REAL(KIND(0.0D0)) FUNCTION WIGNER_3j(JA,JB,JC,MA,MB,MC)
!                                                                      *
!     Calculates the Wigner 3-j symbol by its algebraic formulae as    *
!     displayed in many texts on the theory of angular momentum        *
!     ((5.1) see R. D. Cowan, The Theory of Atomic Structure and       *
!     Spectra; University of California Press, 1981, p. 142).          *
!     The integer arguments ja,... of this function must be the double *
!     of the corresponding quantum numbers in the 3-j symbol, i.e.     *
!     jk = jk' + jk', mk = mk' + mk' in                                *
!                                                                      *
!                          (  ja'   jb'   jc'  )                       *
!                          (                   )                       *
!                          (  ma'   mb'   mc'  )                       *
!                                                                      *
!                                                                      *
!     Written by G. Gaigalas and D. Kato                               *
!                                        NIFS, Japan, September 2013   *
!     The last modification made by G. Gaigalas       June      2018   *
!                                                                      *
!***********************************************************************
!
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE vast_kind_param, ONLY: DOUBLE
      USE CONS_C,          ONLY: ZERO, ONE
      USE FACTS_C
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      USE ittk_I
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER, INTENT(IN) :: JA, JB, JC, MA, MB, MC
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!----------------------------------------------
      INTEGER, DIMENSION(14) :: IK
      INTEGER :: K, KMIN, KMAX, I
      REAL(DOUBLE) :: DELTA, SUMK, SUM
!----------------------------------------------
!
!       Test the triangular condition and that for magnetic quantum numbers
!
      IF(MA+MB+MC /= 0 ) THEN
         WIGNER_3j = ZERO
         RETURN
      ELSE IF(ITTK(JA,JB,JC) == 0) THEN
         WIGNER_3j = ZERO
         RETURN
      ELSE IF(IABS(MA) > JA .OR. IABS(MB) > JB                         &
         .OR. IABS(MC) > JC) THEN
         WIGNER_3j = ZERO
         RETURN
      END IF
      IK( 1) =  JA+JB-JC
      IK( 2) =  JA-JB+JC
      IK( 3) = -JA+JB+JC
      IK( 4) =  JA+JB+JC+2
      IK( 5) =  JA-MA
      IK( 6) =  JA+MA
      IK( 7) =  JB-MB
      IK( 8) =  JB+MB
      IK( 9) =  JC-MC
      IK(10) =  JC+MC
      IK(11) =  JB-JC-MA
      IK(12) =  JA-JC+MB
      IK(13) =  JC-JB+MA
      IK(14) =  JA-JB-MC
      DO I = 1,14
         IF(MOD(IK(I),2) /= 0) THEN
         WRITE(0,*) 'Antras',I,IK(I)
            WIGNER_3j = ZERO
            RETURN
         END IF
         IK(I) = IK(I)/2
      END DO
!
!      Calculate the 3-j delta factor
      DELTA = DEXP(GAM(IK(1)+1))*DEXP(GAM(IK(2)+1))*DEXP(GAM(IK(3)+1)) &
            * DEXP(GAM(IK(5)+1))*DEXP(GAM(IK(6)+1))*DEXP(GAM(IK(7)+1)) &
            * DEXP(GAM(IK(8)+1))*DEXP(GAM(IK(9)+1))*DEXP(GAM(IK(10)+1))&
            /DEXP(GAM(IK(4)+1))
!
!     Find out the intervall of summation  k  and sum up
      KMIN = MAX(0,IK(11),IK(12))
      KMAX = MIN(IK(1),IK(5),IK(8))
      !
      SUM = 0.0
      DO K = KMIN,KMAX
         SUMK = DEXP(GAM(K+1))       * DEXP(GAM(IK(1)-K+1))            &
              * DEXP(GAM(IK( 5)-K+1))* DEXP(GAM(IK(8)-K+1))            &
              * DEXP(GAM(IK(13)+K+1))* DEXP(GAM(K-IK(12)+1))
         SUMK = 1.0/SUMK
         IF(MOD(K,2) == 0) THEN
            SUM = SUM + SUMK
         ELSE
            SUM = SUM - SUMK
         END IF
      END DO
      IF(MOD(IK(14),2) /= 0) THEN
         WIGNER_3j = -SQRT(DELTA) * SUM
      ELSE
         WIGNER_3j =  SQRT(DELTA) * SUM
      END IF
      RETURN
      END FUNCTION WIGNER_3j
