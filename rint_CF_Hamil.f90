!***********************************************************************
!                                                                      *
      REAL(KIND(0.0D0)) FUNCTION RINT_CF_Hamil (K,I,J,RJ)
!                                                                      *
!     The value of RINT_CF is an approximation to:                     *
!          k    k+1                                                    *
!     I ( r  / Rj * ( P (r)*P (r) + Q (r)*Q (r)); 0 to infinity)       *
!          <    >     I     J       I     J                            *
!                                                                      *
!                                                                      *
!     where   I ( G(r) ; Range )  denotes  the  integral  of G(r) over *
!     range.                                                           *
!                                                                      *
!     Call(s) to: [LIB92]: QUAD.                                       *
!                                                                      *
!     Written by G. Gaigalas and D. Kato                               *
!                                        NIFS, Japan, September 2013   *
!     The last modification made by G. Gaigalas       June      2018   *
!                                                                      *
!***********************************************************************
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE vast_kind_param, ONLY: DOUBLE
      USE CONS_C,          ONLY: ZERO
      USE grid_C
      USE tatb_C,          ONLY: MTP, TA
      USE wave_C
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      USE quad_I
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER, INTENT(IN)      :: I, J, K
      REAL(DOUBLE), INTENT(IN) :: RJ
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: L
      REAL(DOUBLE) :: RESULT
!-----------------------------------------------
!
!     Tabulate integrand as required for SUBROUTINE QUAD
      MTP = MIN (MF(I),MF(J))
!
!     Value at first tabulation point is arbitrary
      TA(1) = ZERO
      DO L = 2,MTP
!
!        if r < Rj
         IF (R(L) < RJ) THEN
            TA(L) = (R(L)**K)/(RJ**(K+1))*                             &
                    (PF(L,I)*PF(L,J)+QF(L,I)*QF(L,J))*RP(L)
         ELSE
           TA(L) = (RJ**K)/(R(L)**(K+1))*                              &
                    (PF(L,I)*PF(L,J)+QF(L,I)*QF(L,J))*RP(L)
         END IF
      END DO
!
!     Perform integration
      CALL QUAD (RESULT)
      RINT_CF_Hamil = RESULT
!
      RETURN
      END FUNCTION RINT_CF_Hamil
