!***********************************************************************
!                                                                      *
      SUBROUTINE MATELT_CF_Hamil_NEW (K,I1,I2,APART)
!                                                                      *
!     This routine computes the angular part of the reduced matrix     *
!     elements of the crystal field operator.                          *
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
      USE CONS_C,          ONLY: ZERO
      USE parameter_def,   ONLY: NNNW
      USE orb_C,           ONLY: NAK
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      USE CRE_I
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER, INTENT(IN)       :: K, I1, I2
      REAL(DOUBLE), INTENT(OUT) :: APART
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: KAP1, KAP2, L1, L2
!-----------------------------------------------
!
!   Set KAP1 and KAP2
      KAP1 = NAK(I1)
      KAP2 = NAK(I2)
!
!   Determine the l quantum numbers
      IF (KAP1 > 0) THEN
         L1 =  KAP1
      ELSE
         L1 = -KAP1-1
      ENDIF
!
      IF (KAP2 > 0) THEN
         L2 =  KAP2
      ELSE
         L2 = -KAP2-1
      ENDIF
!
      IF (MOD (L1+K+L2,2) == 0) THEN
!
!   Parity selection rule satisfied
         APART = CRE (IABS(KAP1),(2*K+1)/2,IABS(KAP2))
      ELSE
         APART = ZERO
      ENDIF
      RETURN
      END SUBROUTINE MATELT_CF_Hamil_NEW
