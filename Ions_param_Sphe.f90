!***********************************************************************
!                                                                      *
      SUBROUTINE Ions_param_Sphe                                       &
                         (ion_number,q_ion,R_ion,GTheta_ion,Phy_ion)
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
      USE CONS_C,          ONLY: ONE
      USE CrysPAR_C,       ONLY: QQ, Vector_L, X1, X2, X3
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER, INTENT(IN)       :: ion_number
      REAL(DOUBLE), INTENT(OUT) :: q_ion, R_ion, GTheta_ion, Phy_ion
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      REAL(DOUBLE) :: PI
!-----------------------------------------------
!
      PI    = 4.0D00*ATAN (ONE)
!
      q_ion      = QQ(ion_number)
      R_ion      = Vector_L*X1(ion_number)
      GTheta_ion = X2(ion_number)
      Phy_ion    = X3(ion_number)
!
!     Pervedame angstremus i atominius vienetus
!
      R_ion      = R_ion/0.52917706
!
!     Pervedame laipsnius i radianus
!
      GTheta_ion = GTheta_ion*PI / 180.D00
      Phy_ion    = Phy_ion*PI / 180.D00
!
      RETURN
      END SUBROUTINE Ions_param_Sphe
