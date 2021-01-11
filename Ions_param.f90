!***********************************************************************
!                                                                      *
      SUBROUTINE Ions_param(ion_number,q_ion,R_ion,GTheta_ion,Phy_ion)
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
      USE IOUNIT_C,        ONLY: ISTDE
      USE CrysPAR_C,       ONLY: COORD
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER, INTENT(IN)       :: ion_number
      REAL(DOUBLE), INTENT(OUT) :: q_ion, R_ion, GTheta_ion, Phy_ion
!-----------------------------------------------
!
      IF(COORD == "SPHE") THEN
         CALL Ions_param_Sphe(ion_number,q_ion,R_ion,GTheta_ion,Phy_ion)
      ELSE IF(COORD == "CART") THEN
         CALL Ions_param_Cart(ion_number,q_ion,R_ion,GTheta_ion,Phy_ion)
      ELSE
         WRITE(istde,*) 'Error in definition of corrdinates',COORD
         WRITE(istde,*) 'The file Crystaldata'
         STOP
      END IF
!
      RETURN
      END SUBROUTINE Ions_param
