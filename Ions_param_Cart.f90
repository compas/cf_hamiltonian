!***********************************************************************
!                                                                      *
      SUBROUTINE Ions_param_Cart                                       &
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
      USE CONS_C,          ONLY: ZERO, ONE
      USE IOUNIT_C,        ONLY: ISTDE
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
      REAL(DOUBLE) :: PI, X_coord, Y_coord, Z_coord
!-----------------------------------------------
!
      q_ion   = QQ(ion_number)
      PI      = 4.0D00*ATAN (ONE)
!
      X_coord = Vector_L * X1(ion_number)
      Y_coord = Vector_L * X2(ion_number)
      Z_coord = Vector_L * X3(ion_number)
      R_ion   = DSQRT(X_coord*X_coord+Y_coord*Y_coord+Z_coord*Z_coord)
      IF(X_coord == 0.000000D+00 .AND. Y_coord == 0.000000D+00         &
                                 .AND. Z_coord == 0.000000D+00) THEN
         WRITE (29,*)                                                  &
         'Incorect coordinates of ion =',ion_number
         STOP
      ELSE IF(X_coord == 0.000000D+00.AND.Y_coord == 0.000000D+00) THEN
         IF(Z_Coord > 0.0) THEN
            GTheta_ion = ZERO
         ELSE
            GTheta_ion = 180 * PI / 180.D00
         END IF
!GG         GTheta_ion = DACOS(Z_coord/R_ion)
         Phy_ion = ZERO
      ELSE IF(X_coord == 0.000000D+00.AND.Z_coord == 0.000000D+00) THEN
         GTheta_ion = 90 * PI / 180.D00
         IF(Y_Coord > 0.0) THEN
           Phy_ion = 90 * PI / 180.D00
         ELSE
           Phy_ion = -90 * PI / 180.D00
         END IF
!GG         Phy_ion = ATAN2(Y_coord,X_coord)
      ELSE IF(Y_coord == 0.000000D+00.AND.Z_coord == 0.000000D+00) THEN
         GTheta_ion = 90 * PI / 180.D00
         IF(X_Coord > 0.0) THEN
           Phy_ion = ZERO
         ELSE
           Phy_ion = 180.00D00 * PI / 180.D00
         END IF
!GG         Phy_ion = ATAN2(Y_coord,X_coord)
      ELSE IF(X_coord == 0.000000D+00) THEN
         GTheta_ion = DACOS(Z_coord/R_ion)
         IF(Y_Coord > 0.0) THEN
           Phy_ion = 90 * PI / 180.D00
         ELSE
           Phy_ion = -90 * PI / 180.D00
         END IF
!GG         Phy_ion = ATAN2(Y_coord,X_coord)
      ELSE IF(Y_coord == 0.000000D+00) THEN
         GTheta_ion = DACOS(Z_coord/R_ion)
         IF(X_Coord > 0.0) THEN
           Phy_ion = ZERO
         ELSE
           Phy_ion = 180.00D00 * PI / 180.D00
         END IF
!GG         Phy_ion = ATAN2(Y_coord,X_coord)
      ELSE IF(Z_coord == 0.000000D+00) THEN
         GTheta_ion = 90 * PI / 180.D00
         Phy_ion = ATAN2(Y_coord,X_coord)
      ELSE
         GTheta_ion = DACOS(Z_coord/R_ion)
         Phy_ion = ATAN2(Y_coord,X_coord)
      IF(DABS(R_ion*dcos(GTheta_ion)-Z_coord) > 0.0000001)  THEN
         WRITE(istde,*)'ERROR in Ions_param Z Coord'
         WRITE(istde,'(3(3X,E16.9))')X_coord,Y_coord,Z_coord
         WRITE(istde,'(3(3X,E16.9))')                                  &
             R_ion*dsin(GTheta_ion)*dcos(Phy_ion),                     &
             R_ion*dsin(GTheta_ion)*dsin(Phy_ion),                     &
             R_ion*dcos(GTheta_ion)
             STOP
      ELSE IF(DABS(R_ion*dsin(GTheta_ion)*dsin(Phy_ion)-Y_coord)       &
             > 0.0000001)  THEN
         WRITE(istde,*)'ERROR in Ions_param Y Coord'
         WRITE(istde,'(3(3X,E16.9))')X_coord,Y_coord,Z_coord
         WRITE(istde,'(3(3X,E16.9))')                                  &
             R_ion*dsin(GTheta_ion)*dcos(Phy_ion),                     &
             R_ion*dsin(GTheta_ion)*dsin(Phy_ion),                     &
             R_ion*dcos(GTheta_ion)
             STOP
      ELSE IF(DABS(R_ion*dsin(GTheta_ion)*dcos(Phy_ion)-X_coord)       &
             > 0.0000001)  THEN
         WRITE(istde,*)'ERROR in Ions_param X Coord'
         WRITE(istde,'(3(3X,E16.9))')X_coord,Y_coord,Z_coord
         WRITE(istde,'(3(3X,E16.9))')                                  &
             R_ion*dsin(GTheta_ion)*dcos(Phy_ion),                     &
             R_ion*dsin(GTheta_ion)*dsin(Phy_ion),                     &
             R_ion*dcos(GTheta_ion)
             STOP
      END IF
      END IF
!     Pervedame angstremus i atominius vienetus
!
      R_ion      = R_ion/0.52917706
!
!     Pervedame laipsnius i radianus
!
!      GTheta_ion = GTheta_ion*PI / 180.D 00
!      Phy_ion    = Phy_ion   *PI / 180.D 00
!
      RETURN
      END SUBROUTINE Ions_param_Cart
