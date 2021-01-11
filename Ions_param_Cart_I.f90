      MODULE Ions_param_Cart_I
      INTERFACE
      SUBROUTINE Ions_param_Cart                                       &
                            (ion_number,q_ion,R_ion,GTheta_ion,Phy_ion)
      USE vast_kind_param,        ONLY: DOUBLE
      INTEGER, INTENT(IN)      :: ion_number
      REAL(DOUBLE),INTENT(OUT) :: q_ion, R_ion, GTheta_ion, Phy_ion
      END SUBROUTINE
      END INTERFACE
      END MODULE
