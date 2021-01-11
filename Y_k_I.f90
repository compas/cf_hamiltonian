      MODULE Y_K_I
      INTERFACE
      SUBROUTINE Y_K (K,KM,GTheta,Phi,VAL_REAL,VAL_IM)
      USE vast_kind_param,        ONLY: DOUBLE
      INTEGER, INTENT(IN)       :: K, KM
      REAL(DOUBLE), INTENT(IN)  :: GTheta, Phi
      REAL(DOUBLE), INTENT(OUT) :: VAL_REAL, VAL_IM
      END SUBROUTINE
      END INTERFACE
      END MODULE
