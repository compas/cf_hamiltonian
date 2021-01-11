      MODULE plgndr_I
      INTERFACE
      REAL(KIND(0.0D0)) FUNCTION plgndr(l,m,x)
      USE vast_kind_param, ONLY: DOUBLE
      INTEGER, INTENT(IN)      :: L, M
      REAL(DOUBLE), INTENT(IN) :: x
      END FUNCTION
      END INTERFACE
      END MODULE
