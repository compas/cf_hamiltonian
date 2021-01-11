      MODULE index_I
      INTERFACE
      SUBROUTINE index(n,a,ldown,indx)
      USE vast_kind_param,  ONLY: DOUBLE
      LOGICAL          :: ldown
      INTEGER          :: n, indx(n)
      REAL(DOUBLE)     :: aimx
      REAL(DOUBLE), DIMENSION(N) :: A
      END SUBROUTINE
      END INTERFACE
      END MODULE
