      MODULE CrysPAR_C
      USE vast_kind_param, ONLY:  DOUBLE
      CHARACTER(LEN=4)  :: COORD
      CHARACTER(LEN=10) :: SYMMETRY
      REAL(DOUBLE), DIMENSION(:), pointer :: QQ, X1, X2, X3
      REAL(DOUBLE) :: Vector_L
      INTEGER      :: ion_max_number, Max_state, Min_state
      END MODULE CrysPAR_C
