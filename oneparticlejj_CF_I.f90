      MODULE oneparticlejj_CF_I
      INTERFACE
!
      SUBROUTINE ONEPARTICLEJJ_CF(KA,JA,JB,IA1,IA2,VSHELL)
      USE vast_kind_param, ONLY:  DOUBLE
      USE parameter_def,   ONLY:  NNNW
      INTEGER, INTENT(IN)  :: KA, JA, JB
      INTEGER, INTENT(OUT) :: IA1, IA2
      REAL(DOUBLE), DIMENSION(NNNW), INTENT(OUT) :: VSHELL
      END SUBROUTINE
      END INTERFACE
      END MODULE