!***********************************************************************
!                                                                      *
      SUBROUTINE index(n,a,ldown,indx)
!                                                                      *
!     Sort out the order of array a and store the index in indx        *
!                                                       (a pointer)    *
!     The input array a is unchanged written in the bases of UpDown    *
!                                                                      *
!                                        NIFS, Japan, September 2013   *
!     The last modification made by G. Gaigalas       June      2018   *
!                                                                      *
!***********************************************************************
!
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE vast_kind_param, ONLY: DOUBLE
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      LOGICAL          :: ldown           ! .TRUE. then Big ---> Small
      INTEGER          :: n, indx(n)
      REAL(DOUBLE)     :: aimx
      REAL(DOUBLE), DIMENSION(N) :: A
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER          i, j, ipos, jpos, jhere
!-----------------------------------------------
!
!     Initialize the index array
      DO i = 1, n
         indx(i) = i
      ENDDO

      IF (ldown) THEN
         DO i = 1, n
            ipos = indx(i)
            aimx = a(ipos)
            jhere = i
            DO j = i+1, n
               jpos = indx(j)
               IF(a(jpos) > aimx) THEN
                  aimx = a(jpos)
                  jhere = j
               ENDIF
            ENDDO
            indx(i) = indx(jhere)
            indx(jhere) = ipos
         ENDDO
      ELSE
         DO i = 1, n
            ipos = indx(i)
            aimx = a(ipos)
            jhere = i
            DO j = i+1, n
               jpos = indx(j)
               IF(a(jpos) < aimx) THEN
                  aimx = a(jpos)
                  jhere = j
               ENDIF
            ENDDO
            indx(i) = indx(jhere)
            indx(jhere) = ipos
         ENDDO
      ENDIF
      RETURN
      END SUBROUTINE index
