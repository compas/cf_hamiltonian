!***********************************************************************
!                                                                      *
      SUBROUTINE SETSUM(NAME, NCI)
!                                                                      *
!   Open the  .sum  files on stream 24 and 29                          *
!                                                                      *
!   Call(s) to: [LIB92]: LENGTH, OPENFL.                               *
!                                                                      *
!   Written by Gediminas Gaigalas                               2019   *
!                                                                      *
!***********************************************************************
!
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      USE openfl_I
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER, INTENT(IN) :: NCI
      CHARACTER, INTENT(IN) :: NAME*24
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER   :: K, IERR
!GG      CHARACTER :: FILNAM1*256, FILNAM2*256, DEFNAM*11, FORM*11, STATUS*3
      CHARACTER :: FILNAM1*256, FILNAM2*256, FILNAM3*256, DEFNAM*11, FORM*11, STATUS*3
!-----------------------------------------------
!
!   File  CF_Halimtonian.sum  is FORMATTED
!
      K = INDEX(NAME,' ')
      IF (NCI == 0) THEN
         FILNAM1 = NAME(1:K-1)//'.cCF-Hamil'
!GG_Bk         FILNAM2 = NAME(1:K-1)//'.cB-Const'
         FILNAM3 = NAME(1:K-1)//'.cCFm'
      ELSE
         FILNAM1 = NAME(1:K-1)//'.CF-Hamil'
!GG_Bk         FILNAM2 = NAME(1:K-1)//'.B-Const'
         FILNAM3 = NAME(1:K-1)//'.CFm'
      ENDIF
      FORM = 'FORMATTED'
      STATUS = 'NEW'
!
      CALL OPENFL (29, FILNAM1, FORM, STATUS, IERR)
      IF (IERR /= 0) THEN
         WRITE (6, *) 'Error when opening', FILNAM1
         STOP
      ENDIF
!
      STATUS = 'UNKNOWN'
!GG_Bk      CALL OPENFL (30, FILNAM2, FORM, STATUS, IERR)
      IF (IERR /= 0) THEN
         WRITE (6, *) 'Error when opening', FILNAM2
         STOP
      ENDIF
      CALL OPENFL (31,FILNAM3,FORM,STATUS,IERR)
      IF (IERR .NE. 0) THEN
         PRINT *, 'Error when opening',FILNAM3
         STOP
      ENDIF
      WRITE (29,*) 'CF_Hamiltonian: Execution begins ...'
      WRITE (29,*)
!
      RETURN
      END SUBROUTINE SETSUM
