!***********************************************************************
!                                                                      *
      SUBROUTINE Ions_input
!                                                                      *
!     Written by G. Gaigalas and D. Kato                               *
!                                        NIFS, Japan, September 2013   *
!     The last modification made by G. Gaigalas       June      2018   *
!                                                                      *
!***********************************************************************
!
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE vast_kind_param, ONLY:  DOUBLE
      USE memory_man
      USE IOUNIT_C,        ONLY: ISTDE
      USE CrysPAR_C
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      USE openfl_I
      IMPLICIT NONE
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      CHARACTER(LEN=3)   :: STATUS
      CHARACTER(LEN=9)   :: FORM
      CHARACTER(LEN=10)  :: LAB1, LAB2
      CHARACTER(LEN=20)  :: NAME_CRYSTAL
      CHARACTER(LEN=256) :: FILNAM
      REAL(DOUBLE)       :: Q
      INTEGER            :: IERR, I
!-----------------------------------------------
!
      FILNAM = 'Crystaldata'
      FORM = 'FORMATTED'
      STATUS = 'OLD'
      CALL OPENFL (55,FILNAM,FORM,STATUS,IERR)
      IF (IERR .EQ. 1) THEN
         WRITE(istde,*) 'Error when opening',FILNAM
         STOP
      END IF
!
      READ(55,'(A20,A10,A10,I10,2X,F5.2,5X,A4)')                       &
            NAME_CRYSTAL,LAB1,LAB2,ion_max_number,Vector_L,COORD
!
      CALL ALLOC (QQ,ion_max_number,'QQ','Ions_input')
      CALL ALLOC (X1,ion_max_number,'X1','Ions_input')
      CALL ALLOC (X2,ion_max_number,'X2','Ions_input')
      CALL ALLOC (X3,ion_max_number,'X3','Ions_input')
!
      SYMMETRY = trim(LAB2)
      SYMMETRY =ADJUSTL(SYMMETRY)
      DO I = 1,ion_max_number
         READ(55,'(3X,D13.6,3(3X,D13.6))')X1(I),X2(I),X3(I),Q
         QQ(I) = Q
      END DO
      CLOSE (55)
!
      WRITE (istde,*)
      WRITE (istde,*)
      WRITE (istde,'(A20,A10,A10,I10,2X,F5.2,5X,A4)')                  &
               NAME_CRYSTAL,LAB1,SYMMETRY,ion_max_number,Vector_L,COORD
      WRITE (29,*)
      WRITE (29,*)
      WRITE (29,'(A20,A10,A10,I10,2X,F5.2)')                           &
                   NAME_CRYSTAL,LAB1,SYMMETRY,ion_max_number,Vector_L
      WRITE (31,'(A20,A10,A10,I10,2X,F5.2)')                           &
                   NAME_CRYSTAL,LAB1,SYMMETRY,ion_max_number,Vector_L
      WRITE (istde,*)
      WRITE (istde,*)
      WRITE (istde,*) 'The crystal field is made totaly from',         &
                      ion_max_number,'nearest neighbors.'
      WRITE (29,*)
      WRITE (29,*)    'The crystal field is made totaly from',         &
                      ion_max_number,'nearest neighbors:'
      RETURN
      END SUBROUTINE Ions_input
