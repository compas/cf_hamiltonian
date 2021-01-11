!***********************************************************************
!                                                                      *
      PROGRAM CF_Hamiltonian
!                                                                      *
!     Crystal field splitting program                                  *
!                                                                      *
!     Written by G. Gaigalas and D. Kato                               *
!                                        NIFS, Japan, September 2013   *
!     The last modification made by G. Gaigalas       June      2018   *
!                                                                      *
!***********************************************************************
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE vast_kind_param, ONLY:  DOUBLE
      USE default_C
      USE iounit_C
      USE CrysPAR_C
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      USE getyn_I
      USE setdbg_I
      USE setmc_I
      USE setcon_I
      USE setsum_I
      USE setcsla_I
      USE gethfd_I
      USE getmixblock_I
      USE strsum_I
      USE factt_I
      USE Ions_input_I
      USE CF_Hamil_I
      IMPLICIT NONE
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      LOGICAL           :: YES
      CHARACTER(LEN=24) :: NAME
      INTEGER           :: K, NCI, NCORE, NCOUNT1
!-----------------------------------------------
!
      print *, " "
      print *, "CF_Hamiltonian: Calculate the crystal-field splitting of ion"
      print *, "                in crystal compound     (Fortran 95 version)"
      print *, "                (C) Copyright by    G. Gaigalas and D. Kato,"
      print *, "                (2019)."
      print *, "                Input  files: isodata, name.c, name.w,      "
      print *, "                              name.(c)m"
      print *, "                Output files: name.(c)CF-Hamil, name.(c)CFm "
      print *, " "
      CALL STARTTIME (NCOUNT1, 'CF_Hamiltonian')

      WRITE (ISTDE, *)
      WRITE (ISTDE, *) 'Default settings?'
      YES = GETYN()
      WRITE (ISTDE, *)
      IF (YES) THEN
         NDEF = 0
      ELSE
         NDEF = 1
      ENDIF
!
! Get name of the state (used in files like <name>.c, <name>.s)
!
      DO WHILE(.TRUE.)
         WRITE (ISTDE, *) 'Name of state'
         READ (*, '(A)') NAME
         K = INDEX(NAME,' ')
         IF (K > 1) EXIT
         WRITE (ISTDE, *) 'Names may not start with a blank'
      END DO
!
      WRITE (ISTDE, *)
      WRITE (ISTDE, *) 'Mixing coefficients from a CI calc.?'
      YES = GETYN()
      IF (YES) THEN
         NCI = 0
      ELSE
         NCI = 1
      ENDIF
!
!   Check compatibility of plant substitutions
!
!GG      CALL CHKPLT
!
!   Determine if there is to be any debug printout; this will be
!   made on the  .dbg  file
!
      CALL SETDBG
!
!   Perform machine- and installation-dependent setup
!
      CALL SETMC
!
!   Set up the physical constants
!
      CALL SETCON
!
!   Open the  .sum  file
!
      CALL SETSUM (NAME, NCI)
!
!   Open, check, load data from, and close, the  .csl  file
!
      CALL SETCSLA (NAME, NCORE)
!
!   Get the remaining information
!
      CALL GETHFD (NAME)
!
!   Get the eigenvectors
!
      CALL GETMIXBLOCK (NAME, NCI)
!
      WRITE (istde,*)
      WRITE (istde,*) 'For which State ?'
      WRITE (istde,*) 'Max ?'
      READ(*,'(I10)') Max_state
      WRITE (istde,*) 'Min ?'
      READ(*,'(I10)') Min_state
      IF(Max_state < Min_state) THEN
        WRITE (istde,*)'The error in input'
        WRITE (istde,*)'The Max_state = ',Max_state,                   &
        '  must be larger then  Min_state = ',Min_state
        STOP
      END IF
!
!   Append a summary of the inputs to the  .sum  file
!
      CALL STRSUM
!
!   Set up the table of logarithms of factorials
!
      CALL FACTT
!
!   Proceed with the HFS calculation
!
      CALL Ions_input
      CALL CF_Hamil(NCORE)
!
!   Print completion message
!
      CALL STOPTIME (NCOUNT1, 'CF_Hamiltonian')
!
      STOP
      END PROGRAM CF_Hamiltonian
