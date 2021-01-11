!***********************************************************************
!                                                                      *
      SUBROUTINE CF_Hamil(NCORE)
!                                                                      *
!     This routine controls the main sequence of routine calls for     *
!     the calculation of the crystal field splitting.                  *
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
      USE vast_kind_param, ONLY: DOUBLE
      USE memory_man_CF
      USE parameter_def,   ONLY: NNNW
      USE CONS_C,          ONLY: ZERO
      USE decide_C
      USE CrysPAR_C
      USE jlabl_C,               LABJ=>JLBR, LABP=>JLBP
      USE prnt_C
      USE decide_C
      USE syma_C
      USE OPT6_C
!GG      USE DEF_C,           ONLY: CVAC, EMPAM, RBCM, AUCM
      USE orb_C
      USE foparm_C
      USE npar_C
      USE iounit_C
      USE eigv_C
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      USE itjpo_I
      USE convrt_I
      USE MATELT_CF_Hamil_I
      USE Ions_param_I
      USE Y_K_DK_I
      USE RINT_CF_Hamil_I
      USE ONESCALAR_I
      USE ONEPARTICLEJJ_CF_I
      USE WIGNER_3j_I
      USE INDEX_I
      USE wghtd5gg_I
      IMPLICIT NONE
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
!     Matrix elements smaller than CUTOFF are not accumulated
      REAL(DOUBLE), PARAMETER :: CUTOFF = 1.0D-12,                     &
                                 EPS = 1.0D-12,                        &
                                 Rydberg = 109737.31534D0,             &
                                 Convert_meV = 27211.6
      INTEGER, PARAMETER :: N_MATRIX = 500, K_MAX_MAX = 20,            &
                            MJ_MAX = 41,    LWMAX = 4*N_MATRIX
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER, INTENT(IN) :: NCORE
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      LOGICAL :: ldown              ! .TRUE. then Big ---> Small
      CHARACTER(LEN=4), DIMENSION(N_MATRIX) :: LEV_J, LABP_J
      CHARACTER(LEN=11) :: CNUM
      REAL(DOUBLE) :: VAL_REAL, VAL_IM, PI, APART, B_Constant_R
      REAL(DOUBLE) :: B_Constant_I, q_ion, R_ion, GTheta_ion, Phi_ion
      REAL(DOUBLE) :: V_LOC_REAL, V_LOC_IM, B_Constant, ELEMNT_I
      REAL(DOUBLE) :: ELEMNT_R, APART_3j, CONTR_R, CONTR_I, HFC_MAX
      REAL(DOUBLE) :: PART_R, PART_I, VR_MAX
      INTEGER :: FFMIN, FFMAX, FF, I, I1, I2, I3, I4, KAP1, KAP2
      INTEGER :: LL1, J1, JSHELL_MIN, K_MIN, K_V, K_Value, Iq_MAX
      INTEGER :: K_MAX_G, iq_V, K_MAX, iq_Value, ion_number, IC, IR
      INTEGER :: ITJPOC, ITJPOR, IDIFF, K_MAX_2, I_SPIN_ANGULAR=1
      INTEGER :: I_COUNT_TNSRJJ, KT, MJ_MAX_C, IAA, LCNUM, IND, IA, IB
      INTEGER :: MJ, MJ_MAX_R, MJ_Value, MJ_S, K, IA_GG, MJS_Value, KK
      INTEGER :: LOC1, LOC2, II, J, JJ, I_MATRICA, J_MATRICA, JJII
      INTEGER :: JSHELL_MAX, JI, I_FIRST, LWORK, IB_GG, J2, INFO, LL2
!
!GG MATRICA
      COMPLEX*16, DIMENSION(N_MATRIX,N_MATRIX) :: RMATRICA, VL, VR
      COMPLEX*16, DIMENSION(N_MATRIX)     :: WR
      COMPLEX*16, DIMENSION(4*N_MATRIX)   :: WORK, RWORK
      REAL(DOUBLE), DIMENSION(K_MAX_MAX*NNNW) :: IA_array, IB_array,   &
                                                 K_T_array,TNSRJJ_array
      REAL(DOUBLE), DIMENSION(N_MATRIX)   :: MJ_AR
!GG     :     WI(N_MATRIX),
      REAL(DOUBLE), DIMENSION(NNNW) :: TSHELL
      INTEGER, DIMENSION(N_MATRIX) :: INDX, NUMBER_VR, LEV, LEV_MJ,    &
                                     LEV_NO

      REAL(DOUBLE),  DIMENSION(:,:,:,:), pointer ::  AMELT_R
      REAL(DOUBLE),  DIMENSION(:,:,:,:), pointer ::  AMELT_I
      REAL(DOUBLE),  DIMENSION(:,:,:), pointer ::  HFC_R
      REAL(DOUBLE),  DIMENSION(:,:,:), pointer ::  HFC_I
      REAL(DOUBLE),  DIMENSION(:), pointer ::  RINTME
!GG MATRICA
!-----------------------------------------------
!
      IF(Min_state > NVEC .OR. Max_state > NVEC)THEN
         WRITE (istde,*) 'The NVEC = ',NVEC
         WRITE (istde,*) 'The Min_state = ',Min_state
         WRITE (istde,*) 'The Max_state = ',Max_state
         WRITE (istde,*) 'Error in CF_Hamil   1'
         STOP
      END IF
!
!   Initialise
      PI = 4.0D00 * DATAN (1.0D00)
      DO I2 = 1,N_MATRIX
         DO I1 = 1,N_MATRIX
            RMATRICA(I1,I2) =  CMPLX(ZERO,ZERO)
         END DO
      END DO
      CALL ALLOC                                                       &
              (AMELT_R,K_MAX_MAX,MJ_MAX,NNNW,NNNW,'AMELT_R','CF_Hamil')
      CALL ALLOC                                                       &
              (AMELT_I,K_MAX_MAX,MJ_MAX,NNNW,NNNW,'AMELT_I','CF_Hamil')
      DO I4 = 1,NNNW
         DO I3 = 1,NNNW
            DO I2 = 1,MJ_MAX
               DO I1 = 1,K_MAX_MAX
                  AMELT_R(I1,I2,I3,I4) =  ZERO
                  AMELT_I(I1,I2,I3,I4) =  ZERO
               END DO
            END DO
         END DO
      END DO
!GG_Bk      WRITE (30,'(A)')' The B^k_q parameters'
!GG_Bk      WRITE (30,*) '------------------------------------------',       &
!GG_Bk                  '----------------------------------'
!GG_Bk      WRITE (30,*) 'Summation over subshells       ',                  &
!GG_Bk                '         The  B^k_q  parameters'
!GG_Bk      WRITE (30,*) '   I        K             k   q',                  &
!GG_Bk                '     Real part         Imaginary part'
!GG_Bk      WRITE (30,*) '                               ',                  &
!GG_Bk                '      (cm^-1)             (cm^-1)'
!GG_Bk      WRITE (30,*) '------------------------------------------',       &
!GG_Bk                  '----------------------------------'
      K_MAX_G = 0
      CALL ALLOC (RINTME,ion_max_number,'RINTME','CF_Hamil')
!
!   Summation over the subshells
      DO I = 1,NW
        KAP1 = NAK(I)
!
!   Determine the 2*j quantum numbers
        IF(KAP1 >= 0) THEN
           LL1 = KAP1
        ELSE
           LL1 = -KAP1-1
        END IF
        J1 = 2*LL1 - KAP1/IABS(KAP1)
        IF(NCORE >= I) THEN
           JSHELL_MIN = I
           JSHELL_MAX = I
        ELSE
           JSHELL_MIN = NCORE +1
           JSHELL_MAX = NW
        END IF
        DO J = JSHELL_MIN,JSHELL_MAX
!GG        DO J = 1,NW
          KAP2 = NAK(J)
          WRITE (istde,'(A,2I3)')'Summation over the subshells = ',I,J
!GG_Bk          WRITE (30,'(2X,I3,6X,I3)')I,J
!
!   Determine the 2*j quantum numbers
          IF(KAP2 >= 0) THEN
             LL2 = KAP2
          ELSE
             LL2 = -KAP2-1
          END IF
          J2 = 2*LL2 - KAP2/IABS(KAP2)
          K_MAX = (J1+J2)/2+1
          IF(K_MAX_G < K_MAX) K_MAX_G = K_MAX
          K_MIN = IABS(J1-J2)/2+1
          WRITE (istde,'(10X,4(A,I2))')'2*J1=',J1,'  2*J2=',J2,        &
                          '  K_MIN=',K_MIN-1,'  K_MAX=',K_MAX-1
          WRITE (istde,'(A)')''
          IF(K_MAX > K_MAX_MAX) then
             WRITE (istde,*) '1 K_MAX > K_MAX_MAX',K_MAX
             STOP
          END IF
!
!   Summation over the k
          DO K_V = K_MIN, K_MAX
             K_Value = K_V - 1
             IF(SYMMETRY == 'C2') THEN
                IF(K_Value == 0 ) CYCLE
                IF (MOD(K_Value,2) /= 0) CYCLE
             END IF
             CALL MATELT_CF_Hamil (K_Value,I,J,APART)
!GG             CALL MATELT_CF_Hamil_NEW (K_Value,I,J,APART)
             Iq_MAX = 2*K_Value + 1
             IF (DABS (APART) > EPS) THEN
!
!   Summation over the q
                DO iq_V = 1, Iq_MAX
                   iq_Value = iq_V - K_Value - 1
                   IF(SYMMETRY == 'C2') THEN
                      IF(MOD(iq_Value,2) /= 0 ) CYCLE
                   END IF
                   B_Constant_R = ZERO
                   B_Constant_I = ZERO
!
!   Summation over the ions
                   DO ion_number = 1, ion_max_number
                      CALL Ions_param                                  &
                      (ion_number,q_ion,R_ion,GTheta_ion,Phi_ion)
                      CALL Y_K_DK(K_Value,-iq_Value,                   &
                          GTheta_ion,Phi_ion,VAL_REAL,VAL_IM)
                      IF(MOD(K_Value,2) /= 0) THEN
                          V_LOC_REAL = VAL_IM
                          V_LOC_IM   = VAL_REAL
                          VAL_REAL   = V_LOC_REAL
                          VAL_IM     = V_LOC_IM
                         IF(MOD(K_Value+1,4) /= 0) VAL_REAL=-VAL_REAL
                         IF(MOD(K_Value+1,4) /=0) VAL_IM=-VAL_IM
                      ELSE
                         IF(MOD(K_Value,4) /= 0) VAL_REAL=-VAL_REAL
                         IF(MOD(K_Value,4) /= 0) VAL_IM=-VAL_IM
                      END IF
                      IF(iq_V == 1) THEN
                         RINTME(ion_number)=                           &
                             RINT_CF_Hamil (K_Value,I,J,R_ion)
                      END IF
!
!     B^q_k constant begining
                      B_Constant                                       &
!GG     :                   + ((-1)**(K_value-iq_Value+1))*DBLE(q_ion)
!GG     :                    = ((-1)**(iq_Value+1)) * DBLE(q_ion)
!GG     :                    = ((-1)**(iq_Value))*DBLE(q_ion)
!GG     :                    = ((-1)**(K_Value+iq_Value+1))*DBLE(q_ion)
!GG GERAS     :                    = ((-1)**(iq_Value)) * DBLE(q_ion)
!GG V9     :                    = ((-1)**(K_Value-iq_Value+1)) * DBLE(q_ion)
!GG CS     :                    = ((-1)**(iq_Value)) * DBLE(q_ion)
                          = ((-1)**(K_Value-iq_Value+1)) * DBLE(q_ion) &
                          * SQRT((4.0*PI)/(2.0*K_Value+1.0))
                      B_Constant_R = B_Constant_R                      &
                                + B_Constant*RINTME(ion_number)*VAL_REAL
                      B_Constant_I = B_Constant_I                      &
                                + B_Constant*RINTME(ion_number)*VAL_IM
!
!     B^q_k constant end
                   END DO
!
!GG_Bk                   WRITE (30,'(25X,I3,1X,I3,3X,E15.7,3X,E15.7)')       &
!GG_Bk                   K_Value,iq_Value,                                   &
!GG_Bk                   B_Constant_R*Rydberg*2,                             &
!GG_Bk                   B_Constant_I*Rydberg*2
!
                   AMELT_R(K_V,iq_V,I,J) = APART * B_Constant_R
                   AMELT_I(K_V,iq_V,I,J) = APART * B_Constant_I
                   WRITE (istde,'(A,I3,A,I3,A,E15.7,A,E15.7)')         &
                   'K_Value =',K_Value,'  iq_Value =',iq_Value,        &
                   '  AMELT_R =',                                      &
                   AMELT_R(K_V,iq_V,I,J)*Rydberg*2,                    &
                   '  AMELT_I =',                                      &
                   AMELT_I(K_V,iq_V,I,J)*Rydberg*2
                END DO
             END IF
          END DO
        END DO
      END DO
      CALL DALLOC (RINTME,'RINTME','CF_Hamil')
      CALL DALLOC (QQ,'QQ','CF_Hamil')
      CALL DALLOC (X1,'X1','CF_Hamil')
      CALL DALLOC (X2,'X2','CF_Hamil')
      CALL DALLOC (X3,'X3','CF_Hamil')
      WRITE (istde,*)'Summation over the subshells was completed.'
!
!   Set the parity of the one-body operators
!GG      IPT = 1
!
!   Sweep through the Hamiltonian matrix to determine the
!   diagonal and off-diagonal matrix elements
      WRITE (29,101) CUTOFF,NVEC
      WRITE (29,102)
!
!   Allocate storage for local arrays
      CALL ALLOC (HFC_R,MJ_MAX,MJ_MAX,NVEC*NVEC,'HFC_R','CF_Hamil')
      CALL ALLOC (HFC_I,MJ_MAX,MJ_MAX,NVEC*NVEC,'HFC_I','CF_Hamil')
!
!   Initialise
      DO I3 = 1,NVEC*NVEC
         DO I2 = 1,MJ_MAX
            DO I1 = 1,MJ_MAX
               HFC_R(I1,I2,I3) = ZERO
               HFC_I(I1,I2,I3) = ZERO
            END DO
         END DO
      END DO
!
!   Summation over the first configuration
      DO IC = 1,NCF
        ITJPOC = ITJPO (IC)
!
!   Output IC on the screen to show how far the calculation has preceede
        IF (MOD(IC,100) == 0) THEN
           CALL CONVRT (IC,CNUM,LCNUM)
           WRITE (istde,*)'Column '//CNUM(1:LCNUM)//' complete;'
        ENDIF
!
!   Summation over the second configuration
        DO 22 IR = 1,NCF
          ITJPOR = ITJPO (IR)
          IDIFF = ITJPOC - ITJPOR
!GG          IF((IDIFF .EQ. 0) .AND. (IR. GE. IC)
!GG     :                      .OR. IABS(IDIFF) .GT. 0) THEN
            K_MAX_2 = (ITJPOC-1+ITJPOR-1)/2 + 1
            K_MAX = MIN(K_MAX_G,K_MAX_2)
            IF(K_MAX > K_MAX_MAX) then
               WRITE (istde,*) '2 K_MAX > K_MAX_MAX',K_MAX
               STOP
            END IF
            K_MIN = IABS(ITJPOC-1-ITJPOR+1)/2 + 1
            DO I_SPIN_ANGULAR=1,K_MAX_MAX*NNNW
               TNSRJJ_array(I_SPIN_ANGULAR) = ZERO
               IA_array(I_SPIN_ANGULAR)  = 0
               IB_array(I_SPIN_ANGULAR)  = 0
               K_T_array(I_SPIN_ANGULAR) = -100
            END DO
            I_COUNT_TNSRJJ = 0
            DO K_V = K_MIN, K_MAX
               KT = K_V - 1
               IF(SYMMETRY == 'C2') THEN
                  IF(KT == 0 ) CYCLE
                  IF (MOD(KT,2) /= 0) CYCLE
               END IF
               DO IND =1,NW
                  TSHELL(IND) = ZERO
               END DO
               IF( KT == 0) THEN
                  CALL ONESCALAR(IC,IR,IA,IB,TSHELL)
               ELSE
                  CALL ONEPARTICLEJJ_CF(KT,IC,IR,IA,IB,TSHELL)
               END IF
               IF (IA == IB) THEN
                  DO IAA = 1,NW
                    IF (DABS (TSHELL(IAA)) > EPS) THEN
                       I_COUNT_TNSRJJ = I_COUNT_TNSRJJ + 1
                       IF(I_COUNT_TNSRJJ > K_MAX_MAX*NNNW) THEN
                        WRITE(istde,*)'Error: CF_Hamil I_COUNT_TNSRJJ 1'
                        STOP
                       END IF
                       TNSRJJ_array(I_COUNT_TNSRJJ) = TSHELL(IAA)
                       IA_array(I_COUNT_TNSRJJ)  = IAA
                       IB_array(I_COUNT_TNSRJJ)  = IAA
                       K_T_array(I_COUNT_TNSRJJ) = KT
                    END IF
                  END DO
               ELSE
                 IF (DABS (TSHELL(1)) > EPS) THEN
                   I_COUNT_TNSRJJ = I_COUNT_TNSRJJ + 1
                   IF(I_COUNT_TNSRJJ > K_MAX_MAX*NNNW) THEN
                    WRITE(istde,*) 'Error: CF_Hamil I_COUNT_TNSRJJ 2'
                    STOP
                   END IF
                   TNSRJJ_array(I_COUNT_TNSRJJ) = TSHELL(1)
                   IA_array(I_COUNT_TNSRJJ)  = IA
                   IB_array(I_COUNT_TNSRJJ)  = IB
                   K_T_array(I_COUNT_TNSRJJ) = KT
                 ENDIF
               ENDIF
            END DO
!GG            IF (I_COUNT_TNSRJJ .EQ. 0 ) GO TO 22
            IF (I_COUNT_TNSRJJ == 0 ) CYCLE
            MJ_MAX_C = 2*ITJPOC
            MJ_MAX_R = 2*ITJPOR
            IF(MJ_MAX_C > 2*MJ_MAX) THEN
              WRITE (istde,*) 'MJ_MAX_C > 2*MJ_MAX',MJ_MAX_C
              STOP
            END IF
            IF(MJ_MAX_R > 2*MJ_MAX) THEN
              WRITE (istde,*) 'MJ_MAX_R > 2*MJ_MAX',MJ_MAX_R
              STOP
            END IF
!
!   Summation over the M_J ir M_J'
          DO MJ = 2,MJ_MAX_C,2
             MJ_Value = (MJ-ITJPOC-1)
             DO MJ_S = 2,MJ_MAX_R,2
                MJS_Value = (MJ_S-ITJPOR-1)
                iq_Value = (MJ_Value-MJS_Value)/2
                IF(IABS(iq_Value)+1 > K_MAX) CYCLE
                IF(SYMMETRY == 'C2') THEN
                   IF(MOD(iq_Value,2) /= 0 ) CYCLE
                END IF
                ELEMNT_I = ZERO
                ELEMNT_R = ZERO
                DO K_V = IABS(iq_Value)+1, K_MAX
                   KT = K_V - 1
                   IF(SYMMETRY == 'C2') THEN
                     IF(KT == 0 ) CYCLE
                     IF (MOD(KT,2) /= 0) CYCLE
                   END IF
                   iq_V = iq_Value + KT + 1
!
!   Summation over the subshells
                   APART_3j=WIGNER_3j(ITJPOC-1, 2*KT, ITJPOR-1,        &
                                -MJ_Value, 2*iq_Value, MJS_Value)
                   IF (DABS (APART_3j) < EPS) CYCLE
                   APART_3j=DSQRT(ITJPOC*1.0D00)*APART_3j
                   IF(IABS(MOD((ITJPOC-1-MJ_Value)/2,2)) /= 0)         &
                                       APART_3j = -APART_3j
!
!   Accumulate the contribution from the one-body operators;
                   DO I_SPIN_ANGULAR = 1,I_COUNT_TNSRJJ
                      IF(K_T_array(I_SPIN_ANGULAR) /= KT) CYCLE
                      IA_GG = IA_array(I_SPIN_ANGULAR)
                      IB_GG = IB_array(I_SPIN_ANGULAR)
                      ELEMNT_R = ELEMNT_R                              &
                               + AMELT_R(KT+1,iq_V,IA_GG,IB_GG)        &
                               * APART_3j                              &
                               * TNSRJJ_array(I_SPIN_ANGULAR)
                      ELEMNT_I = ELEMNT_I                              &
                               + AMELT_I(KT+1,iq_V,IA_GG,IB_GG)        &
                               * APART_3j                              &
                               * TNSRJJ_array(I_SPIN_ANGULAR)
                   END DO
                END DO
!
!   Multiply with the configuration expansion coefficients and add the
!   contributions from the matrix elements to obtain total contributions
                IF(DABS(ELEMNT_R) > EPS.OR.DABS(ELEMNT_I) > EPS)THEN
                  DO K = Min_state, Max_state
                    LOC1 = (K-1)*NCF
                    IF((IDIFF == 0) .AND. (IR /= IC)) THEN
                      IF (DABS(EVEC(IC+LOC1)) < CUTOFF .AND.        &
                          DABS(EVEC(IR+LOC1)) <  CUTOFF) CYCLE
                    ELSE
                      IF (DABS(EVEC(IC+LOC1)) < CUTOFF) CYCLE
                    END IF
                    DO KK = Min_state, Max_state
                      LOC2 = (KK-1)*NCF
                      IF((IDIFF == 0) .AND. (IR /= IC)) THEN
                        IF (DABS(EVEC(IR+LOC2)) < CUTOFF .AND.         &
                            DABS(EVEC(IC+LOC2)) < CUTOFF) CYCLE
                      ELSE
                        IF (DABS(EVEC(IR+LOC2)) < CUTOFF) CYCLE
                      END IF
!GG                      IF((IDIFF .EQ. 0) .AND. (IR /= IC)) THEN
!GG                      IF((IDIFF .EQ. 0) .AND. (IR /= IC)
!GG     : THEN                             .OR. IABS(IDIFF) .GT. 0) THEN
!GG                        CONTR_R = ELEMNT_R
!GG     :                          *(EVEC(IC+LOC1) * EVEC(IR+LOC2)
!GG     :                          + EVEC(IR+LOC1) * EVEC(IC+LOC2))
!GG                        CONTR_I = ELEMNT_I
!GG     :                          *(EVEC(IC+LOC1) * EVEC(IR+LOC2)
!GG     :                          + EVEC(IR+LOC1) * EVEC(IC+LOC2))
!GG                     ELSE
                        CONTR_R = ELEMNT_R                             &
                                * EVEC(IC+LOC1) * EVEC(IR+LOC2)
                        CONTR_I = ELEMNT_I                             &
                                * EVEC(IC+LOC1) * EVEC(IR+LOC2)
!GG                      ENDIF
                      HFC_R(MJ/2,MJ_S/2,NVEC*(K-1)+KK)                 &
                              = HFC_R(MJ/2,MJ_S/2,NVEC*(K-1)+KK)       &
                              + CONTR_R
                      HFC_I(MJ/2,MJ_S/2,NVEC*(K-1)+KK)                 &
                              = HFC_I(MJ/2,MJ_S/2,NVEC*(K-1)+KK)       &
                              + CONTR_I
                    END DO
                  END DO
                END IF
             END DO
          END DO
!GG          END IF
   22   CONTINUE
      END DO
      CALL DALLOC (AMELT_R,'AMELT_R','CF_Hamil')
      CALL DALLOC (AMELT_I,'AMELT_I','CF_Hamil')
!
!   Printouts
      HFC_MAX = ZERO
      I_MATRICA = 0
      DO I = Min_state, Max_state
        JJ = IATJPO(I)
        MJ_MAX_C = 2*JJ
        WRITE (29,"(/A)") "The eigenvalues before CF are:"
        WRITE (29,"(/A)") "                     LABJ   IVEC"
        WRITE(29,"(D20.10,1X,A4,I5)") EVAL(I)+EAV,LABJ(JJ),IVEC(I)
        DO MJ = 2,MJ_MAX_C,2
          I_MATRICA = I_MATRICA + 1
          IF(I_MATRICA > N_MATRIX) THEN
             WRITE (istde,*) 'I_MATRICA > N_MATRIX N_MATRIX = ',       &
                              N_MATRIX
             STOP
          END IF
          MJ_Value = (MJ-JJ-1)
!GG MJ_AR(N_MATRIX)
          J_MATRICA = 0
          DO II = Min_state, Max_state
            JJII = IATJPO(II)
            MJ_MAX_R = 2*JJII
            DO MJ_S = 2,MJ_MAX_R,2
              J_MATRICA = J_MATRICA + 1
              MJS_Value = (MJ_S-JJII-1)
              IF (DABS (HFC_R(MJ/2,MJ_S/2,NVEC*(I-1)+II)) < EPS)THEN
                 HFC_R(MJ/2,MJ_S/2,NVEC*(I-1)+II) = ZERO
              END IF
              IF (DABS (HFC_I(MJ/2,MJ_S/2,NVEC*(I-1)+II)) < EPS)THEN
                 HFC_I(MJ/2,MJ_S/2,NVEC*(I-1)+II) = ZERO
              END IF
              IF (I_MATRICA == J_MATRICA) THEN
                 LEV_J(I_MATRICA) = LABJ(JJII)
                 LEV_MJ(I_MATRICA)= MJS_Value
                 LEV_NO(I_MATRICA) = II
                 LEV(I_MATRICA)   = IVEC(II)
                 LABP_J(I_MATRICA) = LABP((IASPAR(II)+3)/2)
              END IF
              IF (DABS(HFC_R(MJ/2,MJ_S/2,NVEC*(I-1)+II)) > EPS         &
                  .OR.                                                 &
                  DABS(HFC_I(MJ/2,MJ_S/2,NVEC*(I-1)+II)) > EPS)THEN
                 IF (I_MATRICA == J_MATRICA) THEN
                    IF(MJ_Value /= MJS_Value) THEN
                       WRITE (istde,*) 'Problems in CF_Hamil'
                       WRITE (istde,*)                                 &
                       'MJ_Value,MJS_Value',MJ_Value,MJS_Value
                      STOP
                    END IF
!GG                    PART_R =HFC_R(MJ/2,MJ_S/2,NVEC*(I-1)+II)+EVAL(I)+EAV
                    PART_R =HFC_R(MJ/2,MJ_S/2,NVEC*(I-1)+II)+EVAL(I)
                    PART_I =HFC_I(MJ/2,MJ_S/2,NVEC*(I-1)+II)
                    RMATRICA(I_MATRICA,J_MATRICA) = CMPLX(PART_R,PART_I)
                    WRITE (29,113) IVEC(I),LABJ(JJ),                   &
                             LABP((IASPAR(I)+3)/2),                    &
                             0.5*MJ_Value,                             &
                             IVEC(II),LABJ(JJII),                      &
                             LABP((IASPAR(II)+3)/2),                   &
                             0.5*MJS_Value,                            &
                             HFC_R(MJ/2,MJ_S/2,NVEC*(I-1)+II),         &
                             HFC_I(MJ/2,MJ_S/2,NVEC*(I-1)+II),         &
                             DBLE(RMATRICA(I_MATRICA,J_MATRICA))
!GG                    LEV_J(I_MATRICA) = LABJ(JJ)
!GG                    LEV_MJ(I_MATRICA)= MJ_Value
!GG                    LEV_NO(I_MATRICA) = I
!GG                    LEV(I_MATRICA)   = IVEC(I)
                 ELSE
                    PART_R = HFC_R(MJ/2,MJ_S/2,NVEC*(I-1)+II)
                    PART_I = HFC_I(MJ/2,MJ_S/2,NVEC*(I-1)+II)
                    RMATRICA(I_MATRICA,J_MATRICA) = CMPLX(PART_R,PART_I)
                    WRITE (29,103) IVEC(I),LABJ(JJ),                   &
                             LABP((IASPAR(I)+3)/2),                    &
                             0.5*MJ_Value,                             &
                             IVEC(II),LABJ(JJII),                      &
                             LABP((IASPAR(II)+3)/2),                   &
                             0.5*MJS_Value,                            &
                             HFC_R(MJ/2,MJ_S/2,NVEC*(I-1)+II),         &
                            HFC_I(MJ/2,MJ_S/2,NVEC*(I-1)+II)
                 END IF
                 IF(DABS(HFC_I(MJ/2,MJ_S/2,NVEC*(I-1)+II))             &
                                            > DABS(HFC_MAX)) THEN
                   HFC_MAX = HFC_I(MJ/2,MJ_S/2,NVEC*(I-1)+II)
                 END IF
              END IF
            END DO
          END DO
        END DO
      END DO
!
      CALL DALLOC (HFC_R,'HFC_R','CF_Hamil')
      CALL DALLOC (HFC_I,'HFC_R','CF_Hamil')
!
      WRITE (29,"(/A,D20.10)")                                         &
       "The max value of imaginary parts is ", HFC_MAX
!GG      CALL DGEEV('N','V',N_MATRIX,RMATRICA,N_MATRIX,
!GG     :        WR,WI,VL,N_MATRIX,VR,N_MATRIX,WORK,4*N_MATRIX,INFO)
!GG      CALL ZGEEV('N','V',N_MATRIX,RMATRICA,N_MATRIX,
!GG     :      WR,VL,N_MATRIX,VR,N_MATRIX,WORK,4*N_MATRIX,
!GG     :      RWORK,INFO)
!
!     Query the optimal workspace.
!
      LWORK = -1
      CALL ZGEEV('N','V',N_MATRIX,RMATRICA,N_MATRIX,                   &
            WR,VL,N_MATRIX,VR,N_MATRIX,WORK,LWORK,                     &
            RWORK,INFO)
!
!     Solve eigenproblem.
!
      LWORK = MIN(LWMAX,INT(WORK(1)))
      CALL ZGEEV('N','V',N_MATRIX,RMATRICA,N_MATRIX,                   &
            WR,VL,N_MATRIX,VR,N_MATRIX,WORK,LWORK,                     &
            RWORK,INFO)
!
!     Check for convergence.
!
      IF(INFO /= 0 ) THEN
         WRITE(29,"(/A)") "The algorithm failed to compute eigenvalues"
         WRITE(29,"(/A)") "in the program ZGEEV. "
         WRITE(29,"(/A,I3)") "INFO = ",INFO
         STOP
      END IF
      DO IND = 1,I_MATRICA
         VR_MAX = ZERO
         NUMBER_VR(IND) = 0
         DO J = 1,I_MATRICA
            IF (DABS(REAL(VR(J,IND))) > DABS(VR_MAX)) THEN
               VR_MAX = REAL(VR(J,IND))
               NUMBER_VR(IND) = J
            END IF
         END DO
      END DO
      WRITE (29,"(/A)") "The eigenvalues are:"
      DO IND = 1,I_MATRICA
         WRITE(29,"(D20.10,2X,D20.10)")DBLE(WR(IND)),AIMAG(WR(IND))
      END DO
!
!GG      ldown = .TRUE.
      ldown = .FALSE.
      CALL INDEX(I_MATRICA,DBLE(WR),ldown,INDX)
!
      WRITE (29,"(/A)") 'Energy levels for ...'
      WRITE (29,105) Rydberg
      WRITE (29,*) 'No - Serial number of the state; ',                &
                  'Pos - Position of the state within the '
      WRITE (29,*) 'J/P block; Splitting is the energy difference ',   &
                  'from lowest level with same J values'
      WRITE (29,*) '------------------------------------------------', &
                   '------------------------------------'
      WRITE (29,*) 'No Pos  J     2*M_J     Energy Total      Levels', &
                   '      Levels   Splitting   Splitting'
      WRITE (29,*) '                              (a.u.)     (cm^-1)', &
                   '       (meV)     (cm^-1)       (meV)'
      WRITE (29,*) '------------------------------------------------', &
                   '------------------------------------'
!GG      IF(ldown == .TRUE.) then
      IF(ldown .EQV.  .TRUE.) then
        J = 1
        IND = INDX(J)
        I_FIRST = 1
        WRITE (29,'(I3,1X,I3,1X,A4,A4,I3,4X,F14.7,F12.2,F12.2)')       &
            NUMBER_VR(IND),LEV(NUMBER_VR(IND)),LEV_J(NUMBER_VR(IND)),  &
            LABP_J(NUMBER_VR(IND)),                                    &
            LEV_MJ(NUMBER_VR(IND)),DBLE(WR(IND))+EAV
        DO J = 2,I_MATRICA
           IND = INDX(J)
           IF(LEV_J(NUMBER_VR(IND)) /= LEV_J(NUMBER_VR(INDX(J-1)))     &
           .OR.                                                        &
           LEV(NUMBER_VR(IND)) /= LEV(NUMBER_VR(INDX(J-1)))) THEN
              I_FIRST = IND
              WRITE (29,'(I3,1X,I3,1X,A4,A4,I3,4X,F14.7,2F12.2)')      &
              NUMBER_VR(IND),LEV(NUMBER_VR(IND)),LEV_J(NUMBER_VR(IND)),&
              LABP_J(NUMBER_VR(IND)),                                  &
              LEV_MJ(NUMBER_VR(IND)),DBLE(WR(IND))+EAV,                &
              (DBLE(WR(INDX(1)))-DBLE(WR(IND)))*Rydberg*2,             &
              (DBLE(WR(INDX(1)))-DBLE(WR(IND)))*Convert_meV
           ELSE
              WRITE (29,'(I3,1X,I3,1X,A4,A4,I3,4X,F14.7,4F12.2)')      &
              NUMBER_VR(IND),LEV(NUMBER_VR(IND)),LEV_J(NUMBER_VR(IND)),&
              LABP_J(NUMBER_VR(IND)),                                  &
              LEV_MJ(NUMBER_VR(IND)),DBLE(WR(IND))+EAV,                &
              (DBLE(WR(INDX(1)))-DBLE(WR(IND)))*Rydberg*2,             &
              (DBLE(WR(INDX(1)))-DBLE(WR(IND)))*Convert_meV,           &
              (DBLE(WR(I_FIRST))-DBLE(WR(IND)))*Rydberg*2,             &
              (DBLE(WR(I_FIRST))-DBLE(WR(IND)))*Convert_meV
           END IF
        ENDDO
      ELSE
        J = 1
        IND = INDX(J)
!GG        I_FIRST = 1
        I_FIRST = IND
        WRITE (29,'(I3,1X,I3,1X,A4,A4,I3,4X,F14.7,F12.2,F12.2)')       &
            NUMBER_VR(IND),LEV(NUMBER_VR(IND)),LEV_J(NUMBER_VR(IND)),  &
            LABP_J(NUMBER_VR(IND)),                                    &
            LEV_MJ(NUMBER_VR(IND)),DBLE(WR(IND))+EAV
        DO J = 2,I_MATRICA
           IND = INDX(J)
           IF(LEV_J(NUMBER_VR(IND)) /= LEV_J(NUMBER_VR(INDX(J-1)))     &
           .OR.                                                        &
           LEV(NUMBER_VR(IND)) /= LEV(NUMBER_VR(INDX(J-1)))) THEN
              I_FIRST = IND
              WRITE (29,'(I3,1X,I3,1X,A4,A4,I3,4X,F14.7,2F12.2)')      &
              NUMBER_VR(IND),LEV(NUMBER_VR(IND)),LEV_J(NUMBER_VR(IND)),&
              LABP_J(NUMBER_VR(IND)),                                  &
              LEV_MJ(NUMBER_VR(IND)),DBLE(WR(IND))+EAV,                &
              (DBLE(WR(IND))-DBLE(WR(INDX(1))))*Rydberg*2,             &
              (DBLE(WR(IND))-DBLE(WR(INDX(1))))*Convert_meV
           ELSE
              WRITE (29,'(I3,1X,I3,1X,A4,A4,I3,4X,F14.7,4F12.2)')      &
              NUMBER_VR(IND),LEV(NUMBER_VR(IND)),LEV_J(NUMBER_VR(IND)),&
              LABP_J(NUMBER_VR(IND)),                                  &
              LEV_MJ(NUMBER_VR(IND)),DBLE(WR(IND))+EAV,                &
              (DBLE(WR(IND))-DBLE(WR(INDX(1))))*Rydberg*2,             &
              (DBLE(WR(IND))-DBLE(WR(INDX(1))))*Convert_meV,           &
              (DBLE(WR(IND))-DBLE(WR(I_FIRST)))*Rydberg*2,             &
              (DBLE(WR(IND))-DBLE(WR(I_FIRST)))*Convert_meV
           END IF
        END DO
      END IF
      WRITE (29,*) '------------------------------------------------', &
                   '------------------------------------'
      WRITE (29,'(A)')' '
      WRITE (29,'(A)')'ASF composition '
      WRITE (29,'(A)')' '
      WRITE (31,'(6X,A,I5)')'The space of CF expansion is ',I_MATRICA
      DO IND = 1,I_MATRICA
        WRITE (29,'(I3,A,I3,A,A4,A,I3,5X,A,A4,5X,A,I3)')IND,           &
        '   Number of Energy level',LEV_NO(IND),                       &
        '   Parity',LABP_J(IND),                                       &
        'Pos',LEV(IND),'J=',LEV_J(IND),'2*M_J=',LEV_MJ(IND)
         WRITE (31,'(A)')' '
        WRITE (31,'(I3,A,I3,A,A4,A,I3,5X,A,A4,5X,A,I3)')IND,           &
        '   Number of Energy level',LEV_NO(IND),                       &
        '   Parity',LABP_J(IND),                                       &
        'Pos',LEV(IND),'J=',LEV_J(IND),'2*M_J=',LEV_MJ(IND)
      END DO
      WRITE (31,'(A)')' '
      WRITE (31,'(A)')'End of the list of CF expansion.'
      DO J = 1,I_MATRICA
         IND = INDX(J)
         WRITE (31,'(A)')' '
         WRITE (31,'(A,I3,A,I3,A,A4,A,I3,5X,A,A4,5X,A,I3)')            &
         'Energy level No. before CF ',LEV_NO(NUMBER_VR(IND)),         &
         ' after CF ',NUMBER_VR(IND),                                  &
         ' Parity',LABP_J(NUMBER_VR(IND)),                             &
         'Pos',LEV(NUMBER_VR(IND)),'J=',LEV_J(NUMBER_VR(IND)),         &
         '2*M_J=',LEV_MJ(NUMBER_VR(IND))
         WRITE (31,'(6X,A,F14.7,2X,I5)')                               &
         'Energy = ',DBLE(WR(IND))+EAV,NUMBER_VR(IND)
         WRITE(31,999)(DBLE(VR(JI,IND)),JI=1,I_MATRICA)
         WRITE (31,'(A)')' '
!
         WRITE (29,'(F14.7)')DBLE(WR(IND))+EAV
         WRITE(29,999)(DBLE(VR(JI,IND)),JI=1,I_MATRICA)
         WRITE (29,'(A)')' '
      END DO
!
      CALL WGHTD5gg(iatjpo(Min_state),iaspar(Min_state),Min_state,     &
                                                            Max_state)
!
  101 FORMAT (//' CUTOFF set to ',1PD22.15,'NVEC = ',I10)
  102 FORMAT (/' Crystal field matrix elements :'//,                   &
      ' Lev1  J Par  MJ  Lev2  J Par  MJ',4X,'CF Real (a.u.)',         &
      6X,'CF Imaginary (a.u.)'/)
  103 FORMAT (1X,1I3,1X,2A4,F4.1,2X,I3,1X,2A4,F4.1,1P,1D20.10,         &
                                      1D20.10)
  113 FORMAT (1X,1I3,1X,2A4,F4.1,2X,I3,1X,2A4,F4.1,1P,1D20.10,         &
                                      1D20.10,1D20.10,1D20.10)
  105 FORMAT (' Rydberg constant is ',F14.5)
!GG  999 FORMAT(5X,I3,5(F10.4),/,9(8X,5(F10.4),/))
  999 FORMAT(10(8X,5(F10.4),/))
      RETURN
      END SUBROUTINE CF_Hamil
