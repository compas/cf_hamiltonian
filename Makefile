.SUFFIXES: .f90 .mod

EXE = CF_Halimtonian
BINDIR = ${GRASP}/bin
GRASPLIB = ${GRASP}/lib
BINFILE = $(BINDIR)/$(EXE)
SRCLIBDIR = ${GRASP}/src/lib
MODDIR = ${SRCLIBDIR}/libmod
MODL92 = ${SRCLIBDIR}/lib9290
MODLRANG90 = ${SRCLIBDIR}/librang90
MODLMCP90 = ${SRCLIBDIR}/libmcp90
GRASPLIBS =-l9290 -lrang90 -lmcp90 -lmod


APP_LIBS = -L ${GRASPLIB} ${GRASPLIBS}

#   Define data types
VASTO = ${MODDIR}/vast_kind_param_M.o

APP_OBJ= \
        memory_man_CF.o opt6_C.o CrysPAR_C.o \
        clrx_I.o engouth_I.o gethfd_I.o setdbg_I.o \
        setsum_I.o strsum_I.o getmixblock_I.o index_I.o\
        CF_Hamil_I.o rint_CF_Hamil_I.o wghtd5gg_I.o \
        matelt_CF_Hamil_I.o matelt_CF_Hamil_NEW_I.o\
        skfun_I.o Ions_input_I.o Ions_param_I.o Ions_param_Cart_I.o \
        Ions_param_Sphe_I.o Y_k_DK_I.o Y_k_I.o wigner_3j_I.o plgndr_I.o \
        oneparticlejj_CF_I.o \
\
        clrx.o engouth.o gethfd.o setdbg.o \
        setsum.o strsum.o getmixblock.o index.o\
        CF_Hamil.o CF_Hamiltonian.o rint_CF_Hamil.o wghtd5gg.o \
        matelt_CF_Hamil.o matelt_CF_Hamil_NEW.o\
        skfun.o Ions_input.o Ions_param.o Ions_param_Cart.o \
        Ions_param_Sphe.o Y_k_DK.o Y_k.o wigner_3j.o plgndr.o \
        oneparticlejj_CF.o


$(EXE): $(APP_OBJ)
	$(FC) -o $(BINFILE) $(FC_LD) $(APP_OBJ) $(APP_LIBS) ${LAPACK_LIBS}
.f90.o:
	$(FC) -c $(FC_FLAGS) $< -I $(MODDIR) -I ${MODL92} -I $(MODLRANG90) -I $(MODLMCP90) \
		-I ${GRASP}/src/appl/CF_Hamiltonian -I $(MODDIR) -o $@

.f.o:
	$(FC) -c $(FC_FLAGS) $< -o $@

clean:
	-rm -f *.o core  *.mod
