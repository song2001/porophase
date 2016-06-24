ALL: testp
CFLAGS =
#FFLAGS = -ffpe-trap=invalid,zero,overflow
CPPFLAGS =
FPPFLAGS =
CLEANFILES = testp

include ${PETSC_DIR}/lib/petsc/conf/variables
include ${PETSC_DIR}/lib/petsc/conf/rules

testp: src/poroPetscMu_para_ele.o chkopts
	${FLINKER} -o testp src/poroPetscMu_para_ele.o ${PETSC_LIB}
	${RM} src/poroPetscMu_para_ele.o

