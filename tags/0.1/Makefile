all: ida 
CFLAGS = -Wall -Wextra -pedantic -std=c99 -O3
FFLAGS =
CPPFLAGS =
FPPFLAGS =

include ${PETSC_DIR}/conf/variables
include ${PETSC_DIR}/conf/rules

OBJ=check.o compare_true.o exact_solve.o form_A_matrix.o hybrid_ida.o \
		initial_guess.o input_params.o io.o main.o regular_ida.o


ida: $(OBJ) chkopts
	${CLINKER} -o ida $(OBJ) ${PETSC_LIB}
	${RM} $(OBJ)
	mkdir -p ../bin
	mv ida ../bin
