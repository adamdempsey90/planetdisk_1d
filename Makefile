EXECUTABLE=pdisk1d
SOURCES=advance_system.c allocate.c boundary.c dTr_ex.c disk.c hdf5.c initialize.c integration.c main.c matrix_ops.c mdot.c move_planet.c read_params.c set_grid.c set_matrix.c set_uw.c steady_state.c total_torque.c trisolve.c interpolation.c
HEADER=defines.h pdisk.h

LAPACKLIB=-llapack -lblas
OPENMPLIB=-lgomp
MATHLIB=-lm
GSLLIB=-lgsl -lgslcblas
HDF5LIB=-lhdf5
WARNFLAGS=-Wall -Wextra
OPENMPFLAG=-fopenmp
HDF5FLAG=-DH5_USE_16_API




INCLIB=
LDLIB=$(GSLLIB)


BIN=bin/
SRC=src/
#IN=inputs/
#PY=src/pyutils/


UNAME := $(shell echo $(USER))

ifeq ($(UNAME),apollo)
CC=gcc
endif

ifeq ($(UNAME),jupiter)
CC=gcc-6
endif
ifeq ($(UNAME),adam)
CC=gcc
endif
ifeq ($(UNAME),zeus)
CC=gcc
OPENMPLIB=
OPENMPFLAG=
INCLIB=-I/usr/local/include/
LDLIB=-L/usr/local/lib/
#CC=gcc-4.9
#INCLIB=-I/usr/local/include/
#LDLIB=-L/usr/local/lib/
endif


ifeq ($(UNAME),helios)
CC=gcc-4.9
INCLIB=-I/usr/local/include/
LDLIB=-L/usr/local/lib/
endif

ifeq ($(UNAME),amd616)
CC=gcc
LDLIB=-L/software/lapack/3.4.0/lib -L/software/gsl/1.16-gcc4.8.3/lib/ -L/software/hdf5/1.8.12-serial/lib/
INCLIB=-I/software/hdf5/1.8.12-serial/include/
endif

LDFLAGS= $(MATHLIB) $(HDF5LIB) $(OPENMPLIB)
CFLAGS=-c -O3 -g $(WARNFLAGS) $(HDF5FLAG) $(OPENMPFLAG)

#!!!!!DO NOT EDIT ANYTHING UNDER THIS LINE!!!!!
OBJECTS=$(SOURCES:.c=.o)
CSOURCES=$(addprefix $(SRC),$(SOURCES))
COBJECTS=$(addprefix $(BIN),$(OBJECTS))
CHEADER=$(addprefix $(SRC),$(HEADER))




all: $(CSOURCES) $(EXECUTABLE)

$(EXECUTABLE): $(COBJECTS)
	$(CC)  $(COBJECTS) $(LDLIB) $(LDFLAGS) -o $@ 

$(BIN)%.o: $(SRC)%.c $(CHEADER)
	$(CC) $(INCLIB) $(CFLAGS) $< -o $@ 


clean:
	rm $(COBJECTS) $(EXECUTABLE)
