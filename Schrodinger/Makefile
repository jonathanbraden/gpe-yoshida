# Default to GNU compiler without OpenMP, no warnings and no profiling
FC = gfortran
OMP = T
PROF = F
WARN = T
BENCH = F
BENCH_FINE = F

# Default FFTW Lib location
FFTW_DIR = /usr/local/
FFTW_INC = -I$(FFTW_DIR)include/
FFTW_LIB = -L$(FFTW_DIR)lib/

# Intel Compiler
ifeq ($(FC),ifort)
	FFLAGS= -r8 -fpp
ifeq ($(OMP),T)
	FFLAGS+= -openmp
endif	
	FOPT= -O3 -ipo -xHost #-fast -no-prec-div -parallel
ifeq ($(WARN),T)
	FFLAGS +=
endif
ifeq ($(PROF),T)
	FFLAGS +=
endif

# GNU Compiler
else ifeq ($(FC),gfortran)
	FFLAGS=-fdefault-real-8 -fdefault-double-8 -cpp -ffree-line-length-none 
ifeq ($(OMP),T)
	FFLAGS += -fopenmp
endif
	FOPT= -O3 #-malign-double -fopt-info-vec-missed=missed.all #-ftree-vectorizer-verbose=5  #-fblas  (check this blas one to replace matmul calls)
ifeq ($(WARN),T)
	FFLAGS+= -Wall -fcheck=array-temps,bounds -pedantic
endif
ifeq ($(PROF),T)
	FFLAGS+= -g -pg
endif
endif

ifeq ($(OMP),T)
	FFTW = -lfftw3_omp -lfftw3
else
	FFTW = -lfftw3
endif

FLIBS = $(FFTW) 
ifeq ($(BLAS),T)
	FLIBS+=-lblas
endif
FLIBS+= -lm

OBJS = constants.o utils.o fftw_mod_wtype.o fast_cheby_2d.o lattice.o equations.o yoshida.o

gp: %: $(OBJS) evolve.o
	$(FC) $(FFLAGS) $(FOPT) $(FFTW_INC) $(MACROS) -o schrodinger evolve.o $(OBJS) $(FFTW_LIB) $(FLIBS) $(THREAD_LIB)


%.o: %.f90
	$(FC) $(FFLAGS) $(FOPT) $(FFTW_INC) $(MACROS) -c $< -o $@ $(FFTW_LIB) $(FLIBS) $(THREAD_LIB)

.PHONY : clean

clean:
	rm -f *.o
	rm -f *.mod
	rm -f gp-2d solve-bg
