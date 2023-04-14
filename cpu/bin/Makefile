DSRC2=../src

# intel ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#F_COMP=ifort
#F_OPTS= -O0 -r8 -check bounds -traceback -check uninit# -prof-gen #-fpe0

# gnu ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
F_COMP=gfortran
F_OPTS= -O0 -fcheck=bounds -fbacktrace -ffree-line-length-200 # -fdefault-real-8

# nvfortran (pgi) ~~~~~~~~~~~~~~~~~~~
# F_COMP=nvfortran
# F_OPTS= -O0 -g -Mbounds -traceback
# useful flags:
# -O[level]        - Optimization level. Where level in 0-4
# -g               - Instructs the compiler to include symbolic debugging information in the object module; 
#                    sets the optimization level to zero unless a -⁠O option is present on the command line.
# -M[no]bounds     - Enables [disables] array bounds checking.
# -[no]-traceback  - Adds debug information for runtime traceback for use with the environment variable NVCOMPILER_TERM. 
# -r4 and -r8      - Interpret SINGLE / DOUBLE PRECISION variables as REAL
# -mp[=gpu]        - Enable OpenMP directives and optionally enable target offload to GPUs.
# -acc[=multicore] - Enable OpenACC directives, defaults to offload to GPUs, use "multicore" to target multicore CPUs.
# -stdpar[=gpu]    - Enable parallelization of DO CONCURRENT to host or GPU
# -Minfo           - The compiler will give feedback messages on what optimization it's applying or not able to apply.
# -help -Minfo     - for the full list of sub-options.


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# TODO - clean or document flags below
#
#F_OPTS= -O0 -ffree-form -fno-range-check -ffree-line-length-none -Wall -fbacktrace -g -fcheck=bounds -fcheck=all -finit-real=nan  -fallow-argument-mismatch
#F_OPTS=  -O0 -fpe2 -ftrapuv -check bounds -traceback -init=snan,arrays  -check uninit  # -ftz  -check uninit -check all -g #-r8 
#F_OPTS= -init=snan,arrays  -check uninit -fpe0 -ftrapuv  -traceback # -ftz  -check uninit -check all -g #-r8 
#F_OPTS=  -O3  -convert big_endian  -fpe2 -ftrapuv -check bounds -traceback -init=snan,arrays  -check uninit  # -ftz  -check uninit -check all -g #-r8 
#F_OPTS= -r4 -check all -g -O0 -check bounds -traceback -ftz -fpe0  -check uninit#-r8 
#F_OPTS= -r4  -g -O0 -ftz -align all -fno-alias -traceback -debug -nolib-inline \
#-fno-inline-functions -assume protect_parens,minus0 -prec-div -prec-sqrt \
#-check bounds -check uninit -fp-stack-check -ftrapuv -warn unused -traceback    \
#-fPIC -fpe0 -fp-model source -heap-arrays 32  -align dcommons
#F_OPTS= -O0 -check all -traceback -ftz -fpe0 #-r8 
#F_OPTS= -O0 -check bounds -traceback #-r8 
#F_OPTS=-O0  -fcheck=bounds -g -fcheck=all -finit-real=nan \
# -fbacktrace -ffpe-trap=invalid,zero,overflow,underflow,denormal #-Wall
#F_OPTS=-O0 -fdefault-real-8 -fcheck=bounds -g -fcheck=all -finit-real=nan \
# -fbacktrace -ffpe-trap=invalid,zero,overflow,underflow,denormal     
#  -finit-real=zero -Wall



LOADER=$(F_COMP)
LOADER_OPTS=-v #-Bstatic -check_bounds -traceback #-Mlarge_arrays
ARCHIVE=ar rs

# Compiler commands.

F_COMMAND = $(F_COMP) -c $(F_OPTS) $(INCLUDES)

# Define archive and executable names.

BASE=wsm
EXE=$(BASE).x
ARC=$(BASE).a

# Define source.
F_SRC=  $(DSRC2)/module_mp_radar.F90 \
        $(DSRC2)/rconstants.F90 \
        $(DSRC2)/module_mp_wsm3.F90 \
        $(DSRC2)/module_mp_wsm5.F90 \
	     $(DSRC2)/module_mp_wsm6.F90 \
	     $(DSRC2)/module_mp_wsm7.F90 \
	     $(DSRC2)/wsm_test.F90 	
	



# Define targets.

all: $(EXE)

$(EXE): $(ARC) 
	$(LOADER) -o $(EXE) $(LOADER_OPTS) $(ARC) $(LIB_IFS)
###	ln -fs $(EXE) $(BASE)
	rm -f r*.o

$(ARC): $(F_SRC)
	$(F_COMMAND) $(?)
	$(ARCHIVE) $(ARC) *.o
#	$(F_COMMAND) $(DSRC)/main.F90
#	rm -f *.o


install:
	ln -fs $(pwd)/$(EXE) ../$(BASE)
#	ln -fs `pwd`/$(EXE) ../test/$(BASE)

clean:
	rm -f $(ARC) $(EXE) $(BASE) *.o  ___* core* *.mod *.a
