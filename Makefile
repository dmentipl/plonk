OPENMP=yes
SYSTEM=gfortran
DOUBLEPRECISION=no
DEBUG=no

ifeq ($(SYSTEM), gfortran)
   FC= gfortran
   FFLAGS= -O3
   DBLFLAGS= -fdefault-real-8 -fdefault-double-8
   DEBUGFLAG= -Wall -Wextra -pedantic -g -frange-check -fcheck=all -fbacktrace \
              -finit-real=NaN #-ffpe-trap=invalid,zero,overflow
   OMPFLAGS= -fopenmp
endif

ifeq ($(DOUBLEPRECISION), yes)
    FFLAGS += ${DBLFLAGS}
endif

ifeq ($(DEBUG),yes)
    FFLAGS += $(DEBUGFLAG)
endif

ifeq ($(OPENMP),yes)
    FFLAGS += $(OMPFLAGS)
endif

splash:
	@echo
	@echo " Compiling SPLASH interpolation routines with f2py"
	@echo
	@echo "   Fortran compiler: $(FC)"
	@echo "   Fortran flags:    $(FFLAGS)"
	@echo
	@echo " See Makefile for compile time options"
	@echo
	cd plonk/visualization; \
	f2py -m _splash -c --fcompiler=$(FC) --opt='$(FFLAGS)' splash.F90; \
	cd -

clean:
	@echo
	@echo " Removing compiled Splash library"
	@echo
	@if [ -f plonk/visualization/_splash.*.so ]; then \
		rm -r plonk/visualization/_splash.*.so; fi
	@if [ -d plonk/visualization/_splash.*.so.dSYM ]; then \
		rm -r plonk/visualization/_splash.*.so.dSYM; fi
	@if [ -d plonk/visualization/__pycache__ ]; then \
		rm -r plonk/visualization/__pycache__; fi
