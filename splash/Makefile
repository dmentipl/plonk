#
# Compile Splash Fortran source to produce libsplash.so.
# Daniel Mentiplay, 2018
#

OPENMP = yes
DOUBLE = no
DEBUG  = no

FC      = gfortran
FCFLAGS = -O3 -fPIC
FLFLAGS = -shared

ifeq ($(DOUBLE), yes)
    FCFLAGS += -fdefault-real-8 -fdefault-double-8
endif

ifeq ($(DEBUG), yes)
    FCFLAGS += -Wall -Wextra -pedantic -g -frange-check -fcheck=all \
			   -fbacktrace -finit-real=NaN
endif

ifeq ($(OPENMP), yes)
    FCFLAGS += -fopenmp
endif

%.o : %.f90
	$(FC) -c $(FCFLAGS) $< -o $@

%.o : %.F90
	$(FC) -c $(FCFLAGS) $< -o $@

FSOURCES = splash.F90 libsplash.f90

OBJECTS1 = $(FSOURCES:.f90=.o)
OBJECTS  = $(OBJECTS1:.F90=.o)

LIBRARY_DIR = $(CONDA_PREFIX)/lib

libsplash: $(OBJECTS)
	$(FC) $(FCFLAGS) $(FLFLAGS) $(OBJECTS) -o $@.so
	mv $@.so $(LIBRARY_DIR)
	install_name_tool -id @rpath/$@.so $(LIBRARY_DIR)/$@.so

.PHONY: clean
clean:
	rm -f *.o *.mod

.PHONY: cleanall
cleanall:
	rm -f *.o *.mod *.so libsplash.c
