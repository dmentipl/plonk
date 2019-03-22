#
# Compile Splash Fortran source to produce libsplash.so.
# Daniel Mentiplay, 2018
#

UNAME=${shell uname}

OPENMP = yes
DOUBLE = no
DEBUG  = no

FC      = gfortran
FCFLAGS = -O3
FLFLAGS = -fPIC -shared

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

FSOURCES = sort.f90 asciiutils.f90 kernels.f90 timing.f90 globaldata.f90 \
           labels.f90 interpolation.f90 interpolate3D_projection.F90 \
           interpolate3D_opacity.f90 interpolate3D_xsec.f90 libsplash.f90

OBJECTS1 = $(FSOURCES:.f90=.o)
OBJECTS  = $(OBJECTS1:.F90=.o)

LIBRARY_DIR = $(CONDA_PREFIX)/lib

default: libsplash

libsplash: $(OBJECTS)
	$(FC) $(FCFLAGS) $(FLFLAGS) $(OBJECTS) -o $@.so
	mv $@.so $(LIBRARY_DIR)
ifeq ($(UNAME), Darwin)
	install_name_tool -id @rpath/$@.so $(LIBRARY_DIR)/$@.so
endif

.PHONY: clean
clean:
	rm -f *.o *.mod

.PHONY: cleanall
cleanall:
	rm -f *.o *.mod *.so libsplash.c
