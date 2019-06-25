#
# Compile Splash Fortran source to produce libsplash.so.
#
# The main targets are:
#   - libsplash: to build the shared object library
#   - install: to build and install the library to INSTALL_DIR
#   - conda: to build the library when building the Plonk Conda package
#
# Daniel Mentiplay, 2019
#

# Name of the shared object library
LIBSPLASH = libsplash.so

# Install directory for Anaconda
# Specify INSTALL_DIR if you want to install elsewhere
ifeq (X$(CONDA_PREFIX), X)
CONDA_PREFIX = $(HOME)/anaconda
endif
INSTALL_DIR = $(CONDA_PREFIX)/lib

UNAME = ${shell uname}

OPENMP = yes
DOUBLE = no
DEBUG  = no

FC      = gfortran
FCFLAGS = -O3 -fPIC
FLFLAGS = -shared

# DOUBLE will break Cython
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

default: libsplash

.PHONY: print
print:
	@echo
	@echo " Compiling Splash modules"
	@echo

.PHONY: libsplash
libsplash: print $(OBJECTS)
	@echo
	@echo " Linking Splash into a shared object library"
	@echo
	$(FC) $(FCFLAGS) $(FLFLAGS) $(OBJECTS) -o $(LIBSPLASH)

.PHONY: install
install: libsplash
	@echo
	@echo " Copying libsplash.so to INSTALL_DIR=$(INSTALL_DIR)"
	@echo
	cp -f $(LIBSPLASH) $(INSTALL_DIR)
ifeq ($(UNAME), Darwin)
	@echo
	@echo " Fixing macOS library path issues with install_name_tool"
	@echo
	install_name_tool -id @rpath/$(LIBSPLASH) $(INSTALL_DIR)/$(LIBSPLASH)
endif

.PHONY: conda
conda: libsplash
	@echo
	@echo " Copying libsplash.so to PREFIX/lib=$(PREFIX)/lib"
	@echo " Only use this Makefile target for building the Conda package"
	@echo
	cp -f $(LIBSPLASH) $(PREFIX)/lib

.PHONY: clean
clean:
	@echo
	@echo " Removing Fortran build files"
	@echo
	@rm -f *.o *.mod *.so *.dSYM
