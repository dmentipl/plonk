#
# Compile Splash Fortran source to produce libsplash.so.
#
# The main targets are:
#   - build:   to build the shared object library
#   - install: to build and install the library to INSTALL_DIR
#   - conda:   to build the library when building the Plonk Conda package
#
# Daniel Mentiplay, 2019
#

SHELL := /bin/bash

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
DEBUG  = no

FC      = gfortran
FCFLAGS = -O3 -fPIC
FLFLAGS = -shared

FC_BIN = $(shell command -v $(FC))
FC_VERSION = $(shell $(FC) --version | head -1 | awk '{print $$NF}')

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

help:
	@echo "Usage:"
	@echo "  make <target>"
	@echo
	@echo "  build            Build libsplash"
	@echo "  clean            Cleanup build files"
	@echo "  conda            Build libsplash for Plonk Conda package"
	@echo "  help             Display this help"
	@echo "  install          Build libsplash and install to INSTALL_DIR"

.PHONY: print
print:
	@echo
	@echo " Compiling Splash modules"
	@echo " Fortran compiler binary: $(FC_BIN)"
	@echo " Fortran compiler version: $(FC_VERSION)"
	@echo

.PHONY: build
build: print $(OBJECTS)
	@echo
	@echo " Linking Splash into a shared object library"
	@echo
	$(FC) $(FCFLAGS) $(FLFLAGS) $(OBJECTS) -o $(LIBSPLASH)

.PHONY: install
install: build
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
conda: build
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
