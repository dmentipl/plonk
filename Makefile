OPENMP          = yes
SYSTEM          = gfortran
DOUBLEPRECISION = no
DEBUG           = no

ifeq ($(SYSTEM), gfortran)
    FC        = gfortran
    FFLAGS    = -O3
    DBLFLAGS  = -fdefault-real-8 -fdefault-double-8
    DEBUGFLAG = -Wall -Wextra -pedantic -g -frange-check -fcheck=all \
                -fbacktrace -finit-real=NaN #-ffpe-trap=invalid,zero,overflow
    OMPFLAGS  = -fopenmp
endif

ifdef KERNEL
    ifeq ($(KERNEL), cubic)
        FFLAGS += -DKERNEL=1
    else ifeq ($(KERNEL), quartic)
        FFLAGS += -DKERNEL=2
    else ifeq ($(KERNEL), quintic)
        FFLAGS += -DKERNEL=3
    else ifeq ($(KERNEL), quartic2h)
        FFLAGS += -DKERNEL=4
    else ifeq ($(KERNEL), wendlandc2)
        FFLAGS += -DKERNEL=5
    else ifeq ($(KERNEL), wendlandc4)
        FFLAGS += -DKERNEL=6
    else ifeq ($(KERNEL), wendlandc6)
        FFLAGS += -DKERNEL=7
    else
        FFLAGS += -DKERNEL=99
    endif
endif

ifeq ($(DOUBLEPRECISION), yes)
    FFLAGS += ${DBLFLAGS}
	@echo " DOUBLEPRECISION=yes will probably fail"
endif

ifeq ($(DEBUG),yes)
    FFLAGS += $(DEBUGFLAG)
endif

ifeq ($(OPENMP),yes)
    FFLAGS += $(OMPFLAGS)
endif

splash:
	@echo
	@echo " Compiling SPLASH interpolation routines"
	@echo
	@echo "   Fortran compiler: $(FC)"
	@echo "   Fortran flags:    $(FFLAGS)"
	@echo
	@echo " See Makefile for compile time options"
	@echo
	cd plonk/visualization/splash; \
	$(FC) $(FFLAGS) -c splash.F90; \
	$(FC) $(FFLAGS) -c splash_wrapper.f90; \
	$(FC) -shared -fPIC $(FFLAGS) -o libsplashwrapper.so *.o; \
	python setup.py build_ext --inplace; \
	cd -

clean:
	@rm plonk/visualization/splash/*.o;     \
	rm plonk/visualization/splash/*.mod;    \
	rm plonk/visualization/splash/*.so;     \
	rm -r plonk/visualization/splash/build; \
	rm plonk/visualization/splash/splash_wrapper.c
