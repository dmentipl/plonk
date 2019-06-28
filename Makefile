#
# Makefile to call Splash library Makefile.
# And to clean up Cython build.
#
# Daniel Mentiplay, 2019
#

.PHONY: build install conda clean

build:
	@cd splash/fortran; $(MAKE) $(MAKECMDGOALS)

install:
	@cd splash/fortran; $(MAKE) $(MAKECMDGOALS)

conda:
	@cd splash/fortran; $(MAKE) $(MAKECMDGOALS)

clean:
	@echo
	@echo " Removing Cython build files"
	@rm -f splash/splash.c splash/splash.*.so
	@cd splash/fortran; $(MAKE) $(MAKECMDGOALS)
