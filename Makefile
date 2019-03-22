#
# Makefile to call Splash library Makefile.
# And to clean up Cython build.
#
# Daniel Mentiplay, 2019
#

default:
	@cd splash/fortran; $(MAKE) $(MAKECMDGOALS)

%::
	@cd splash/fortran; $(MAKE) $(MAKECMDGOALS)

clean:
	@echo
	@echo " Removing Cython build files"
	@rm -f splash/splash.c splash/splash.*.so
	@cd splash/fortran; $(MAKE) $(MAKECMDGOALS)
