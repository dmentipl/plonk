splash:
	@echo " Compiling SPLASH interpolation routines"
	cd plonk/visualization; \
	f2py -m _splash -c splash.F90; \
	cd -
