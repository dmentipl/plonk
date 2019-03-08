Plonk
=====

Phantom analysis and visualization but with Python.

Install
-------

To use visualization tools you first must compile Splash Fortran interpolation routines. To do this, type `make` in the main repository directory.

You must also set the `LD_LIBRARY_PATH` (for Linux) or `DYLD_LIBRARY_PATH` (for
macOS). For example:

```
export DYLD_LIBRARY_PATH=DYLD_LIBRARY_PATH:~/plonk/plonk/visualization/splash
```
