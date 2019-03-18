Plonk
=====

Phantom analysis and visualization but with Python. For usage see `examples` folder.

Install
-------

To install Plonk:

```
python setup.py install
```

This will compile the Splash Fortran subroutines and link them into a library
called `libsplash.so`. It then uses Cython to compile them into a Python module.

For Fortran-Python compatibility we use Cython and Fortran 2003+ `iso_c_binding` module for C interoperability.

Requirements
------------

Plonk has Python requirements listed in `requirements.txt`.

The easiest way to satisfy these requirements is to use Anaconda (https://anaconda.org/).

In addition, Plonk requires a recent (2003+) gfortran (https://gcc.gnu.org/wiki/GFortran) to compile the Splash subroutines.
