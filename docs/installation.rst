================
Installing Plonk
================

To install Plonk::

   python setup.py install

------------
Requirements
------------

Plonk has Python requirements listed in ``requirements.txt``. The recommended way to satisfy these requirements is to use `Anaconda <https://anaconda.org/>`_.

In addition, Plonk requires a recent (2003+) version of the GCC Fortran compiler `gfortran <https://gcc.gnu.org/wiki/GFortran>`_ to compile the Splash subroutines.

^^^^^^
Splash
^^^^^^

This will compile the Splash Fortran subroutines and link them into a library called ``libsplash.so``. It then uses Cython to compile them into a Python module.

For Fortran-Python compatibility we use Cython and Fortran 2003+ ``iso_c_binding`` module for C interoperability.
