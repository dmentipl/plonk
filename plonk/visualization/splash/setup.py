'''
setup.py
'''

import os.path
import subprocess

from setuptools import setup, Extension
from Cython.Distutils import build_ext
import numpy as np

OPENMP          = True
SYSTEM          = 'gfortran'
DOUBLEPRECISION = False
DEBUG           = False
KERNEL          = None

if SYSTEM == 'gfortran':
    FC        = 'gfortran'
    FFLAGS    = '-O3'
    DBLFLAGS  = '-fdefault-real-8 -fdefault-double-8'
    DEBUGFLAG = '-Wall -Wextra -pedantic -g -frange-check -fcheck=all \
                -fbacktrace -finit-real=NaN #-ffpe-trap=invalid,zero,overflow'
    OMPFLAGS  = '-fopenmp'
else:
    raise ValueError('Compiler not available; try gfortran')

if KERNEL:
    if KERNEL == 'cubic':
        FFLAGS += ' -DKERNEL=1'
    elif KERNEL == 'quartic':
        FFLAGS += ' -DKERNEL=2'
    elif KERNEL == 'quintic':
        FFLAGS += ' -DKERNEL=3'
    elif KERNEL == 'quartic2h':
        FFLAGS += ' -DKERNEL=4'
    elif KERNEL == 'wendlandc2':
        FFLAGS += ' -DKERNEL=5'
    elif KERNEL == 'wendlandc4':
        FFLAGS += ' -DKERNEL=6'
    elif KERNEL == 'wendlandc6':
        FFLAGS += ' -DKERNEL=7'
    else:
        FFLAGS += ' -DKERNEL=99'
else:
    FFLAGS += ' -DKERNEL=99'

if DOUBLEPRECISION:
    FFLAGS += ' ' + DBLFLAGS
    print("DOUBLEPRECISION=yes will probably fail")

if DEBUG:
    FFLAGS += ' ' + DEBUGFLAG

if OPENMP:
    FFLAGS += ' ' + OMPFLAGS

print()
print(" Compiling SPLASH interpolation routines")
print()
print(f"   Fortran compiler: {FC}")
print(f"   Fortran flags:    {FFLAGS}")
print()
print(" See Makefile for compile time options")
print()

local_dir = os.path.dirname(os.path.abspath(__file__))

ffiles = ['splash.F90', 'splash_wrapper.f90']
fc_command = FC + ' ' + FFLAGS + ' -c '
fl_command = FC + ' -shared -fPIC ' + FFLAGS + ' -o libsplashwrapper.so *.o'

for ffile in ffiles:
    subprocess.run(fc_command + ffile, shell=True)

subprocess.run(fl_command, shell=True)

ext_modules = [Extension("_splash",
                         ["splash_wrapper.pyx"],
                         libraries=['splashwrapper', 'gfortran'],
                         library_dirs=[local_dir],
                         runtime_library_dirs=[local_dir],
                        )]

setup(
    name = 'Splash',
    cmdclass = {'build_ext': build_ext},
    include_dirs = [np.get_include()],
    ext_modules = ext_modules
)
