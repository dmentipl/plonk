'''
setup.py
'''

import os.path
import subprocess

from setuptools import setup, Extension
from Cython.Distutils import build_ext
from numpy import get_include

OPENMP = True
DOUBLE = False
DEBUG  = False
KERNEL = 'cubic'

FSOURCES = ['splash.F90', 'splash_wrapper.f90']

FC      = 'gfortran'
FCFLAGS = ['-O3']
FCLIBS  = ['-shared', '-fPIC', '-o libsplashwrapper.so', '*.o']

KERNELS = ['cubic', 'quartic', 'quintic', 'quartic2h', 'wendlandc2',
           'wendlandc4', 'wendlandc6']

if KERNEL:
    for idx, K in enumerate(KERNELS):
        if K == KERNEL:
            FCFLAGS += [f'-DKERNEL={idx+1}']
else:
    FCFLAGS += ['-DKERNEL=99']

if DOUBLE:
    FCFLAGS += ['-fdefault-real-8', '-fdefault-double-8']
    raise ValueError('DOUBLE=yes will probably fail')

if DEBUG:
    FCFLAGS += ['-Wall', '-Wextra', '-pedantic', '-g', '-frange-check',
                '-fcheck=all', '-fbacktrace', '-finit-real=NaN']

if OPENMP:
    FCFLAGS += ['-fopenmp']

print()
print(' Compiling SPLASH interpolation routines')
print()
print(f'   Fortran compiler: {FC}')
print(f'   Fortran flags:    {" ".join(FCFLAGS)}')
print()
print(' See Makefile for compile time options')
print()

for file in FSOURCES:
    subprocess.run(' '.join([FC] + FCFLAGS + ['-c', file]), shell=True)

subprocess.run(' '.join([FC] + FCFLAGS + FCLIBS), shell=True)

local_dir = os.path.dirname(os.path.abspath(__file__))

ext_modules = [Extension('_splash',
                         ['splash_wrapper.pyx'],
                         libraries=['splashwrapper', 'gfortran'],
                         library_dirs=[local_dir],
                         runtime_library_dirs=[local_dir],
                        )]

setup(
    name = 'Splash',
    cmdclass = {'build_ext': build_ext},
    include_dirs = [get_include()],
    ext_modules = ext_modules
)
