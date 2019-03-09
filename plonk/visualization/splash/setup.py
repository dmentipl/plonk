'''
setup.py

Daniel Mentiplay, 2019.
'''

import os.path
import platform
import subprocess

from setuptools import setup, Extension
from Cython.Distutils import build_ext
from numpy import get_include

OPENMP = True
DOUBLE = False
DEBUG  = False
KERNEL = 'cubic'

FSOURCES  = ['splash.F90', 'splash_wrapper.f90'] # Fortran sources in order
LIBSPLASH = 'libsplashwrapper.so'                # Splash library name

local_dir = os.path.dirname(os.path.abspath(__file__))

FC      = 'gfortran'
FCFLAGS = ['-O3', '-fPIC']
FLFLAGS = ['-shared', '-o '+LIBSPLASH, '*.o']

KERNELS = ['cubic', 'quartic', 'quintic', 'quartic2h', 'wendlandc2',
           'wendlandc4', 'wendlandc6']

if KERNEL:
    kernel_set = False
    for idx, K in enumerate(KERNELS):
        if K == KERNEL:
            kernel_set = True
            FCFLAGS += [f'-DKERNEL={idx+1}']
            break
    if not kernel_set:
        FCFLAGS += ['-DKERNEL=99']
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
print(' See setup.py for compile time options')
print()

# Compile
for file in FSOURCES:
    subprocess.run(' '.join([FC] + FCFLAGS + ['-c', file]), shell=True)

# Link to make library
subprocess.run(' '.join([FC] + FCFLAGS + FLFLAGS), shell=True)

# Fix for macOS systems
extra_link_args = list()
if platform.system() == 'Darwin':
    subprocess.run('install_name_tool -id @rpath/'+LIBSPLASH+' ' + \
                   local_dir+'/'+LIBSPLASH, shell=True)
    extra_link_args.append('-Wl,-rpath,'+local_dir)

# Make Splash library available to Python as _splash
ext_modules = [Extension('_splash',
                         ['splash_wrapper.pyx'],
                         libraries=['splashwrapper', 'gfortran'],
                         library_dirs=[local_dir],
                         runtime_library_dirs=[local_dir],
                         extra_link_args=extra_link_args,
                        )]

setup(
    name = 'Splash',
    cmdclass = {'build_ext': build_ext},
    include_dirs = [get_include()],
    ext_modules = ext_modules
)
