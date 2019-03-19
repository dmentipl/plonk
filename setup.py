'''
Install Plonk.
'''

import os
import platform
import subprocess

from setuptools import setup, Extension
from Cython.Distutils import build_ext
from numpy import get_include

# ---------------------------------------------------------------------------- #
# --- Splash ---

FORCE_RECOMPILE_SPLASH = False

LIBSPLASH      = 'libsplash.so'
LIBSPLASH_DIR  = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                              'splash')
LIBSPLASH_PATH = os.path.join(LIBSPLASH_DIR, LIBSPLASH)

def compile_splash():
    '''
    Compile Splash Fortran source.
    '''

    OPENMP = True
    DOUBLE = False
    DEBUG  = False
    KERNEL = 'cubic'

    FSOURCES  = ['splash.F90', 'libsplash.f90'] # Fortran sources in order

    FC      = 'gfortran'
    FCFLAGS = ['-O3', '-fPIC']
    FLFLAGS = ['-shared']

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
        source_file_path = os.path.join(LIBSPLASH_DIR, file)
        object_file_path = os.path.join(LIBSPLASH_DIR, file[:-4]+'.o')
        compile_command = ' '.join([FC] + FCFLAGS +
                                   ['-J' + LIBSPLASH_DIR] +
                                   ['-c', source_file_path] +
                                   ['-o', object_file_path])
        print(compile_command)
        subprocess.run(compile_command, shell=True)

    # Link to make library
    deps = os.path.join(LIBSPLASH_DIR, '*.o')
    compile_command = ' '.join([FC] + FCFLAGS + FLFLAGS + [deps] +
                               ['-o', LIBSPLASH_PATH])
    print(compile_command)
    subprocess.run(compile_command, shell=True)

#--- Compile Splash source

if not os.path.isfile(LIBSPLASH_PATH) or FORCE_RECOMPILE_SPLASH:
    compile_splash()

# Fix for macOS systems
extra_link_args = list()
if platform.system() == 'Darwin':
    extra_link_args.append('-Wl,-rpath,'+LIBSPLASH_DIR)
    if not os.path.isfile(LIBSPLASH_PATH):
        subprocess.run('install_name_tool -id @rpath/'+LIBSPLASH+' ' + \
                       LIBSPLASH_PATH, shell=True)

#--- Splash extension module

ext_modules = [Extension('splash',
                         [os.path.join(LIBSPLASH_DIR, 'libsplash.pyx')],
                         libraries=['splash', 'gfortran'],
                         library_dirs=[LIBSPLASH_DIR],
                         runtime_library_dirs=[LIBSPLASH_DIR],
                         extra_link_args=extra_link_args,
                        )]

# ---------------------------------------------------------------------------- #
# --- Plonk ---

cmdclass = dict()
cmdclass.update({'build_ext': build_ext})
include_dirs = [get_include()]

with open(os.path.join(os.path.dirname(__file__), 'requirements.txt')) as fp:
    install_requires = fp.read()

setup(
    name='Plonk',
    version='0.1',
    author='Daniel Mentiplay',
    author_email='d.mentiplay@gmail.com',
    url='https://github.com/dmentipl/plonk',
    description='Phantom analysis and visualization but with Python.',
    packages=['plonk', 'plonk.analysis', 'plonk.tests', 'plonk.visualization'],
    license='MIT',
    install_requires=install_requires,
    cmdclass=cmdclass,
    include_dirs=include_dirs,
    ext_modules=ext_modules,
)
