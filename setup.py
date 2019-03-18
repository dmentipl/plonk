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

libsplash_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                             'splash')

def compile_splash():
    '''
    Compile Splash source.
    '''

    OPENMP = True
    DOUBLE = False
    DEBUG  = False
    KERNEL = 'cubic'

    FSOURCES  = ['splash.F90', 'libsplash.f90'] # Fortran sources in order
    LIBSPLASH = 'libsplash.so'                  # Splash library name

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
        source_file_path = os.path.join(libsplash_dir, file)
        object_file_path = os.path.join(libsplash_dir, file[:-4]+'.o')
        compile_command = ' '.join([FC] + FCFLAGS +
                                   ['-J' + libsplash_dir] +
                                   ['-c', source_file_path] +
                                   ['-o', object_file_path])
        print(compile_command)
        subprocess.run(compile_command, shell=True)

    # Link to make library
    file_path = os.path.join(libsplash_dir, LIBSPLASH)
    deps = os.path.join(libsplash_dir, '*.o')
    compile_command = ' '.join([FC] + FCFLAGS + FLFLAGS + [deps] + ['-o', file_path])
    print(compile_command)
    subprocess.run(compile_command, shell=True)

    # Fix for macOS systems
    extra_link_args = list()
    if platform.system() == 'Darwin':
        subprocess.run('install_name_tool -id @rpath/'+LIBSPLASH+' ' + \
                       libsplash_dir+'/'+LIBSPLASH, shell=True)
        extra_link_args.append('-Wl,-rpath,'+libsplash_dir)

    return extra_link_args

#--- Compile Splash source

extra_link_args = compile_splash()

#--- Splash extension module

ext_modules = [Extension('splash',
                         [os.path.join(libsplash_dir, 'libsplash.pyx')],
                         libraries=['splash', 'gfortran'],
                         library_dirs=[libsplash_dir],
                         runtime_library_dirs=[libsplash_dir],
                         extra_link_args=extra_link_args,
                        )]

# ---------------------------------------------------------------------------- #
# --- Plonk ---

with open(os.path.join(os.path.dirname(__file__), 'requirements.txt')) as fp:
    install_requires = fp.read()

setup(
    name='Plonk',
    version='0.1',
    author='Daniel Mentiplay',
    author_email='d.mentiplay@gmail.com',
    packages=['plonk', 'plonk.analysis', 'plonk.visualization'],
    install_requires=install_requires,
    url='https://github.com/dmentipl/plonk',
    license='MIT',
    description='Phantom analysis and visualization but with Python.',
    cmdclass = {'build_ext': build_ext},
    include_dirs = [get_include()],
    ext_modules = ext_modules
)
