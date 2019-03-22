'''
Install Plonk.
'''

import os
import sysconfig

from setuptools import setup, Extension
from Cython.Distutils import build_ext
from numpy import get_include

# ---------------------------------------------------------------------------- #
# --- Splash ---

LIBRARY_DIR = os.path.dirname(sysconfig.get_path('stdlib'))
SOURCES     = ['splash/splash.pyx']
LIBRARIES   = ['splash', 'gfortran']

ext_modules = [Extension('splash.splash',
                         sources=SOURCES,
                         libraries=LIBRARIES,
                         library_dirs=[LIBRARY_DIR],
                         runtime_library_dirs=[LIBRARY_DIR],
                        )]

# ---------------------------------------------------------------------------- #
# --- Plonk ---

cmdclass = dict()
cmdclass.update({'build_ext': build_ext})
include_dirs = [get_include()]

install_requires = \
'''
cython
h5py
matplotlib
numpy
pandas
'''

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