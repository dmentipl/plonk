"""
Install Plonk.
"""

import pathlib
import sysconfig

from setuptools import Extension, setup

from Cython.Distutils import build_ext
from numpy import get_include

# ---------------------------------------------------------------------------- #
# --- Splash ---

LIBRARY_DIR = pathlib.Path(sysconfig.get_path('stdlib')).parent
SOURCES = ['splash/splash.pyx']
LIBRARIES = ['splash', 'gfortran']

ext_modules = [
    Extension(
        'splash.splash',
        sources=SOURCES,
        libraries=LIBRARIES,
        library_dirs=[LIBRARY_DIR],
        runtime_library_dirs=[LIBRARY_DIR],
    )
]

# ---------------------------------------------------------------------------- #
# --- Plonk ---

cmdclass = dict()
cmdclass.update({'build_ext': build_ext})
include_dirs = [get_include()]

install_requires = ['cython', 'h5py', 'matplotlib', 'numpy', 'pandas']

description = (
    'Smoothed particle hydrodynamics analysis and visualization with Python.'
)

setup(
    name='Plonk',
    version='0.1',
    author='Daniel Mentiplay',
    author_email='d.mentiplay@gmail.com',
    url='https://github.com/dmentipl/plonk',
    description=description,
    packages=['plonk', 'plonk.analysis', 'plonk.tests', 'plonk.visualization'],
    license='MIT',
    install_requires=install_requires,
    cmdclass=cmdclass,
    include_dirs=include_dirs,
    ext_modules=ext_modules,
)
