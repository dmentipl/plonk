"""
Install Plonk.
"""

import io
import pathlib
import re
import sysconfig

from Cython.Distutils import build_ext
from numpy import get_include
from setuptools import Extension, setup

# ---------------------------------------------------------------------------- #
# --- Splash ---

LIBRARY_DIR = str(pathlib.Path(sysconfig.get_path('stdlib')).parent)
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

__version__ = re.search(
    r'__version__\s*=\s*[\'"]([^\'"]*)[\'"]',  # It excludes inline comment too
    io.open('plonk/__init__.py', encoding='utf_8_sig').read(),
).group(1)

cmdclass = dict()
cmdclass.update({'build_ext': build_ext})
include_dirs = [get_include()]

install_requires = ['h5py', 'matplotlib', 'numpy', 'pandas', 'sympy']

description = 'Smoothed particle hydrodynamics analysis and visualization with Python.'

setup(
    name='plonk',
    version=__version__,
    author='Daniel Mentiplay',
    author_email='d.mentiplay@gmail.com',
    url='https://github.com/dmentipl/plonk',
    description=description,
    packages=['plonk', 'plonk.analysis', 'plonk.core', 'plonk.visualization'],
    license='MIT',
    install_requires=install_requires,
    cmdclass=cmdclass,
    include_dirs=include_dirs,
    ext_modules=ext_modules,
)
