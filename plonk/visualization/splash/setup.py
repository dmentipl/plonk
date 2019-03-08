'''
setup.py
'''

import os.path

from setuptools import setup, Extension
from Cython.Distutils import build_ext
import numpy as np

here = os.path.dirname(os.path.abspath(__file__))

ext_modules = [Extension("_splash",
                         ["splash_wrapper.pyx"],
                         libraries=['splashwrapper', 'gfortran'],
                         library_dirs=[here],
                         runtime_library_dirs=[here],
                        )]

setup(
    name = 'Splash',
    cmdclass = {'build_ext': build_ext},
    include_dirs = [np.get_include()],
    ext_modules = ext_modules
)
