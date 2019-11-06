"""Plonk setup.py."""

import io
import re

from setuptools import setup

# ---------------------------------------------------------------------------- #
# --- Plonk ---

__version__ = re.search(
    r'__version__\s*=\s*[\'"]([^\'"]*)[\'"]',  # It excludes inline comment too
    io.open('plonk/__init__.py', encoding='utf_8_sig').read(),
).group(1)

install_requires = [
    'KDEpy',
    'astropy',
    'h5py',
    'matplotlib',
    'numpy',
    'pandas',
    'pint',
    'skimage>=0.16',
    'scipy',
]
packages = [
    'plonk',
    'plonk.analysis',
    'plonk.simulation',
    'plonk.snap',
    'plonk.snap.readers',
    'plonk.utils',
    'plonk.visualization',
]

description = 'Smoothed particle hydrodynamics analysis and visualization with Python.'

setup(
    name='plonk',
    version=__version__,
    author='Daniel Mentiplay',
    author_email='d.mentiplay@gmail.com',
    url='https://github.com/dmentipl/plonk',
    description=description,
    packages=packages,
    license='MIT',
    install_requires=install_requires,
)
