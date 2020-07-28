"""Plonk setup.py."""

import io
import pathlib
import re

from setuptools import find_packages, setup

__version__ = re.search(
    r'__version__\s*=\s*[\'"]([^\'"]*)[\'"]',
    io.open('plonk/__init__.py', encoding='utf_8_sig').read(),
).group(1)

install_requires = [
    'h5py',
    'matplotlib',
    'numba',
    'numpy',
    'pandas',
    'pint>=0.10.1',
    'scipy',
]

packages = find_packages(exclude=['tests'])
description = 'Smoothed particle hydrodynamics analysis and visualization with Python.'
long_description = (pathlib.Path(__file__).parent / 'README.md').read_text()

setup(
    name='plonk',
    version=__version__,
    author='Daniel Mentiplay',
    author_email='d.mentiplay@gmail.com',
    url='https://github.com/dmentipl/plonk',
    description=description,
    long_description=long_description,
    long_description_content_type='text/markdown',
    packages=packages,
    license='MIT',
    install_requires=install_requires,
    classifiers=[
        "Development Status :: 5 - Production/Stable",
        "License :: OSI Approved :: MIT License",
        "Programming Language :: Python :: 3.6",
        "Programming Language :: Python :: 3.7",
        "Programming Language :: Python :: 3.8",
        "Operating System :: OS Independent",
        "Topic :: Scientific/Engineering :: Astronomy",
        "Topic :: Scientific/Engineering :: Visualization",
    ],
)
