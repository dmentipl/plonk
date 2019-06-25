#
# build.sh
# Install script for conda package.
#
# Daniel Mentiplay, 2019
#

#--- Compile Splash Fortran source

make conda

#--- Plonk Python package

python setup.py install --single-version-externally-managed --record=record.txt
