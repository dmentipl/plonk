'''
splash.py

Daniel Mentiplay, 2019.
'''

try:
    from ._splash import splash
except ImportError:
    print('Cannot import splash module.')
    print('Need to compile Splash Fortran interpolation subroutines.')
    print('See README.md.')
