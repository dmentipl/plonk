'''
libsplash.py

Daniel Mentiplay, 2019.
'''

import subprocess

compileSplash = 'gfortran -c splash.f90'
f90wrapSplash = 'f90wrap -m splash splash.f90'
f2pySplash    = 'f2py-f90wrap -m _splash splash.o -c f90wrap_splash.f90'

subprocess.run(compileSplash, shell=True)
subprocess.run(f90wrapSplash, shell=True)
subprocess.run(f2pySplash,    shell=True)

with open('splash.py', 'r') as file:
    lines = file.readlines()

lines[0] = 'from .libsplash import _splash\n'

with open('splash.py', 'w') as file:
    file.writelines(lines)

subprocess.run('mv splash.py ..', shell=True)
