'''
libsplash.py

Daniel Mentiplay, 2019.
'''

import subprocess

files = ['interpolate3D_projection.F90']

for file in files:
    subprocess.run('gfortran -c ' + file, shell=True)

subprocess.run('f90wrap -m interpolate3D_projection ' \
              + 'interpolate3D_projection.f90', shell=True)
subprocess.run('f2py-f90wrap -m _interpolate3D_projection *.o -c ' \
               + 'f90wrap_interpolate3D_projection.f90', shell=True)

with open('interpolate3D_projection.py', 'r') as file:
    lines = file.readlines()

lines[0] = 'from .libsplash import _interpolate3D_projection\n'

with open('interpolate3D_projection.py', 'w') as file:
    file.writelines(lines)

subprocess.run('mv interpolate3D_projection.py ..', shell=True)
subprocess.run('rm *.mod *.o f90wrap_*.f90', shell=True)
