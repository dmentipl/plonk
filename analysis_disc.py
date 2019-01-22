'''
analysis_disc.py

Daniel Mentiplay, 2019.
'''

# TODO: really just a testing script at the moment.

import numpy as np

from dump import Dump

dump = Dump('disc_00000.ascii')

udist = dump.units['dist']
utime = dump.units['time']
umass = dump.units['mass']

#--- Surface density

mpGas = dump.massParticles['gas']
cylindricalRadius = np.linalg.norm(dump.position['gas'][:,0:2], axis=1)

nbins = 300

rIn = 10
rOut = 200
dR = (rOut - rIn)/(nbins - 1)

radii = np.linspace(rIn, rOut, nbins)
sigma = np.empty_like(radii)

for idx, R in enumerate(radii):

    npart1 = len(np.where(cylindricalRadius < R + dR/2)[0])
    npart2 = len(np.where(cylindricalRadius < R - dR/2)[0])
    npart = abs(npart1 - npart2)
    area = np.pi * ( (R + dR/2)**2 - (R - dR/2)**2 )
    sigma[idx] = mpGas * npart / area

#--- Angular momentum

position = dump.position['gas']
momentum = dump.massParticles['gas'] * dump.velocity['gas']

angularMomentum = np.cross(position, momentum)
angularMomentumTotal = np.sum(angularMomentum, axis=0)
angularMomentumTotalMagnitude = np.linalg.norm(angularMomentumTotal)

angmom = angularMomentumTotalMagnitude * umass * udist**2 / utime
print(f'Total angular momentum in gas disc is {angmom:g} g.cm^2/s')

nSinks = dump.nParticles['sink']
position = dump.position['sink']
momentum = np.empty_like(position)
for i in range(nSinks):
    momentum[i] = dump.massParticles['sink'][i] * dump.velocity['sink'][i]

angularMomentum = np.cross(position, momentum)
angularMomentumTotalMagnitude = np.sum(angularMomentum, axis=1)
angularMomentumTotalMagnitudeCGS = angularMomentumTotalMagnitude * umass * udist**2 / utime
for idx, angmom in enumerate(angularMomentumTotalMagnitudeCGS):
    print(f'Total angular momentum in sink {idx} is {angmom:g} g.cm^2/s')
