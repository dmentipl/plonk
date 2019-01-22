'''
analysis_disc.py

Daniel Mentiplay, 2019.
'''

# TODO: really just a testing script at the moment.

import numpy as np

from dump import Dump
from utils import density_from_smoothing_length

#--- Read dump file

dump = Dump('disc_00000.ascii')

#--- Units

udist = dump.units['dist']
utime = dump.units['time']
umass = dump.units['mass']

#--- Particle properties

massGasParticle = dump.massParticles['gas']
position        = dump.position['gas']
velocity        = dump.velocity['gas']
momentum        = massGasParticle * velocity
angularMomentum = np.cross(position, momentum)
smoothingLength = dump.smoothingLength['gas']

#--- Radial binning

nBins = 150
rIn = 10
rOut = 200

dR = (rOut - rIn)/(nBins - 1)
radius = np.linspace(rIn, rOut, nBins)
cylindricalRadius = np.linalg.norm(dump.position['gas'][:,0:2], axis=1)

#--- Surface density

surfaceDensity  = np.empty_like(radius)
averageDensity  = np.empty_like(radius)
meanHeight      = np.empty_like(radius)
scaleHeight     = np.empty_like(radius)

for idx, R in enumerate(radius):

    indicies = np.where((cylindricalRadius < R + dR/2) & \
                        (cylindricalRadius > R - dR/2))[0]

    npart = len(indicies)

    area = np.pi * ( (R + dR/2)**2 - (R - dR/2)**2 )
    surfaceDensity[idx] = massGasParticle * npart / area

    meanHeight[idx]  = np.sum(position[indicies, 2])/npart
    scaleHeight[idx] = np.sqrt(np.sum( ( position[indicies, 2] \
                                       - meanHeight[idx] )**2 ) / (npart - 1))

    averageDensity[idx] = np.sum( density_from_smoothing_length( \
            smoothingLength[indicies], massGasParticle ) ) / npart




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
