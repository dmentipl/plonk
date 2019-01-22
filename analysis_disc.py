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

unitDist = dump.units['dist']
unitTime = dump.units['time']
unitMass = dump.units['mass']

unitDens = unitMass / unitDist**3
unitSurfaceDens = unitMass / unitDist**2

#--- Particle properties

massParticleGas = dump.massParticles['gas']
massParticleDust = dump.massParticles['dust']
massParticleDust[1] = massParticleDust[0]  # TODO: hack for broken splash to ascii
massParticleSink = dump.massParticles['sink']

positionGas = dump.position['gas']
velocityGas = dump.velocity['gas']
momentumGas = massParticleGas * velocityGas

angularMomentumGas = np.cross(positionGas, momentumGas)
smoothingLengthGas = dump.smoothingLength['gas']

# TODO: get from .setup or HDF5 file
grainDens = np.array([3., 3.]) / unitDens
grainSize = np.array([0.01, 0.1]) / unitDist

#--- Radial binning

nBins = 150
rIn = 10
rOut = 200

dR = (rOut - rIn)/(nBins - 1)
radius = np.linspace(rIn, rOut, nBins)

cylindricalRadiusGas = np.linalg.norm(dump.position['gas'][:,0:2], axis=1)

nDustTypes = len(dump.nParticles['dust'])
cylindricalRadiusDust = list()
for i in range(nDustTypes):
    cylindricalRadiusDust.append(np.linalg.norm(dump.position['dust'][i][:,0:2], axis=1))

#--- Surface density

surfaceDensityGas  = np.empty_like(radius)
averageDensityGas  = np.empty_like(radius)
meanHeightGas      = np.empty_like(radius)
scaleHeightGas     = np.empty_like(radius)

surfaceDensityDust = [np.empty_like(radius) for i in range(nDustTypes)]

for idx, R in enumerate(radius):

    area = np.pi * ( (R + dR/2)**2 - (R - dR/2)**2 )

    indiciesGas = np.where((cylindricalRadiusGas < R + dR/2) & \
                        (cylindricalRadiusGas > R - dR/2))[0]

    nPartGas = len(indiciesGas)

    surfaceDensityGas[idx] = massParticleGas * nPartGas / area

    meanHeightGas[idx]  = np.sum(positionGas[indiciesGas, 2])/nPartGas
    scaleHeightGas[idx] = np.sqrt(np.sum( ( positionGas[indiciesGas, 2] \
                                       - meanHeightGas[idx] )**2 ) / (nPartGas - 1))

    averageDensityGas[idx] = np.sum( density_from_smoothing_length( \
            smoothingLengthGas[indiciesGas], massParticleGas ) ) / nPartGas

    for i in range(nDustTypes):

        indiciesDust = np.where((cylindricalRadiusDust[i] < R + dR/2) & \
                            (cylindricalRadiusDust[i] > R - dR/2))[0]

        nPartDust = len(indiciesDust)

        surfaceDensityDust[i][idx] = massParticleDust[i] * nPartDust / area



#--- Angular momentum

position = dump.position['gas']
momentum = dump.massParticles['gas'] * dump.velocity['gas']

angularMomentum = np.cross(position, momentum)
angularMomentumTotal = np.sum(angularMomentum, axis=0)
angularMomentumTotalMagnitude = np.linalg.norm(angularMomentumTotal)

angmom = angularMomentumTotalMagnitude * unitMass * unitDist**2 / unitTime
print(f'Total angular momentum in gas disc is {angmom:g} g.cm^2/s')

nSinks = dump.nParticles['sink']
position = dump.position['sink']
momentum = np.empty_like(position)
for i in range(nSinks):
    momentum[i] = dump.massParticles['sink'][i] * dump.velocity['sink'][i]

angularMomentum = np.cross(position, momentum)
angularMomentumTotalMagnitude = np.sum(angularMomentum, axis=1)
angularMomentumTotalMagnitudeCGS = angularMomentumTotalMagnitude * unitMass * unitDist**2 / unitTime
for idx, angmom in enumerate(angularMomentumTotalMagnitudeCGS):
    print(f'Total angular momentum in sink {idx} is {angmom:g} g.cm^2/s')
