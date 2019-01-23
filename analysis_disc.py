'''
analysis_disc.py

Daniel Mentiplay, 2019.
'''

# TODO: really just a testing script at the moment.

import numpy as np

from dump import Dump
from utils import density_from_smoothing_length

#--- Options

midplaneSlice = False
nRadialBins = 150
minPart = 5

#--- Parameters

gamma = 1    # TODO: read from dump
rIn   = 10   # TODO: read from dump
rOut  = 200  # TODO: read from dump

#--- Read dump file

dump = Dump('disc_00000.ascii')
isFullDump = bool(dump.dumpType == 'full')

#--- Units

unitDist = dump.units['dist']
unitTime = dump.units['time']
unitMass = dump.units['mass']

unitMomen = unitMass * unitDist / unitTime
unitAngMomen = unitMass * unitDist**2 / unitTime
unitDens = unitMass / unitDist**3
unitSurfaceDens = unitMass / unitDist**2

#--- Gas particle properties

massParticleGas = dump.massParticles['gas']

smoothingLengthGas = dump.smoothingLength['gas']
positionGas = dump.position['gas']

if isFullDump:
    velocityGas = dump.velocity['gas']
    momentumGas = massParticleGas * velocityGas
    angularMomentumGas = np.cross(positionGas, momentumGas)

#--- Dust particle properties

nDustTypes = len(dump.nParticles['dust'])
massParticleDust = dump.massParticles['dust']
massParticleDust[1] = massParticleDust[0]  # TODO: hack for broken splash to ascii

# TODO: get from .setup or HDF5 file
grainDens = np.array([3., 3.]) / unitDens
grainSize = np.array([0.01, 0.1]) / unitDist

smoothingLengthDust = dump.smoothingLength['dust']
positionDust = dump.position['dust']

if isFullDump:
    velocityDust = dump.velocity['dust']
    momentumDust = list()
    angularMomentumDust = list()
    for idx in range(nDustTypes):
        momentumDust.append(massParticleDust[idx] * velocityDust[idx])
        angularMomentumDust.append(np.cross(positionDust[idx], momentumDust[idx]))

#--- Sink particle properties

nSinks = dump.nParticles['sink']
massParticleSink = dump.massParticles['sink']

smoothingLengthSink = dump.smoothingLength['sink']
positionSink = dump.position['sink']

if isFullDump:
    velocitySink = dump.velocity['sink']
    momentumSink = list()
    angularMomentumSink = list()
    for idx in range(nSinks):
        momentumSink.append(massParticleSink[idx] * velocitySink[idx])
        angularMomentumSink.append(np.cross(positionSink[idx], momentumSink[idx]))

#--- Radial binning

dR = (rOut - rIn)/(nRadialBins - 1)
radius = np.linspace(rIn, rOut, nRadialBins)

cylindricalRadiusGas = np.linalg.norm(positionGas[:, 0:2], axis=1)
heightGas = positionGas[:, 2]

cylindricalRadiusDust = list()
heightDust = list()
for idx in range(nDustTypes):
    cylindricalRadiusDust.append( \
            np.linalg.norm(dump.position['dust'][idx][:, 0:2], axis=1) )
    heightDust.append(positionDust[idx][:, 2])

#--- Calculate radially binned quantities

meanSmoothingLengthGas = np.empty_like(radius)
surfaceDensityGas      = np.empty_like(radius)
midplaneDensityGas     = np.empty_like(radius)
scaleHeightGas         = np.empty_like(radius)

if isFullDump:
    meanAngularMomentumGas = np.empty_like([radius, radius, radius])
    tilt  = np.empty_like(radius)
    twist = np.empty_like(radius)

meanSmoothingLengthDust = [np.empty_like(radius) for i in range(nDustTypes)]
surfaceDensityDust      = [np.empty_like(radius) for i in range(nDustTypes)]
midplaneDensityDust     = [np.empty_like(radius) for i in range(nDustTypes)]
scaleHeightDust         = [np.empty_like(radius) for i in range(nDustTypes)]
Stokes                  = [np.empty_like(radius) for i in range(nDustTypes)]

if isFullDump:
    meanAngularMomentumDust = [np.empty_like([radius, radius, radius])
                               for i in range(nDustTypes)]


for idxi, R in enumerate(radius):

    #--- Gas

    area = np.pi * ( (R + dR/2)**2 - (R - dR/2)**2 )

    indiciesGas = np.where((cylindricalRadiusGas < R + dR/2) & \
                           (cylindricalRadiusGas > R - dR/2))[0]

    nPartGas = len(indiciesGas)

    meanSmoothingLengthGas[idxi] = np.sum( smoothingLengthGas[indiciesGas] ) \
                                         / nPartGas

    surfaceDensityGas[idxi] = massParticleGas * nPartGas / area

    meanHeightGas = np.sum( heightGas[indiciesGas] ) / nPartGas

    scaleHeightGas[idxi] = np.sqrt( np.sum( \
                           (heightGas[indiciesGas] - meanHeightGas)**2 ) \
                                 / (nPartGas - 1) )

    if isFullDump:

        meanAngularMomentumGas[:, idxi] = np.sum( angularMomentumGas[indiciesGas],
                                                  axis=0 ) / nPartGas

        magnitudeAngularMomentumGas = np.linalg.norm(meanAngularMomentumGas[:, idxi])

        tilt[idxi] = np.arccos( meanAngularMomentumGas[2, idxi] \
                              / magnitudeAngularMomentumGas )

        twist[idxi] = np.arctan2(
            meanAngularMomentumGas[1, idxi] / magnitudeAngularMomentumGas,
            meanAngularMomentumGas[0, idxi] / magnitudeAngularMomentumGas )

    if midplaneSlice:

        frac = 1/2

        indiciesMidplaneGas = np.where(
            (cylindricalRadiusGas < R + dR/2) &
            (cylindricalRadiusGas > R - dR/2) &
            (heightGas < meanHeightGas + frac * scaleHeightGas[idxi]) &
            (heightGas > meanHeightGas - frac * scaleHeightGas[idxi])
            )[0]

        midplaneDensityGas[idxi] = np.sum(
            density_from_smoothing_length(
                smoothingLengthGas[indiciesMidplaneGas], massParticleGas ) ) \
                / nPartGas

    else:

        midplaneDensityGas[idxi] = surfaceDensityGas[idxi] / np.sqrt(2*np.pi) \
                                 / scaleHeightGas[idxi]

    #--- Dust

    for idxj in range(nDustTypes):

        indiciesDust = np.where((cylindricalRadiusDust[idxj] < R + dR/2) & \
                                (cylindricalRadiusDust[idxj] > R - dR/2))[0]

        nPartDust = len(indiciesDust)

        surfaceDensityDust[idxj][idxi] = massParticleDust[idxj] \
                                       * nPartDust / area

        if nPartDust > minPart:

            meanSmoothingLengthDust[idxj][idxi] = np.sum(
                smoothingLengthDust[idxj][indiciesDust] ) / nPartDust

            if isFullDump:
                meanAngularMomentumDust[idxj][:, idxi] = np.sum(
                    angularMomentumDust[idxj][indiciesDust], axis=0 ) / nPartDust

            meanHeightDust = np.sum(heightDust[idxj][indiciesDust]) / nPartDust

            scaleHeightDust[idxj][idxi] = np.sqrt(
                np.sum(
                    ( heightDust[idxj][indiciesDust] - meanHeightDust )**2
                    ) / (nPartDust - 1) )
        else:

            meanSmoothingLengthDust[idxj][idxi] = np.nan

            if isFullDump:
                meanAngularMomentumDust[idxj][:, idxi] = np.nan

            scaleHeightDust[idxj][idxi] = np.nan

        midplaneDensityDust[idxj][idxi] = surfaceDensityDust[idxj][idxi] \
                                        / np.sqrt(2*np.pi) \
                                        / scaleHeightDust[idxj][idxi]

        midplaneDensityDust[idxj][idxi] = np.nan_to_num(
            midplaneDensityDust[idxj][idxi] )

        Stokes[idxj][idxi] = \
            np.sqrt(gamma*np.pi/8) * grainDens[idxj] * grainSize[idxj] \
            / ( scaleHeightGas[idxi] \
            * (midplaneDensityGas[idxi] + midplaneDensityDust[idxj][idxi]) )
