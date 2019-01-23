'''
analysis_disc.py

Daniel Mentiplay, 2019.
'''

# TODO: really just a testing script at the moment.

import numpy as np
from numpy.linalg import norm

from constants import constants
from dump import Dump
from utils import density_from_smoothing_length

#--- Options

midplaneSlice    = False
numberRadialBins = 150
minPart          = 5

#--- Parameters

gamma = 1    # TODO: read from dump
rIn   = 10   # TODO: read from dump, or set manually, or calculate
rOut  = 200  # TODO: read from dump, or set manually, or calculate

#--- Dump file name

dumpFileName = 'disc_00000.ascii'  # TODO: get dump filename as input
                                   # TODO: read multiple dumpfiles

# ---------------------------------------------------------------------------- #

def calculate_radially_binned_quantities( nRadialBins=None,
                                          radiusIn=None,
                                          radiusOut=None,
                                          cylindricalRadius=None,
                                          height=None,
                                          smoothingLength=None,
                                          massParticle=None,
                                          angularMomentum=None ):
    '''
    Calculate averaged radially binned quantities:
        - radial bins
        - surface density
        - midplane density
        - scale height
        - smoothing length
        - angular momentum
        - tilt
        - twist
        - psi
        - eccentricity
    '''

    if nRadialBins is None:
        raise ValueError('Need nRadialBins')

    if radiusIn is None:
        raise ValueError('Need radiusIn')

    if radiusOut is None:
        raise ValueError('Need radiusOut')

    if cylindricalRadius is None:
        raise ValueError('Need cylindricalRadius')

    if height is None:
        raise ValueError('Need height')

    if smoothingLength is None:
        raise ValueError('Need smoothingLength')

    if massParticle is None:
        raise ValueError('Need massParticle')

    dR         = (radiusOut - radiusIn) / (nRadialBins - 1)
    radialBins = np.linspace(rIn, rOut, nRadialBins)

    meanSmoothingLength = np.empty_like(radialBins)
    surfaceDensity      = np.empty_like(radialBins)
    midplaneDensity     = np.empty_like(radialBins)
    scaleHeight         = np.empty_like(radialBins)

    if angularMomentum is not None:

        meanAngularMomentum      = np.empty_like(3*[radialBins])
        magnitudeAngularMomentum = np.empty_like(radialBins)
        meanTilt                 = np.empty_like(radialBins)
        meanTwist                = np.empty_like(radialBins)
        meanPsi                  = np.empty_like(radialBins)
        meanEccentricity         = np.empty_like(radialBins)

    else:

        meanAngularMomentum      = None
        magnitudeAngularMomentum = None
        meanTilt                 = None
        meanTwist                = None
        meanPsi                  = None
        meanEccentricity         = None


    for index, R in enumerate(radialBins):

        area = np.pi * ( (R + dR/2)**2 - (R - dR/2)**2 )

        indicies = np.where((cylindricalRadius < R + dR/2) & \
                            (cylindricalRadius > R - dR/2))[0]

        nPart = len(indicies)

        surfaceDensity[index] = massParticle * nPart / area

        if nPart > minPart:

            meanSmoothingLength[index] = np.sum(
                smoothingLength[indicies] ) / nPart

            meanHeight = np.sum( height[indicies] ) / nPart

            scaleHeight[index] = np.sqrt( np.sum(
                (height[indicies] - meanHeight)**2 ) / (nPart - 1) )

            if angularMomentum is not None:

                meanAngularMomentum[:, index] = np.sum(
                    angularMomentum[indicies], axis=0 ) / nPart

                magnitudeAngularMomentum[index] = \
                        norm(meanAngularMomentum[:, index])

                meanTilt[index] = np.arccos( meanAngularMomentum[2, index] \
                                           / magnitudeAngularMomentum[index] )

                meanTwist[index] = np.arctan2(
                    meanAngularMomentum[1, index] / magnitudeAngularMomentum[index],
                    meanAngularMomentum[0, index] / magnitudeAngularMomentum[index] )

        else:

            meanSmoothingLength[index] = np.nan

            scaleHeight[index] = np.nan

            if angularMomentum is not None:

                meanAngularMomentum[:, index] = np.nan
                meanTilt[index]               = np.nan
                meanTwist[index]              = np.nan

        if midplaneSlice:

            frac = 1/2

            indiciesMidplane = np.where(
                (cylindricalRadius < R + dR/2) &
                (cylindricalRadius > R - dR/2) &
                (height < meanHeight + frac * scaleHeight[index]) &
                (height > meanHeight - frac * scaleHeight[index])
                )[0]

            midplaneDensity[index] = np.sum(
                density_from_smoothing_length(
                    smoothingLength[indiciesMidplane], massParticle ) ) / nPart

        else:

            midplaneDensity[index] = surfaceDensity[index] / np.sqrt(2*np.pi) \
                                                           / scaleHeight[index]

        midplaneDensity[index] = np.nan_to_num( midplaneDensity[index] )

    if angularMomentum is not None:

        unitAngularMomentum = meanAngularMomentum/magnitudeAngularMomentum

        meanPsi = radialBins * np.sqrt(
            np.gradient(unitAngularMomentum[0], dR)**2 + \
            np.gradient(unitAngularMomentum[1], dR)**2 + \
            np.gradient(unitAngularMomentum[2], dR)**2 )

    return ( radialBins,
             surfaceDensity,
             midplaneDensity,
             meanSmoothingLength,
             scaleHeight,
             meanAngularMomentum,
             meanTilt,
             meanTwist,
             meanPsi,
             meanEccentricity )

# ---------------------------------------------------------------------------- #

def calculate_eccentricity( massParticle,
                            position,
                            velocity,
                            angularMomentum ):
    '''
    Calculate eccentricity.
    '''

    gravitationalParameter = constants.G / ( unitDist**3 / unitTime**2 / unitMass )
    gravitationalParameter *= stellarMass

    specificKineticEnergy = 1/2 * norm(velocity, axis=1)**2
    specificGravitationalEnergy = - gravitationalParameter / norm(position, axis=1)
    specificEnergy = specificKineticEnergy + specificGravitationalEnergy

    specificAngularMomentum = norm(angularMomentum, axis=1) \
                            / massParticle

    term = 2 * specificEnergy * specificAngularMomentum**2 \
         / gravitationalParameter**2

    eccentricity = np.sqrt( 1 + term )

    return eccentricity

# ---------------------------------------------------------------------------- #

if __name__ == '__main__':

#--- Read dump file

    dump = Dump(dumpFileName)
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
    else:
        velocityGas = None
        momentumGas = None
        angularMomentumGas = None

#--- Dust particle properties

    nDustTypes = len(dump.nParticles['dust'])
    massParticleDust = dump.massParticles['dust']

    # TODO: hack for broken splash to ascii
    massParticleDust[1] = massParticleDust[0]

    # TODO: hack; get from .setup or HDF5 file
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
            angularMomentumDust.append(np.cross(positionDust[idx],
                                                momentumDust[idx]))
    else:
        velocityDust = None
        momentumDust = list()
        angularMomentumDust = list()
        for idx in range(nDustTypes):
            momentumDust.append(None)
            angularMomentumDust.append(None)

#--- Sink particle properties

    nSinks = dump.nParticles['sink']
    massParticleSink = dump.massParticles['sink']

    stellarMass = massParticleSink[0]

    smoothingLengthSink = dump.smoothingLength['sink']
    positionSink = dump.position['sink']

    if isFullDump:
        velocitySink = dump.velocity['sink']
        momentumSink = list()
        angularMomentumSink = list()
        for idx in range(nSinks):
            momentumSink.append(massParticleSink[idx] * velocitySink[idx])
            angularMomentumSink.append(np.cross(positionSink[idx],
                                                momentumSink[idx]))

#--- Gas eccentricity

    eccentricityGas = calculate_eccentricity( massParticleGas,
                                              positionGas,
                                              velocityGas,
                                              angularMomentumGas )

#--- Radially bin gas

    cylindricalRadiusGas = norm(positionGas[:, 0:2], axis=1)
    heightGas = positionGas[:, 2]

    vals = calculate_radially_binned_quantities( numberRadialBins,
                                                 rIn,
                                                 rOut,
                                                 cylindricalRadiusGas,
                                                 heightGas,
                                                 smoothingLengthGas,
                                                 massParticleGas,
                                                 angularMomentumGas )

    radialBinsDisc         = vals[0]
    surfaceDensityGas      = vals[1]
    midplaneDensityGas     = vals[2]
    meanSmoothingLengthGas = vals[3]
    scaleHeightGas         = vals[4]
    meanAngularMomentumGas = vals[5]
    tiltGas                = vals[6]
    twistGas               = vals[7]
    psiGas                 = vals[8]
    meanEccentricityGas    = vals[9]

#--- Radially bin dust

    cylindricalRadiusDust   = list()
    heightDust              = list()

    surfaceDensityDust      = list()
    midplaneDensityDust     = list()
    meanSmoothingLengthDust = list()
    scaleHeightDust         = list()
    meanAngularMomentumDust = list()
    tiltDust                = list()
    twistDust               = list()
    psiDust                 = list()
    eccentricityDust        = list()
    meanEccentricityDust    = list()

    for idx in range(nDustTypes):

        #--- Dust eccentricity

        eccentricityDust.append(
            calculate_eccentricity( massParticleDust[idx],
                                    positionDust[idx],
                                    velocityDust[idx],
                                    angularMomentumDust[idx] ) )

        cylindricalRadiusDust.append(
            norm(dump.position['dust'][idx][:, 0:2], axis=1) )

        heightDust.append(positionDust[idx][:, 2])

        vals = calculate_radially_binned_quantities( numberRadialBins,
                                                     rIn,
                                                     rOut,
                                                     cylindricalRadiusDust[idx],
                                                     heightDust[idx],
                                                     smoothingLengthDust[idx],
                                                     massParticleDust[idx],
                                                     angularMomentumDust[idx] )

        surfaceDensityDust.     append( vals[1] )
        midplaneDensityDust.    append( vals[2] )
        meanSmoothingLengthDust.append( vals[3] )
        scaleHeightDust.        append( vals[4] )
        meanAngularMomentumDust.append( vals[5] )
        tiltDust.               append( vals[6] )
        twistDust.              append( vals[7] )
        psiDust.                append( vals[8] )
        meanEccentricityDust.   append( vals[9] )

    #--- Stokes

    Stokes = [np.empty_like(radialBinsDisc) for i in range(nDustTypes)]

    for idxi in range(len(radialBinsDisc)):
        for idxj in range(nDustTypes):

            Stokes[idxj][idxi] = \
                np.sqrt(gamma*np.pi/8) * grainDens[idxj] * grainSize[idxj] \
                / ( scaleHeightGas[idxi] \
                * (midplaneDensityGas[idxi] + midplaneDensityDust[idxj][idxi]) )
