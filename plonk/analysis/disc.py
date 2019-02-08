'''
disc.py

Phantom analysis for dusty discs.

Daniel Mentiplay, 2019.
'''

import numpy as np
from numpy.linalg import norm

from ..constants import constants
from ..ParticleData import density_from_smoothing_length

# ---------------------------------------------------------------------------- #

def disc_analysis(dump, radiusIn, radiusOut, numberRadialBins,
                  midplaneSlice=False, minParticleAverage=5):
    '''
    Perform disc analysis.

    Arguments:
        dump             : Dump object
        radiusIn         : Inner disc radius for radial binning
        radiusOut        : Outer disc radius for radial binning
        numberRadialBins : Number of radial bins

    Optional:
        midplaneSlice      : Calculate midplane density by taking a slice
        minParticleAverage : Minimum number of particles to compute averages
    '''

    # TODO: add to docs

    particleData = dump.ParticleData
    sinkData = dump.SinkData
    parameters = dump.Parameters

    isFullDump = 'vx' in particleData.columns

    nDustSmall = parameters['ndustsmall']
    nDustLarge = parameters['ndustlarge']

    containsSmallDust = bool(nDustSmall > 0)
    containsLargeDust = bool(nDustLarge > 0)
    containsDust = bool(containsSmallDust or containsLargeDust)

    if containsDust:
        grainSize = parameters['grainsize']
        grainDens = parameters['graindens']

    nSinks = parameters['nptmass']
    containsSinks = bool(nSinks > 0)

#--- Units

    units = dump.Units

    uDist = units['distance']
    uTime = units['time']
    uMass = units['mass']

#--- Extra particle properties

    if isFullDump:

        particleData['px'] = particleData['m'] * particleData['vx']
        particleData['py'] = particleData['m'] * particleData['vy']
        particleData['pz'] = particleData['m'] * particleData['vz']

        angularMomentum = np.cross(particleData[['x', 'y', 'z']],
                                   particleData[['vx', 'vy', 'vz']])

        particleData['lx'] = angularMomentum[:, 0]
        particleData['ly'] = angularMomentum[:, 1]
        particleData['lz'] = angularMomentum[:, 2]

#--- Sink properties

    if containsSinks:

        # TODO: check if sink[0] is really the star; check if binary
        stellarMass = sinkData['m'][0]

        gravitationalParameter = constants.gravitationalConstant \
                               / ( uDist**3 / uTime**2 / uMass ) \
                               * stellarMass
    else:

        raise Exception('Must have at least one sink to do disc analysis')

    if isFullDump:

        sinkData['px'] = sinkData['m'] * sinkData['vx']
        sinkData['py'] = sinkData['m'] * sinkData['vy']
        sinkData['pz'] = sinkData['m'] * sinkData['vz']

        angularMomentum = np.cross(sinkData[['x', 'y', 'z']],
                                   sinkData[['vx', 'vy', 'vz']])

        sinkData['lx'] = angularMomentum[:, 0]
        sinkData['ly'] = angularMomentum[:, 1]
        sinkData['lz'] = angularMomentum[:, 2]

#--- Eccentricity

    if isFullDump:

        particleData['e'] = \
            _calculate_eccentricity( sinkData['m'],
                                     sinkData[['x', 'y', 'z']],
                                     sinkData[['vx', 'vy', 'vz']],
                                     sinkData[['lx', 'ly', 'lz']],
                                     gravitationalParameter )

        sinkData['e'] = \
            _calculate_eccentricity( sinkData['m'],
                                     sinkData[['x', 'y', 'z']],
                                     sinkData[['vx', 'vy', 'vz']],
                                     sinkData[['lx', 'ly', 'lz']],
                                     gravitationalParameter )

#--- Calculate radially binned quantities

    radialAverages = _calculate_radially_binned_quantities( numberRadialBins,
                                                            radiusIn,
                                                            radiusOut,
                                                            particleData,
                                                            isFullDump,
                                                            parameters,
                                                            midplaneSlice,
                                                            minParticleAverage )

#--- Stokes

    gamma = parameters.eos['gamma']

    Stokes = [np.full_like(radialBinsDisc, np.nan) for i in range(nDustLarge)]

    for idxi in range(len(radialBinsDisc)):
        for idxj in range(nDustLarge):

            Stokes[idxj][idxi] = \
                np.sqrt(gamma*np.pi/8) * grainDens[idxj] * grainSize[idxj] \
                / ( scaleHeightGas[idxi] \
                * (midplaneDensityGas[idxi] + midplaneDensityDust[idxj][idxi]) )

#--- Return

    return radialAverages

# ---------------------------------------------------------------------------- #

def _calculate_radially_binned_quantities( numberRadialBins=None,
                                           radiusIn=None,
                                           radiusOut=None,
                                           particleData=None,
                                           isFullDump=None,
                                           parameters=None,
                                           midplaneSlice=None,
                                           minParticleAverage=None ):
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

    if numberRadialBins is None:
        raise ValueError('Need numberRadialBins')

    if radiusIn is None:
        raise ValueError('Need radiusIn')

    if radiusOut is None:
        raise ValueError('Need radiusOut')

    if midplaneSlice and parameters is None:
        raise ValueError('"parameters" required to calculate midplane slice')

    for index, R in enumerate(radialBins):

        area = np.pi * ( (R + dR/2)**2 - (R - dR/2)**2 )

        indicies = np.where((cylindricalRadius < R + dR/2) & \
                            (cylindricalRadius > R - dR/2))[0]

        nPart = len(indicies)

        surfaceDensity[index] = massParticle * nPart / area

        if nPart > minParticleAverage:

            meanSmoothingLength[index] = np.sum(
                smoothingLength[indicies] ) / nPart

            meanHeight = np.sum( height[indicies] ) / nPart

            scaleHeight[index] = np.sqrt( np.sum(
                (height[indicies] - meanHeight)**2 ) / (nPart - 1) )

            if useVelocities:

                meanAngularMomentum[:, index] = np.sum(
                    angularMomentum[indicies], axis=0 ) / nPart

                magnitudeAngularMomentum[index] = \
                        norm(meanAngularMomentum[:, index])

                meanTilt[index] = np.arccos( meanAngularMomentum[2, index] \
                                           / magnitudeAngularMomentum[index] )

                meanTwist[index] = np.arctan2(
                    meanAngularMomentum[1, index] / magnitudeAngularMomentum[index],
                    meanAngularMomentum[0, index] / magnitudeAngularMomentum[index] )

                meanEccentricity[index] = np.sum( eccentricity[indicies] ) / nPart

        else:

            meanSmoothingLength[index] = np.nan

            scaleHeight[index] = np.nan

            if useVelocities:

                meanAngularMomentum[:, index] = np.nan
                meanTilt[index]               = np.nan
                meanTwist[index]              = np.nan
                meanEccentricity[index]       = np.nan

        if midplaneSlice:

            frac = 1/2

            indiciesMidplane = np.where(
                (cylindricalRadius < R + dR/2) &
                (cylindricalRadius > R - dR/2) &
                (height < meanHeight + frac * scaleHeight[index]) &
                (height > meanHeight - frac * scaleHeight[index])
                )[0]

            hfact = parameters.numerical['hfact']

            midplaneDensity[index] = np.sum(
                density_from_smoothing_length(
                    smoothingLength[indiciesMidplane], massParticle , hfact)
                ) / nPart

        else:

            midplaneDensity[index] = surfaceDensity[index] / np.sqrt(2*np.pi) \
                                                           / scaleHeight[index]

        midplaneDensity[index] = np.nan_to_num( midplaneDensity[index] )

    if useVelocities:

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

def _calculate_eccentricity( massParticle,
                             position,
                             velocity,
                             angularMomentum,
                             gravitationalParameter ):
    '''
    Calculate eccentricity.
    '''

    specificKineticEnergy = 1/2 * norm(velocity, axis=1)**2
    specificGravitationalEnergy = - gravitationalParameter / norm(position, axis=1)
    specificEnergy = specificKineticEnergy + specificGravitationalEnergy

    specificAngularMomentum = norm(angularMomentum, axis=1) \
                            / massParticle

    term = 2 * specificEnergy * specificAngularMomentum**2 \
         / gravitationalParameter**2

    eccentricity = np.sqrt( 1 + term )

    return eccentricity
