'''
disc.py

Daniel Mentiplay, 2019.
'''

import numpy as np
from numpy.linalg import norm

from ..constants import constants
from ..dumps import iGas, iDust
from ..utils import density_from_smoothing_length

#--- Options

midplaneSlice    = False  # Calculate midplane density by taking a slice
numberRadialBins = 150    # Number of radial bins
minPart          = 5      # Minimum number of particles to compute averages

#--- Parameters

rIn   = 10   # TODO: get as input (or calculate from data?)
rOut  = 200  # TODO: get as input (or calculate from data?)

# ---------------------------------------------------------------------------- #

def disc_analysis(dump):
    '''
    Perform disc analysis.
    '''
    # TODO: add docs

    isFullDump = bool(dump.dumpType == 'full')
    arrays = dump.arrays
    parameters = dump.parameters
    containsDust = dump.containsDust

#--- Units

    units = dump.units

    uDist = units['distance']
    uTime = units['time']
    uMass = units['mass']

#--- Gas particle properties

    massParticleGas = parameters['massoftype'][iGas - 1]

    smoothingLengthGas = arrays.smoothingLength['gas']
    positionGas = arrays.position['gas']

    if isFullDump:
        velocityGas = arrays.velocity['gas']
        momentumGas = massParticleGas * velocityGas
        angularMomentumGas = np.cross(positionGas, momentumGas)
    else:
        velocityGas = None
        momentumGas = None
        angularMomentumGas = None

#--- Dust particle properties

    if containsDust:

        if 'ndustlarge' in parameters:

            nDustLarge = parameters['ndustlarge']
            massParticleDust = parameters['massoftype'] \
                    [iDust - 1 : iDust + nDustLarge - 1]

        else:
            nDustLarge = 0

        grainDens = parameters['graindens'][0: nDustLarge]
        grainSize = parameters['grainsize'][0: nDustLarge]

        smoothingLengthDust = arrays.smoothingLength['dust']
        positionDust = arrays.position['dust']

        if isFullDump:
            velocityDust = arrays.velocity['dust']
            momentumDust = list()
            angularMomentumDust = list()
            for idx in range(nDustLarge):
                momentumDust.append(massParticleDust[idx] * velocityDust[idx])
                angularMomentumDust.append(np.cross(positionDust[idx],
                                                    momentumDust[idx]))
        else:
            velocityDust = None
            momentumDust = list()
            angularMomentumDust = list()
            for idx in range(nDustLarge):
                momentumDust.append(None)
                angularMomentumDust.append(None)

#--- Sink particle properties

    nSinks = parameters['nptmass']
    massParticleSink = arrays.massSink

    # TODO: check if sink[0] is really the star; check if binary
    stellarMass = massParticleSink[0]

    gravitationalParameter = constants.G / ( uDist**3 / uTime**2 / uMass ) \
                           * stellarMass

    smoothingLengthSink = arrays.smoothingLength['sink']
    positionSink = arrays.position['sink']

    if isFullDump:

        velocitySink = arrays.velocity['sink']
        momentumSink = list()
        angularMomentumSink = list()

        for idx in range(nSinks):
            momentumSink.append(massParticleSink[idx] * velocitySink[idx])
            angularMomentumSink.append(np.cross(positionSink[idx],
                                                momentumSink[idx]))

        eccentricitySink = _calculate_eccentricity( massParticleSink,
                                                    positionSink,
                                                    velocitySink,
                                                    angularMomentumSink,
                                                    gravitationalParameter )

#--- Gas eccentricity

    if isFullDump:

        eccentricityGas = _calculate_eccentricity( massParticleGas,
                                                   positionGas,
                                                   velocityGas,
                                                   angularMomentumGas,
                                                   gravitationalParameter )

    else:

        eccentricityGas = None

#--- Radially bin gas

    cylindricalRadiusGas = norm(positionGas[:, 0:2], axis=1)
    heightGas = positionGas[:, 2]

    vals = _calculate_radially_binned_quantities( numberRadialBins,
                                                  rIn,
                                                  rOut,
                                                  cylindricalRadiusGas,
                                                  heightGas,
                                                  smoothingLengthGas,
                                                  massParticleGas,
                                                  angularMomentumGas,
                                                  eccentricityGas )

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

    for idx in range(nDustLarge):

        #--- Dust eccentricity

        if isFullDump:

            eccentricityDust.append(
                _calculate_eccentricity( massParticleDust[idx],
                                         positionDust[idx],
                                         velocityDust[idx],
                                         angularMomentumDust[idx],
                                         gravitationalParameter ) )

        else:

            eccentricityDust.append(None)

        cylindricalRadiusDust.append(
            norm(arrays.position['dust'][idx][:, 0:2], axis=1) )

        heightDust.append(positionDust[idx][:, 2])

        vals = _calculate_radially_binned_quantities( numberRadialBins,
                                                      rIn,
                                                      rOut,
                                                      cylindricalRadiusDust[idx],
                                                      heightDust[idx],
                                                      smoothingLengthDust[idx],
                                                      massParticleDust[idx],
                                                      angularMomentumDust[idx],
                                                      eccentricityDust[idx] )

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

    gamma = parameters['gamma']

    Stokes = [np.empty_like(radialBinsDisc) for i in range(nDustLarge)]

    for idxi in range(len(radialBinsDisc)):
        for idxj in range(nDustLarge):

            Stokes[idxj][idxi] = \
                np.sqrt(gamma*np.pi/8) * grainDens[idxj] * grainSize[idxj] \
                / ( scaleHeightGas[idxi] \
                * (midplaneDensityGas[idxi] + midplaneDensityDust[idxj][idxi]) )

#--- Package into dictionary

    surfaceDensity              = dict()
    surfaceDensity['gas']       = surfaceDensityGas
    surfaceDensity['dust']      = surfaceDensityDust

    midplaneDensity             = dict()
    midplaneDensity['gas']      = midplaneDensityGas
    midplaneDensity['dust']     = midplaneDensityDust

    meanSmoothingLength         = dict()
    meanSmoothingLength['gas']  = meanSmoothingLengthGas
    meanSmoothingLength['dust'] = meanSmoothingLengthDust
    meanSmoothingLength['sink'] = smoothingLengthSink

    scaleHeight                 = dict()
    scaleHeight['gas']          = scaleHeightGas
    scaleHeight['dust']         = scaleHeightDust

    meanAngularMomentum         = dict()
    meanAngularMomentum['gas']  = meanAngularMomentumGas
    meanAngularMomentum['dust'] = meanAngularMomentumDust
    meanAngularMomentum['sink'] = angularMomentumSink

    tilt                        = dict()
    tilt['gas']                 = tiltGas
    tilt['dust']                = tiltDust

    twist                       = dict()
    twist['gas']                = twistGas
    twist['dust']               = twistDust

    psi                         = dict()
    psi['gas']                  = psiGas
    psi['dust']                 = psiDust

    meanEccentricity            = dict()
    meanEccentricity['gas']     = meanEccentricityGas
    meanEccentricity['dust']    = meanEccentricityDust
    meanEccentricity['sink']    = eccentricitySink


#--- Return

    return ( radialBinsDisc,
             surfaceDensity,
             midplaneDensity,
             meanSmoothingLength,
             scaleHeight,
             meanAngularMomentum,
             tilt,
             twist,
             psi,
             meanEccentricity,
             Stokes )

# ---------------------------------------------------------------------------- #

def _calculate_radially_binned_quantities( nRadialBins=None,
                                           radiusIn=None,
                                           radiusOut=None,
                                           cylindricalRadius=None,
                                           height=None,
                                           smoothingLength=None,
                                           massParticle=None,
                                           angularMomentum=None,
                                           eccentricity=None ):
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
    radialBins = np.linspace(radiusIn, radiusOut, nRadialBins)

    meanSmoothingLength = np.empty_like(radialBins)
    surfaceDensity      = np.empty_like(radialBins)
    midplaneDensity     = np.empty_like(radialBins)
    scaleHeight         = np.empty_like(radialBins)

    if angularMomentum is not None:

        if eccentricity is None:
            raise ValueError('Need eccentricity')

        useVelocities = True

        meanAngularMomentum      = np.empty_like(3*[radialBins])
        magnitudeAngularMomentum = np.empty_like(radialBins)
        meanTilt                 = np.empty_like(radialBins)
        meanTwist                = np.empty_like(radialBins)
        meanPsi                  = np.empty_like(radialBins)
        meanEccentricity         = np.empty_like(radialBins)

    else:

        useVelocities = False

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

            midplaneDensity[index] = np.sum(
                density_from_smoothing_length(
                    smoothingLength[indiciesMidplane], massParticle ) ) / nPart

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
