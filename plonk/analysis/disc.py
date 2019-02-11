'''
disc.py

Phantom analysis for dusty discs.

Daniel Mentiplay, 2019.
'''

import numpy as np
from numpy.linalg import norm
import pandas as pd

from ..constants import constants

# ---------------------------------------------------------------------------- #

def disc_analysis(radiusIn=None,
                  radiusOut=None,
                  numberRadialBins=None,
                  particleData=None,
                  sinkData=None,
                  parameters=None,
                  units=None,
                  minParticleAverage=None):
    '''
    Perform disc analysis.

    TODO: add more docs

    Arguments:
        radiusIn         : Inner disc radius for radial binning
        radiusOut        : Outer disc radius for radial binning
        numberRadialBins : Number of radial bins

    Optional:
        minParticleAverage : Minimum number of particles to compute averages
    '''

    # TODO: add to docs

#--- Dump type

    isFullDump = 'vx' in particleData.columns

#--- Units

    uDist = units['distance']
    uTime = units['time']
    uMass = units['mass']

#--- Dust

    nDustSmall = parameters['ndustsmall']
    nDustLarge = parameters['ndustlarge']

    containsSmallDust = bool(nDustSmall > 0)
    containsLargeDust = bool(nDustLarge > 0)
    containsDust = bool(containsSmallDust or containsLargeDust)

    if containsDust:
        grainSize = parameters['grainsize']
        grainDens = parameters['graindens']

#--- Momentum and angular momentum

    if isFullDump:

        particleData['px'] = particleData['m'] * particleData['vx']
        particleData['py'] = particleData['m'] * particleData['vy']
        particleData['pz'] = particleData['m'] * particleData['vz']

        angularMomentum = np.cross(particleData[['x', 'y', 'z']],
                                   particleData[['vx', 'vy', 'vz']])

        particleData['lx'] = angularMomentum[:, 0]
        particleData['ly'] = angularMomentum[:, 1]
        particleData['lz'] = angularMomentum[:, 2]

    else:

        particleData['px'] = np.nan
        particleData['py'] = np.nan
        particleData['pz'] = np.nan
        particleData['lx'] = np.nan
        particleData['ly'] = np.nan
        particleData['lz'] = np.nan

#--- Sinks

    nSinks = parameters['nptmass']
    containsSinks = bool(nSinks > 0)

    if containsSinks:

        # TODO: check if sink[0] is really the star; check if binary
        stellarMass = sinkData['m'][0]
        print('Assuming the first sink particle is the central star')

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

    else:

        sinkData['px'] = np.nan
        sinkData['py'] = np.nan
        sinkData['pz'] = np.nan
        sinkData['lx'] = np.nan
        sinkData['ly'] = np.nan
        sinkData['lz'] = np.nan

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

    else:

        particleData['e'] = np.nan
        sinkData['e']     = np.nan

#--- Calculate radially binned quantities

    particleData['R'] = norm(particleData[['x', 'y']], axis=1)

    radialData = _calculate_radially_binned_quantities( numberRadialBins,
                                                        radiusIn,
                                                        radiusOut,
                                                        particleData,
                                                        minParticleAverage )

#--- Stokes

    gamma = parameters['gamma']

    # for idxi in range(len(radialBinsDisc)):
    #     for idxj in range(nDustLarge):

    #         Stokes[idxj][idxi] = \
    #             np.sqrt(gamma*np.pi/8) * grainDens[idxj] * grainSize[idxj] \
    #             / ( scaleHeightGas[idxi] \
    #             * (midplaneDensityGas[idxi] + midplaneDensityDust[idxj][idxi]) )

#--- Return

    return radialData, particleData, sinkData

# ---------------------------------------------------------------------------- #

def _calculate_radially_binned_quantities( numberRadialBins=None,
                                           radiusIn=None,
                                           radiusOut=None,
                                           particleData=None,
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

    if particleData is None:
        raise ValueError('Need particleData')

    if minParticleAverage is None:
        minParticleAverage = 5

    radialBinWidth = (radiusOut - radiusIn) / (numberRadialBins - 1)
    radialBins     = np.linspace(radiusIn, radiusOut, numberRadialBins)

    radialData = pd.DataFrame(radialBins, columns=['R'])
    radialData['area'] = np.pi * ( (radialBins + radialBinWidth/2)**2 \
                                 - (radialBins - radialBinWidth/2)**2 )

    radialData = radialData.reindex(
        columns = radialData.columns.tolist() \
        + ['sigma', 'h', 'H', 'lx', 'ly', 'lz', 'l', 'tilt', 'twist', 'e'] )

    for index, radius in enumerate(radialBins):

        radiusLeft  = radius - radialBinWidth/2
        radiusRight = radius + radialBinWidth/2

        particlesInRadialBin = \
            particleData.loc[ (particleData['R'] >= radiusLeft) \
                            & (particleData['R'] <= radiusRight) ]

        radialData['sigma'].iloc[index] = particlesInRadialBin['m'].sum() \
                                        / radialData['area'].iloc[index]

        if len(particlesInRadialBin) > minParticleAverage:

            radialData['h'].iloc[index] = particlesInRadialBin['h'].mean()
            radialData['H'].iloc[index] = particlesInRadialBin['z'].std()

            radialData['lx'].iloc[index] = particlesInRadialBin['lx'].mean()
            radialData['ly'].iloc[index] = particlesInRadialBin['ly'].mean()
            radialData['lz'].iloc[index] = particlesInRadialBin['lz'].mean()

            radialData['e'].iloc[index] = particlesInRadialBin['e'].mean()


    radialData['l'] = norm(radialData[['lx','ly','lz']].values, axis=1)

    radialData['tilt'] = np.arccos( radialData['lz'] / radialData['l'] )

    radialData['twist'] = np.arctan2( radialData['ly'] / radialData['l'],
                                      radialData['lx'] / radialData['l'] )

    radialData['rho'] = radialData['sigma'] / radialData['H'] / np.sqrt(2*np.pi)

    radialData['psi'] = radialBins * np.sqrt(
        np.gradient( radialData['lx']/radialData['l'], radialBinWidth )**2 + \
        np.gradient( radialData['ly']/radialData['l'], radialBinWidth )**2 + \
        np.gradient( radialData['lz']/radialData['l'], radialBinWidth )**2 )

    return radialData

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
