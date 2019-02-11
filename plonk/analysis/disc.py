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

#--- Spherical and cylindrical distance

    particleData['r'] = norm(particleData[['x', 'y', 'z']], axis=1)
    particleData['R'] = norm(particleData[['x', 'y']], axis=1)

#--- Momentum and angular momentum

    if isFullDump:

        particleData['|v|'] = norm(particleData[['vx', 'vy', 'vz']], axis=1)

        particleData['px'] = particleData['m'] * particleData['vx']
        particleData['py'] = particleData['m'] * particleData['vy']
        particleData['pz'] = particleData['m'] * particleData['vz']

        angularMomentum = np.cross(particleData[['x', 'y', 'z']],
                                   particleData[['px', 'py', 'pz']])

        particleData['Lx'] = angularMomentum[:, 0]
        particleData['Ly'] = angularMomentum[:, 1]
        particleData['Lz'] = angularMomentum[:, 2]
        particleData['|L|'] = norm(angularMomentum, axis=1)
        particleData['|l|'] = particleData['|L|'] / particleData['m']

    else:

        particleData['|v|'] = np.nan
        particleData['px']  = np.nan
        particleData['py']  = np.nan
        particleData['pz']  = np.nan
        particleData['Lx']  = np.nan
        particleData['Ly']  = np.nan
        particleData['Lz']  = np.nan
        particleData['|L|'] = np.nan
        particleData['|l|'] = np.nan

#--- Sinks

    nSinks = parameters['nptmass']
    containsSinks = bool(nSinks > 0)

    if containsSinks:

        sinkData['r'] = norm(sinkData[['x', 'y', 'z']], axis=1)
        sinkData['R'] = norm(sinkData[['x', 'y']], axis=1)

        # TODO: check if sink[0] is really the star; check if binary
        stellarMass = sinkData['m'][0]
        print('Assuming the first sink particle is the central star')

        gravitationalParameter = constants.gravitationalConstant \
                               / ( uDist**3 / uTime**2 / uMass ) \
                               * stellarMass
    else:

        raise Exception('Must have at least one sink to do disc analysis')

    if isFullDump:

        sinkData['|v|'] = norm(sinkData[['vx', 'vy', 'vz']], axis=1)

        sinkData['px'] = sinkData['m'] * sinkData['vx']
        sinkData['py'] = sinkData['m'] * sinkData['vy']
        sinkData['pz'] = sinkData['m'] * sinkData['vz']

        angularMomentum = np.cross(sinkData[['x', 'y', 'z']],
                                   sinkData[['px', 'py', 'pz']])

        sinkData['Lx'] = angularMomentum[:, 0]
        sinkData['Ly'] = angularMomentum[:, 1]
        sinkData['Lz'] = angularMomentum[:, 2]
        sinkData['|L|'] = norm(angularMomentum, axis=1)
        sinkData['|l|'] = sinkData['|L|'] / sinkData['m']

    else:

        sinkData['|v|'] = np.nan
        sinkData['px']  = np.nan
        sinkData['py']  = np.nan
        sinkData['pz']  = np.nan
        sinkData['Lx']  = np.nan
        sinkData['Ly']  = np.nan
        sinkData['Lz']  = np.nan
        sinkData['|L|'] = np.nan
        sinkData['|l|'] = np.nan

#--- Eccentricity

    if isFullDump:

        specificKineticEnergy = 1/2 * particleData['|v|']**2
        specificGravitationalEnergy = \
            - gravitationalParameter / particleData['r']
        specificEnergy = specificKineticEnergy + specificGravitationalEnergy
        term = 2 * specificEnergy * particleData['|l|']**2 \
             / gravitationalParameter**2

        particleData['e'] = np.sqrt( 1 + term )

        specificKineticEnergy = 1/2 * sinkData['|v|']**2
        specificGravitationalEnergy = \
            - gravitationalParameter / sinkData['r']
        specificEnergy = specificKineticEnergy + specificGravitationalEnergy
        term = 2 * specificEnergy * sinkData['|l|']**2 \
             / gravitationalParameter**2

        sinkData['e'] = np.sqrt( 1 + term )

    else:

        particleData['e'] = np.nan
        sinkData['e']     = np.nan

#--- Calculate radially binned quantities

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
        + ['sigma', 'h', 'H', 'Lx', 'Ly', 'Lz', 'l', 'tilt', 'twist', 'e'] )

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

            radialData['Lx'].iloc[index] = particlesInRadialBin['Lx'].mean()
            radialData['Ly'].iloc[index] = particlesInRadialBin['Ly'].mean()
            radialData['Lz'].iloc[index] = particlesInRadialBin['Lz'].mean()

            radialData['e'].iloc[index] = particlesInRadialBin['e'].mean()


    radialData['l'] = norm(radialData[['Lx','Ly','Lz']].values, axis=1)

    radialData['tilt'] = np.arccos( radialData['Lz'] / radialData['l'] )

    radialData['twist'] = np.arctan2( radialData['Ly'] / radialData['l'],
                                      radialData['Lx'] / radialData['l'] )

    radialData['rho'] = radialData['sigma'] / radialData['H'] / np.sqrt(2*np.pi)

    radialData['psi'] = radialBins * np.sqrt(
        np.gradient( radialData['Lx']/radialData['l'], radialBinWidth )**2 + \
        np.gradient( radialData['Ly']/radialData['l'], radialBinWidth )**2 + \
        np.gradient( radialData['Lz']/radialData['l'], radialBinWidth )**2 )

    return radialData
