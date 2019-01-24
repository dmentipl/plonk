'''
arrays.py

Daniel Mentiplay, 2019.
'''

import numpy as np

from itypes import iTypes

positionIndex = slice(0, 3)
massIndex = 3
smoothingLengthIndex = 4
velocityIndex = slice(6, 9)

class Arrays:
    '''
    Arrays class represents the particle arrays from Phantom output.

    Arguments:
        data     : 2d numpy array of particle data
        dumpType : 'full' or 'small'
    '''

    def __init__(self, data, dumpType):

        if isinstance(data, np.ndarray):
            if data.ndim != 2:
                raise ValueError('data must be 2d numpy array')
        else:
            raise ValueError('data must be 2d numpy array')

        if dumpType in ['F', 'f', 'Full', 'full']:
            isFullDump = True
        elif dumpType in ['s', 'S', 'Small', 'small']:
            isFullDump = False
        else:
            raise ValueError('Cannot determine dump type')

        iGas = iTypes.iGas
        iSink = iTypes.iSink
        iDust = iTypes.iDust
        maxDustTypes = iTypes.maxDustTypes

        #--- Number of particles

        nParticles = dict()

        itype = data[:,-1]

        nParticles['gas'] = len(np.where(itype == iGas)[0])
        nParticles['sink'] = len(np.where(itype == iSink)[0])

        nDust = list()
        for i in range(iDust, maxDustTypes+1):
            l = len(np.where(itype == i)[0])
            if l > 0:
                nDust.append(l)
        nParticles['dust'] = nDust
        nDustTypes = len(nDust)

        #--- Mass of particles

        massParticles = dict()

        massParticles['gas'] = data[np.where(itype == iGas)[0][0], massIndex]

        massSink = list()
        for i in np.where(itype == iSink)[0]:
            massSink.append(data[i, massIndex])
        massParticles['sink'] = massSink

        massDust = list()
        for i in range(iDust, iDust + nDustTypes):
            massDust.append(data[np.where(itype == i)[0][0], massIndex])
        massParticles['dust'] = massDust

        #--- Position

        position = dict()

        position['gas'] = data[np.where(itype == iGas), positionIndex][0]
        position['sink'] = data[np.where(itype == iSink), positionIndex][0]

        positionDust = list()
        for i in range(iDust, iDust + nDustTypes):
            positionDust.append(data[np.where(itype == i)[0], positionIndex])
        position['dust'] = positionDust

        #--- Velocity

        if isFullDump:

            velocity = dict()

            velocity['gas'] = data[np.where(itype == iGas), velocityIndex][0]
            velocity['sink'] = data[np.where(itype == iSink), velocityIndex][0]

            velocityDust = list()
            for i in range(iDust, iDust + nDustTypes):
                velocityDust.append(data[np.where(itype == i)[0], velocityIndex])
            velocity['dust'] = velocityDust

        #--- Smoothing length

        smoothingLength = dict()

        smoothingLength['gas'] = data[np.where(itype == iGas),
                                      smoothingLengthIndex][0]
        smoothingLength['sink'] = data[np.where(itype == iSink),
                                       smoothingLengthIndex][0]

        smoothingLengthDust = list()
        for i in range(iDust, iDust + nDustTypes):
            smoothingLengthDust.append(data[np.where(itype == i)[0],
                                            smoothingLengthIndex])
        smoothingLength['dust'] = smoothingLengthDust


        #--- Dust fraction

        if isFullDump:
            dustFracIndexStart = 9
        else:
            dustFracIndexStart = 6

        dustFracIndex = slice(dustFracIndexStart, dustFracIndexStart + nDustTypes)

        #--- Add to class

        self.nParticles = nParticles
        self.massParticles = massParticles
        self.position = position
        self.smoothingLength = smoothingLength
        if isFullDump:
            self.velocity = velocity
        self.dustFrac = data[np.where(itype == iGas), dustFracIndex][0]
