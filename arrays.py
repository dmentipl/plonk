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
        data         : 2d numpy array of particle data
        dumpType     : 'full' or 'small'

    Optional:
        nDustSmall   : number of small dust types
        nDustLarge   : number of large dust types
    '''

    def __init__(self, data, dumpType, nDustSmall=None, nDustLarge=None):

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

        if nDustSmall is not None:
            if isinstance(nDustSmall, int):
                if nDustSmall > 0:
                    containsSmallDust = True
                elif nDustSmall == 0:
                    containsSmallDust = False
                else:
                    raise ValueError('nDustSmall must be > 0')
            else:
                raise ValueError('nDustSmall must be integer')
        else:
            containsSmallDust = False
            nDustSmall = 0

        if nDustLarge is not None:
            if isinstance(nDustLarge, int):
                if nDustLarge > 0:
                    containsLargeDust = True
                elif nDustLarge == 0:
                    containsLargeDust = False
                else:
                    raise ValueError('nDustLarge must be > 0')
            else:
                raise ValueError('nDustLarge must be integer')
        else:
            containsLargeDust = False
            nDustLarge = 0

        containsDust = bool(containsSmallDust or containsLargeDust)
        nDustTypes = nDustSmall + nDustLarge

        iGas  = iTypes.iGasSplash
        iSink = iTypes.iSinkSplash
        iDust = iTypes.iDustSplash

        itype = data[:,-1]

        #--- Position

        position = dict()

        position['gas'] = data[np.where(itype == iGas), positionIndex][0]
        position['sink'] = data[np.where(itype == iSink), positionIndex][0]

        if containsLargeDust:

            positionDust = list()
            for i in range(iDust, iDust + nDustLarge):
                positionDust.append(data[np.where(itype == i)[0], positionIndex])
            position['dust'] = positionDust

        #--- Velocity

        if isFullDump:

            velocity = dict()

            velocity['gas'] = data[np.where(itype == iGas), velocityIndex][0]
            velocity['sink'] = data[np.where(itype == iSink), velocityIndex][0]

            if containsLargeDust:

                velocityDust = list()
                for i in range(iDust, iDust + nDustLarge):
                    velocityDust.append(data[np.where(itype == i)[0], velocityIndex])
                velocity['dust'] = velocityDust

        #--- Smoothing length

        smoothingLength = dict()

        smoothingLength['gas'] = data[np.where(itype == iGas),
                                      smoothingLengthIndex][0]
        smoothingLength['sink'] = data[np.where(itype == iSink),
                                       smoothingLengthIndex][0]

        if containsLargeDust:

            smoothingLengthDust = list()
            for i in range(iDust, iDust + nDustLarge):
                smoothingLengthDust.append(data[np.where(itype == i)[0],
                                                smoothingLengthIndex])
            smoothingLength['dust'] = smoothingLengthDust

        #--- Dust fraction

        if containsDust:

            if isFullDump:
                dustFracIndexStart = 9
            else:
                dustFracIndexStart = 6

            dustFracIndex = slice(dustFracIndexStart, dustFracIndexStart + nDustTypes)
            dustFrac = data[np.where(itype == iGas), dustFracIndex][0]

        #--- Sink masses

        massSink = list()
        for i in np.where(itype == iSink)[0]:
            massSink.append(data[i, massIndex])

        #--- Add to object

        self.position = position
        self.smoothingLength = smoothingLength

        if isFullDump:
            self.velocity = velocity

        if containsDust:
            self.dustFrac = dustFrac

        if massSink:
            self.massSink = massSink
