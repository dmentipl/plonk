'''
dump.py

Daniel Mentiplay, 2019.
'''

import numpy as np

from constants import constants
from itypes import iTypes

positionIndex = slice(0, 3)
massIndex = 3
smoothingLengthIndex = 4
velocityIndex = slice(6, 9)

class Dump:
    '''
    Dump class represents a dump file (ascii only).
    '''

    def __init__(self, filename):

        self.filename = filename
        header, data = self._read_file()
        self._get_from_header(header)
        self._get_from_data(data)

    def _read_file(self):

        #--- Header

        with open(self.filename, 'r') as file:
            header = list()
            for line in file:
                line = line.rstrip('\n')
                if '#' in line:
                    if line == '#':
                        pass
                    else:
                        header.append(line[2:])
                else:
                    break

        #--- Dump type

        if 'v_x' in header[8].split():
            self.dumpType = 'full'
        else:
            self.dumpType = 'small'

        #--- Data

        data = np.loadtxt(self.filename)

        #--- Return

        return header, data

    def _get_from_header(self, header):

        #--- Column labels

        columnLabels = header[-1].split('  ')
        columnLabels = [col for col in columnLabels if col != '']
        columnLabels = [col[1:] if col[0] == ' ' else col for col in columnLabels]

        nColumns = len(columnLabels)

        #--- Units

        units = dict()

        self.time = float(header[2].split()[0])

        timeUnit = float(header[2].split()[1])
        timeUnitLabel = header[1].split('(')[1].split()[0].split(')')[0]

        if timeUnitLabel in ['yr', 'yrs', 'years', 'year']:
            timeUnit *= constants.year
        else:
            raise ValueError('Cannot determine time unit')

        units['time'] = timeUnit

        distUnit = float(header[6].split()[0])
        distUnitLabel = header[7].split()[0][1:-1]

        if distUnitLabel in ['au']:
            distUnit *= constants.au
        elif distUnitLabel in ['cm']:
            pass
        else:
            raise ValueError('Cannot determine dist unit')

        units['dist'] = distUnit

        massUnit = float(header[6].split()[3])
        massUnitLabel = header[7].split()[3][1:-1]

        if massUnitLabel in ['g']:
            pass
        else:
            raise ValueError('Cannot determine mass unit')

        units['mass'] = massUnit

        self.units = units

    def _get_from_data(self, data):

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

        self.nParticles = nParticles

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

        self.massParticles = massParticles

        #--- Position

        position = dict()

        position['gas'] = data[np.where(itype == iGas), positionIndex][0]
        position['sink'] = data[np.where(itype == iSink), positionIndex][0]

        positionDust = list()
        for i in range(iDust, iDust + nDustTypes):
            positionDust.append(data[np.where(itype == i)[0], positionIndex])
        position['dust'] = positionDust

        self.position = position

        #--- Velocity

        if self.dumpType == 'full':

            velocity = dict()

            velocity['gas'] = data[np.where(itype == iGas), velocityIndex][0]
            velocity['sink'] = data[np.where(itype == iSink), velocityIndex][0]

            velocityDust = list()
            for i in range(iDust, iDust + nDustTypes):
                velocityDust.append(data[np.where(itype == i)[0], velocityIndex])
            velocity['dust'] = velocityDust

            self.velocity = velocity

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

        self.smoothingLength = smoothingLength

        #--- Dust fraction

        if self.dumpType == 'full':
            dustFracIndexStart = 9
        else:
            dustFracIndexStart = 6

        dustFracIndex = slice(dustFracIndexStart, dustFracIndexStart + nDustTypes)
        self.dustFrac = data[np.where(itype == iGas), dustFracIndex][0]
