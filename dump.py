'''
dump.py

Daniel Mentiplay, 2019.
'''

import numpy as np

from constants import constants
from itypes import iTypes

class Dump:
    '''
    Dump class represents a dump file (ascii only).
    '''

    # pylint: disable=too-many-instance-attributes

    def __init__(self, filename):

        self.filename = filename
        self._read_file()
        self._get_from_header()
        self._get_from_data()

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

        self.header = header

        #--- Data

        self.data = np.loadtxt(self.filename)

    def _get_from_header(self):

        header = self.header

        #--- Column labels

        columnLabels = header[-1].split('  ')
        columnLabels = [col for col in columnLabels if col != '']
        columnLabels = [col[1:] if col[0] == ' ' else col for col in columnLabels]
        self.columnLabels = columnLabels

        self.nColumns = len(self.columnLabels)

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

        distUnit = float(self.header[6].split()[0])
        distUnitLabel = self.header[7].split()[0][1:-1]

        if distUnitLabel in ['au']:
            distUnit *= constants.au
        elif distUnitLabel in ['cm']:
            pass
        else:
            raise ValueError('Cannot determine dist unit')

        units['dist'] = distUnit

        massUnit = float(self.header[6].split()[3])
        massUnitLabel = self.header[7].split()[3][1:-1]

        if massUnitLabel in ['g']:
            pass
        else:
            raise ValueError('Cannot determine mass unit')

        units['mass'] = massUnit

        velocityUnit = float(self.header[6].split()[6])
        velocityUnitLabel = self.header[7].split()[6][1:-1]

        if velocityUnitLabel in ['cm/s']:
            pass
        else:
            raise ValueError('Cannot determine velocity unit')

        units['velocity'] = velocityUnit

        self.units = units

    def _get_from_data(self):

        data = self.data

        iGas = iTypes.iGas
        iSink = iTypes.iSink
        iDust = iTypes.iDust
        maxDustTypes = iTypes.maxDustTypes

        #--- Number of particles

        self.nPart = data.shape[0]

        itype = data[:,-1]

        self.nGas = len(np.where(itype == iGas)[0])
        self.nSink = len(np.where(itype == iSink)[0])

        self.nDust = list()
        for i in range(iDust, maxDustTypes+1):
            l = len(np.where(itype == i)[0])
            if l > 0:
                self.nDust.append(l)

        self.nDustTypes = len(self.nDust)
