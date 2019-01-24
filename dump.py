'''
dump.py

Daniel Mentiplay, 2019.
'''

import numpy as np

from arrays import Arrays
from header import Header

class Dump:
    '''
    Dump class represents output from a Phantom calculation at one time.

    Two files are required:
        (i)  the dump file in ascii format, e.g. 'disc_00000.ascii'
        (ii) the header file in showheader format, e.g. 'disc_00000.header'

    Arguments:
        filePrefix : e.g. 'disc_00000' using the example above
    '''

    def __init__(self, filePrefix):

        self.filePrefix = filePrefix

        #--- Header

        headerFileName = self.filePrefix + '.header'
        header = Header(headerFileName)

        self.parameters = header.parameters
        self.units = header.units
        self.dumpType = header.dumpType
        self.containsDust = header.containsDust
        self.time = self.parameters['time']

        nDustSmall = 0
        nDustLarge = 0
        if 'ndustsmall' in self.parameters:
            nDustSmall = self.parameters['ndustsmall']
        if 'ndustlarge' in self.parameters:
            nDustLarge = self.parameters['ndustlarge']

        #--- Arrays

        dumpFileName = self.filePrefix + '.ascii'
        data = np.loadtxt(dumpFileName)
        arrays = Arrays(data, self.dumpType, nDustSmall, nDustLarge)
        self.arrays = arrays
