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

        self.parameters = header.header

        self.dumpType = header.dumpType
        self.containsDust = header.containsDust

        if 'ndustsmall' in header.header:
            nDustSmall = header.header['ndustsmall']
        else:
            nDustSmall = 0

        if 'ndustlarge' in header.header:
            nDustLarge = header.header['ndustlarge']
        else:
            nDustLarge = 0

        #--- Data

        dumpFileName = self.filePrefix + '.ascii'
        data = np.loadtxt(dumpFileName)
        arrays = Arrays(data, self.dumpType, nDustSmall, nDustLarge)
        self.arrays = arrays
