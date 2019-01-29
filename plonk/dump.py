'''
dump.py

Daniel Mentiplay, 2019.
'''

import h5py
import numpy as np

from .arrays import Arrays
from .header import Header

class Dump:
    '''
    Dump class represents output from a Phantom calculation at one time.

    Two file types can be read: ASCII or HDF5.

    1. For HDF5 only one file is required and it must have the extension h5,
    e.g. 'disc_00000.h5'.

    2. For ASCII dumps two files are required:
        (i)  the dump file in ASCII format, e.g. 'disc_00000.ascii'
        (ii) the header file in "showheader" format, e.g. 'disc_00000.header'

    Arguments:
        filename : e.g. 'disc_00000.h5' for HDF5 or 'disc_00000.ascii' for ASCII
    '''

    def __init__(self, filename):

        self.filename = filename
        filePrefix    = filename.split('.')[0]
        fileExtension = filename.split('.')[-1]

        if fileExtension == 'h5':
            dumpFileFormat = 'HDF5'
        elif fileExtension == 'ascii':
            dumpFileFormat = 'ASCII'

        #--- Header

        headerFileName = filePrefix + '.header'
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

        if dumpFileFormat == 'HDF5':
            data = self._read_hdf5_dump()

        elif dumpFileFormat == 'ASCII':
            data = np.loadtxt(filename)

        arrays = Arrays(data, self.dumpType, nDustSmall, nDustLarge)
        self.arrays = arrays


    def _read_hdf5_dump(self):

        f = h5py.File(self.filename, 'r')

        arraysGroup = f['arrays']

        xyzh = arraysGroup['xyzh_label'].value
        vxyzu = arraysGroup['vxyzu_label'].value
        itype = arraysGroup['itype'].value

        f.close()

        data = np.stack( ( xyzh[:, 0], xyzh[:, 1], xyzh[:, 2],
                           np.ones_like(itype), xyzh[:, 3],
                           np.ones_like(itype),  vxyzu[:, 0], vxyzu[:, 1],
                           vxyzu[:, 2], itype ) )

        return data
