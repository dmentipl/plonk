'''
PhantomDump.py

Daniel Mentiplay, 2019.
'''

import os

import h5py
import numpy as np

from .ParticleData import ParticleData
from .SinkData import SinkData
from .units import Units
from .utils import print_warning, print_error

# ---------------------------------------------------------------------------- #

#--- Reading Phantom ACSII files

iGas         = 1
iSink        = 3
iDust        = 7
iDustSplash  = 8

positionIndex        = slice(0, 3)
massIndex            = 3
smoothingLengthIndex = 4
velocityIndex        = slice(6, 9)

# ---------------------------------------------------------------------------- #

class PhantomDump:
    '''
    Phantom dump.
    '''

    def __init__(self, filename):

        self.ParticleData = ParticleData()
        self.SinkData     = SinkData()
        self.Parameters   = None
        self.Units        = None

        self._read_dump(filename)

    #--- Read dump

    def _read_dump(self, filename):
        '''
        Read a Phantom dump file.

        Two file types can be read: ASCII or HDF5.

        1. For HDF5 only one file is required and it must have the extension h5,
        e.g. 'disc_00000.h5'.

        2. For ASCII dumps two files are required:
            (i)  the dump file in ASCII format, e.g. 'disc_00000.ascii'
            (ii) the header file in "showheader" format, e.g.
            'disc_00000.header'

        Arguments:
            filename : e.g. 'disc_00000.h5' for HDF5 or 'disc_00000.ascii' for
            ASCII
        '''

        # TODO: keep working on HDF5 Phantom output

        exists = os.path.isfile(filename)

        if not exists:
            raise FileNotFoundError('Cannot find dump file')

        filePrefix    = filename.split('.')[0]
        fileExtension = filename.split('.')[-1]

        if fileExtension == 'h5':
            dumpFileFormat = 'HDF5'
            print_error('HDF5 dump reader not fully implemented')

        elif fileExtension == 'ascii':
            dumpFileFormat = 'ASCII'

        else:
            raise ValueError('Cannot determine dump file format')

        isFullDump = self._read_header(filePrefix, dumpFileFormat)
        self._read_arrays(filePrefix, dumpFileFormat, isFullDump)

    #--- Read header

    def _read_header(self, filePrefix, dumpFileFormat):

        if dumpFileFormat == 'HDF5':
            fileExtension = 'h5'

        elif dumpFileFormat == 'ASCII':
            fileExtension = 'header'

        headerFileName = filePrefix + '.' + fileExtension

        exists = os.path.isfile(headerFileName)
        if not exists:
            raise FileNotFoundError('Cannot find header file: ' +
                                    headerFileName)

        if dumpFileFormat == 'HDF5':
            header = _read_header_from_hdf5(headerFileName)

        elif dumpFileFormat == 'ASCII':
            header, isFullDump = _read_header_from_showheader(headerFileName)

        self.Parameters = header

        self.Units = Units(header['udist'], header['umass'],
                           header['utime']).units

        return isFullDump

    #--- Read arrays

    def _read_arrays(self, filePrefix, dumpFileFormat, isFullDump):

        if dumpFileFormat == 'HDF5':
            fileExtension = 'h5'

        elif dumpFileFormat == 'ASCII':
            fileExtension = 'ascii'

        dumpFileName = filePrefix + '.' + fileExtension

        exists = os.path.isfile(dumpFileName)
        if not exists:
            raise FileNotFoundError('Cannot find dump file: ' + dumpFileName)

        nDustSmall = self.Parameters['ndustsmall']
        nDustLarge = self.Parameters['ndustlarge']

        containsSmallDust = bool(nDustSmall > 0)
        containsLargeDust = bool(nDustLarge > 0)

        nDustTypes = nDustSmall + nDustLarge
        containsDust = bool(containsSmallDust or containsLargeDust)

        nSinks = self.Parameters['nptmass']
        containsSinks = bool(nSinks > 0)

        if dumpFileFormat == 'HDF5':

            f = h5py.File(dumpFileName, 'r')
            arrays = f['arrays']

            self.ParticleData.particleType     = arrays['itype'].value
            self.ParticleData.position         = arrays['xyzh_label'].value[:, 0:3]
            self.ParticleData.smoothingLength  = arrays['xyzh_label'].value[:, 3]

            if isFullDump:
                self.ParticleData.velocity     = arrays['vxyzu_label'].value[:, 0:3]
            else:
                self.ParticleData.velocity     = None

            if containsSinks:
                self.SinkData.mass            = arrays['xyzmh_ptmass_label'].value[:, 3]
                self.SinkData.position        = arrays['xyzmh_ptmass_label'].value[:, 0:3]
                self.SinkData.accretionRadius = arrays['xyzmh_ptmass_label'].value[:, 4]
                self.SinkData.velocity        = arrays['vxyz_ptmass_label'].value
            else:
                self.SinkData = None

            # TODO: add dustfrac
            if containsDust:
                print_warning('HDF5 dump reader cannot read dustfrac yet')

            f.close()

        elif dumpFileFormat == 'ASCII':

            data = np.loadtxt(dumpFileName)

            itype = data[:, -1]
            particleIndicies = np.where(itype != iSink)
            sinkIndicies     = np.where(itype == iSink)

            self.ParticleData.particleType    = data[particleIndicies, -1][0]
            self.ParticleData.position        = data[particleIndicies, positionIndex][0]
            self.ParticleData.particleMass    = data[particleIndicies, massIndex][0]
            self.ParticleData.smoothingLength = data[particleIndicies, smoothingLengthIndex][0]

            if isFullDump:
                self.ParticleData.velocity = data[particleIndicies, velocityIndex][0]
            else:
                self.ParticleData.velocity = None

            if containsDust:

                if isFullDump:
                    dustFracIndexStart = 9
                else:
                    dustFracIndexStart = 6

                dustFracIndex = slice(dustFracIndexStart,
                                      dustFracIndexStart + nDustTypes)

                self.ParticleData.dustFrac = data[particleIndicies, dustFracIndex][0]

            else:
                self.ParticleData.dustFrac = None

            if containsSinks:

                self.SinkData.position        = data[sinkIndicies, positionIndex][0]
                self.SinkData.mass            = data[sinkIndicies, massIndex][0]
                self.SinkData.accretionRadius = data[sinkIndicies, smoothingLengthIndex][0]

                if isFullDump:
                    self.SinkData.velocity = data[sinkIndicies, velocityIndex][0]
                else:
                    self.SinkData.velocity = None

            else:

                self.SinkData = None

#--- Functions

def _read_header_from_hdf5(headerFileName):

    f = h5py.File(headerFileName, 'r')

    headerGroup = f['header']

    header = dict()
    for key in headerGroup.keys():
        header[key] = headerGroup[key].value

    f.close()

    return header

def _read_header_from_showheader(headerFileName):

    with open(headerFileName, 'r') as file:

        keys = list()
        values = list()

        firstLine = True

        for line in file:

            if firstLine:

                if line[0] == 'F':
                    isFullDump = True
                elif line[0] == 'S':
                    isFullDump = False
                else:
                    raise ValueError('Cannot determine dump type')

                firstLine = False
                continue

            keys.append(line.rstrip('\n').split()[0])
            value = line.rstrip('\n').split()[1]

            if '.' in value:
                value = float(value)
            else:
                value = int(value)

            values.append(value)

    newKeys = list()
    newValues = list()

    prevKey = None
    firstMultipleKey = True
    idxj = 0

    for idxi, key in enumerate(keys):

        if key == prevKey:

            if firstMultipleKey:
                newValues[idxj-1] = list()
                newValues[idxj-1].append(values[idxi-1])
                newValues[idxj-1].append(values[idxi])
                firstMultipleKey = False
            else:
                newValues[idxj-1].append(values[idxi])

        else:

            firstMultipleKey = True
            idxj += 1
            newKeys.append(key)
            newValues.append(values[idxi])

        prevKey = key

    header = dict(zip(newKeys, newValues))

    for key in header:
        if isinstance(header[key], list):
            header[key] = np.array(header[key])

    return header, isFullDump
