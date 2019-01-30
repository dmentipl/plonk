'''
dumps.py

Daniel Mentiplay, 2019.
'''

import os

import h5py
import numpy as np

from .particles import Gas, Dust
from .parameters import Parameters
from .sinks import Sinks
from .units import Units

# ---------------------------------------------------------------------------- #

#--- Reading Phantom ACSII files

iGas  = 1
iSink = 3
iDust = 8

positionIndex = slice(0, 3)
massIndex = 3
smoothingLengthIndex = 4
velocityIndex = slice(6, 9)

# ---------------------------------------------------------------------------- #

class Dump:
    '''
    Phantom dump.
    '''

    # pylint: disable=too-many-instance-attributes

    def __init__(self):

        self.gas = Gas()
        self.dust = list()
        self.sinks = Sinks()
        self.parameters = Parameters()
        self.time = 0.0
        self.filename = ''
        self.fullDump = True

#--- Read dump

    def read_dump(self, filename):
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

        exists = os.path.isfile(filename)

        if exists:
            self.filename = filename
        else:
            raise FileNotFoundError

        filePrefix    = filename.split('.')[0]
        fileExtension = filename.split('.')[-1]

        if fileExtension == 'h5':
            dumpFileFormat = 'HDF5'

        elif fileExtension == 'ascii':
            dumpFileFormat = 'ASCII'

        self._read_header(filePrefix, dumpFileFormat)
        self._read_arrays(filePrefix, dumpFileFormat)

#--- Read header

    def _read_header(self, filePrefix, dumpFileFormat):

        if dumpFileFormat == 'HDF5':
            fileExtension = 'h5'

        elif dumpFileFormat == 'ASCII':
            fileExtension = 'header'

        headerFileName = filePrefix + '.' + fileExtension

        exists = os.path.isfile(headerFileName)
        if not exists:
            raise FileNotFoundError

        if dumpFileFormat == 'HDF5':
            header = _read_header_from_hdf5(headerFileName)

        elif dumpFileFormat == 'ASCII':
            header = self._read_header_from_showheader(headerFileName)

        self.parameters.dust['nDustSmall']   = header.get('ndustsmall')
        self.parameters.dust['nDustLarge']   = header.get('ndustlarge')
        self.parameters.dust['grainSize']    = header.get('grainsize')
        self.parameters.dust['grainDens']    = header.get('graindens')

        self.parameters.eos['ieos']          = header.get('ieos')
        self.parameters.eos['gamma']         = header.get('gamma')
        self.parameters.eos['polyk']         = header.get('polyk')
        self.parameters.eos['qfacdisc']      = header.get('qfacdisc')

        self.parameters.numerical['hfact']   = header.get('hfact')
        self.parameters.numerical['tolh']    = header.get('tolh')
        self.parameters.numerical['C_cour']  = header.get('C_cour')
        self.parameters.numerical['C_force'] = header.get('C_force')
        self.parameters.numerical['alpha']   = header.get('alpha')

        self.parameters.sinks['nSinks']      = header.get('nptmass')

        self.parameters.units = Units(header['udist'], header['umass'],
                                      header['utime']).units

    def _read_header_from_showheader(self, headerFileName):

        with open(headerFileName, 'r') as file:

            keys = list()
            values = list()

            firstLine = True

            for line in file:

                if firstLine:

                    if line[0] == 'F':
                        self.fullDump = True
                    elif line[0] == 'S':
                        self.fullDump = False
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

        return header

#--- Read arrays

    def _read_arrays(self, filePrefix, dumpFileFormat):

        if dumpFileFormat == 'HDF5':
            fileExtension = 'h5'

        elif dumpFileFormat == 'ASCII':
            fileExtension = 'ascii'

        dumpFileName = filePrefix + '.' + fileExtension

        exists = os.path.isfile(dumpFileName)
        if not exists:
            raise FileNotFoundError

        if dumpFileFormat == 'HDF5':

            f = h5py.File(dumpFileName, 'r')

            arraysGroup = f['arrays']

            xyzh = arraysGroup['xyzh_label'].value
            vxyzu = arraysGroup['vxyzu_label'].value
            itype = arraysGroup['itype'].value

            f.close()

            data = np.stack( ( xyzh[:, 0], xyzh[:, 1], xyzh[:, 2],
                               np.ones_like(itype), xyzh[:, 3],
                               np.ones_like(itype), vxyzu[:, 0], vxyzu[:, 1],
                               vxyzu[:, 2], itype ) )

        elif dumpFileFormat == 'ASCII':

            data = np.loadtxt(dumpFileName)

        self._put_arrays_into_objects(data)


    def _put_arrays_into_objects(self, data):

        if isinstance(data, np.ndarray):
            if data.ndim != 2:
                raise ValueError('data must be 2d numpy array')
        else:
            raise ValueError('data must be 2d numpy array')

        nDustSmall = self.parameters.dust['nDustSmall']
        nDustLarge = self.parameters.dust['nDustLarge']

        containsSmallDust = bool(nDustSmall > 0)
        containsLargeDust = bool(nDustLarge > 0)

        nDustTypes = nDustSmall + nDustLarge
        containsDust = bool(containsSmallDust or containsLargeDust)

        itype = data[:,-1]

        containsSinks = self.parameters.sinks['nSinks']

        #--- Gas

        gas = Gas()

        gas.number = len(data[np.where(itype == iGas), massIndex][0])

        gas.mass = data[np.where(itype == iGas), massIndex][0][0]

        gas.position = data[np.where(itype == iGas), positionIndex][0]

        gas.smoothingLength = data[np.where(itype == iGas),
                                   smoothingLengthIndex][0]

        if self.fullDump:
            gas.velocity = data[np.where(itype == iGas), velocityIndex][0]

        if containsDust:

            if self.fullDump:
                dustFracIndexStart = 9
            else:
                dustFracIndexStart = 6
            dustFracIndex = slice(dustFracIndexStart,
                                  dustFracIndexStart + nDustTypes)

            gas.dustFrac = data[np.where(itype == iGas), dustFracIndex][0]

        self.gas = gas

        #--- Dust

        if containsLargeDust:

            dust = list()

            for i in range(nDustLarge):

                dust.append(Dust())

                itype_ = i + iDust

                dust[i].number = len(data[np.where(itype == iDust),
                                          massIndex][0])

                dust[i].mass = data[np.where(itype == iGas), massIndex][0][0]

                dust[i].position = data[np.where(itype == itype_)[0],
                                        positionIndex]

                dust[i].smoothingLength = data[np.where(itype == itype_)[0],
                                               smoothingLengthIndex]

                if self.fullDump:
                    dust[i].velocity = data[np.where(itype == itype_)[0],
                                            velocityIndex]

            self.dust = dust

        #--- Sinks

        if containsSinks:

            sinks = Sinks()

            sinks.number = len(data[np.where(itype == iSink), massIndex][0])

            sinks.mass = data[np.where(itype == iSink), massIndex][0]

            sinks.position = data[np.where(itype == iSink),
                                  positionIndex][0]

            sinks.smoothingLength = data[np.where(itype == iSink),
                                         smoothingLengthIndex][0]

            if self.fullDump:

                sinks.velocity = data[np.where(itype == iSink),
                                      velocityIndex][0]

            self.sinks = sinks

#--- Functions

def _read_header_from_hdf5(headerFileName):

    f = h5py.File(headerFileName, 'r')

    headerGroup = f['header']

    header = dict()
    for key in headerGroup.keys():
        header[key] = headerGroup[key].value

    f.close()

    warning = '''
WARNING: For now HDF5 header read cannot determine if full/small dump
     For now assume full dump.
    '''
    print(warning)

    return header
