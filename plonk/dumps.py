'''
dump.py

Daniel Mentiplay, 2019.
'''

import h5py
import numpy as np

from .itypes import iTypes
from .units import Units

positionIndex = slice(0, 3)
massIndex = 3
smoothingLengthIndex = 4
velocityIndex = slice(6, 9)

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
        iDust = iTypes.iDustSplash

        itype = data[:,-1]

        #--- Position

        position = dict()

        position['gas'] = data[np.where(itype == iGas), positionIndex][0]

        if containsLargeDust:

            positionDust = list()
            for i in range(iDust, iDust + nDustLarge):
                positionDust.append(data[np.where(itype == i)[0], positionIndex])
            position['dust'] = positionDust

        #--- Velocity

        if isFullDump:

            velocity = dict()

            velocity['gas'] = data[np.where(itype == iGas), velocityIndex][0]

            if containsLargeDust:

                velocityDust = list()
                for i in range(iDust, iDust + nDustLarge):
                    velocityDust.append(data[np.where(itype == i)[0], velocityIndex])
                velocity['dust'] = velocityDust

        #--- Smoothing length

        smoothingLength = dict()

        smoothingLength['gas'] = data[np.where(itype == iGas),
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


        #--- Add to object

        self.position = position
        self.smoothingLength = smoothingLength

        if isFullDump:
            self.velocity = velocity

        if containsDust:
            self.dustFrac = dustFrac

class Header:
    '''
    Header class represents the header from a dump file (in the format produced
    by the phantom showheader utility program).

    Arguments:
        filename : e.g. disc_00000.header
    '''

    def __init__(self, filename):

        self.filename = filename

        fileExtension = filename.split('.')[-1]

        if fileExtension  == 'h5':
            self._read_hdf5_file()

        elif fileExtension == 'header':
            self._read_showheader_file()

        else:
            raise ValueError('Cannot read header from file')

    def _read_showheader_file(self):

        with open(self.filename, 'r') as file:

            keys = list()
            values = list()

            firstLine = True

            for line in file:

                if firstLine:

                    if line[0] == 'F':
                        self.dumpType = 'full'
                    elif line[0] == 'S':
                        self.dumpType = 'small'
                    else:
                        raise ValueError('Cannot determine dump type')

                    if 'dust' in line:
                        self.containsDust = True
                    else:
                        self.containsDust = False

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

        parameters = dict(zip(newKeys, newValues))

        for key in parameters:
            if isinstance(parameters[key], list):
                parameters[key] = np.array(parameters[key])

        self.parameters = parameters

        units = Units(parameters['udist'], parameters['umass'], parameters['utime'])
        self.units = units.units

    def _read_hdf5_file(self):

        f = h5py.File(self.filename, 'r')

        headerGroup = f['header']

        parameters = dict()
        for key in headerGroup.keys():
            parameters[key] = headerGroup[key].value

        f.close()

        self.parameters = parameters

        units = Units(parameters['udist'], parameters['umass'], parameters['utime'])
        self.units = units.units

        warning = '''
Warning: hdf5 header read cannot determine if full/small dump or if the dump
contains dust. For now assume full dump and no dust.
        '''
        print(warning)

        self.dumpType = 'full'
        self.containsDust = False

class ParticleTypes:
    '''
    ParticleTypes class represents integer labels for particle types.
    '''

    # pylint: disable=too-many-instance-attributes

    def __init__(self):

        self.iGas  = 1
        self.iSink = 3
        self.iDust = 7

        self.iGasSplash  = 1
        self.iSinkSplash = 3
        self.iDustSplash = 8

        self.iGasLabel  = 'gas'
        self.iSinkLabel = 'sink'
        self.iDustLabel = 'dust'
