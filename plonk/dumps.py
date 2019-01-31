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
from .utils import print_warning, print_error

# ---------------------------------------------------------------------------- #

#--- Reading Phantom ACSII files

iGas         = 1
iSink        = 3
iDustPhantom = 7
iDustSplash  = 8

positionIndex        = slice(0, 3)
massIndex            = 3
smoothingLengthIndex = 4
velocityIndex        = slice(6, 9)

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

        # TODO: keep working on HDF5 Phantom output

        exists = os.path.isfile(filename)

        if exists:
            self.filename = filename
        else:
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
            raise FileNotFoundError('Cannot find header file: ' +
                                    headerFileName)

        if dumpFileFormat == 'HDF5':
            header = _read_header_from_hdf5(headerFileName)

        elif dumpFileFormat == 'ASCII':
            header = self._read_header_from_showheader(headerFileName)

        self.parameters.particles['npartoftype'] = header.get('npartoftype')
        self.parameters.particles['massoftype']  = header.get('massoftype')

        self.parameters.dust['nDustSmall']       = header.get('ndustsmall')
        self.parameters.dust['nDustLarge']       = header.get('ndustlarge')
        self.parameters.dust['grainSize']        = header.get('grainsize')
        self.parameters.dust['grainDens']        = header.get('graindens')

        self.parameters.eos['ieos']              = header.get('ieos')
        self.parameters.eos['gamma']             = header.get('gamma')
        self.parameters.eos['polyk']             = header.get('polyk')
        self.parameters.eos['qfacdisc']          = header.get('qfacdisc')

        self.parameters.numerical['hfact']       = header.get('hfact')
        self.parameters.numerical['tolh']        = header.get('tolh')
        self.parameters.numerical['C_cour']      = header.get('C_cour')
        self.parameters.numerical['C_force']     = header.get('C_force')
        self.parameters.numerical['alpha']       = header.get('alpha')

        self.parameters.sinks['nSinks']          = header.get('nptmass')

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
            raise FileNotFoundError('Cannot find dump file: ' + dumpFileName)

        nDustSmall = self.parameters.dust['nDustSmall']
        nDustLarge = self.parameters.dust['nDustLarge']

        containsSmallDust = bool(nDustSmall > 0)
        containsLargeDust = bool(nDustLarge > 0)

        nDustTypes = nDustSmall + nDustLarge
        containsDust = bool(containsSmallDust or containsLargeDust)

        nSinks = self.parameters.sinks['nSinks']
        containsSinks = bool(nSinks > 0)

        if dumpFileFormat == 'HDF5':

            f = h5py.File(dumpFileName, 'r')
            arrays = f['arrays']

            particleType            = arrays['itype'].value

            particlePosition        = arrays['xyzh_label'].value[:, 0:3]
            particleSmoothingLength = arrays['xyzh_label'].value[:, 3]

            if self.fullDump:
                particleVelocity    = arrays['vxyzu_label'].value[:, 0:3]

            if containsSinks:

                sinkPosition        = arrays['xyzmh_ptmass_label'].value[:, 0:3]
                sinkMass            = arrays['xyzmh_ptmass_label'].value[:, 3]
                sinkAccretionRadius = arrays['xyzmh_ptmass_label'].value[:, 4]

                sinkVelocity        = arrays['vxyz_ptmass_label'].value

            # TODO: add dustfrac
            if containsDust:
                print_warning('HDF5 dump reader cannot read dustfrac yet')

            f.close()

        elif dumpFileFormat == 'ASCII':

            data = np.loadtxt(dumpFileName)

            particleType            = data[:, -1]

            particlePosition        = data[:, positionIndex]
            particleSmoothingLength = data[:, smoothingLengthIndex]

            if self.fullDump:
                particleVelocity    = data[:, velocityIndex]

            if containsDust:

                if self.fullDump:
                    dustFracIndexStart = 9
                else:
                    dustFracIndexStart = 6

                dustFracIndex = slice(dustFracIndexStart,
                                      dustFracIndexStart + nDustTypes)

                dustFrac            = data[:, dustFracIndex]

            if containsSinks:

                sinkPosition        = data[np.where(particleType == iSink),
                                           positionIndex][0]

                sinkMass            = data[np.where(particleType == iSink),
                                           massIndex][0]

                sinkAccretionRadius = data[np.where(particleType == iSink),
                                           smoothingLengthIndex][0]

                if self.fullDump:
                    sinkVelocity    = data[np.where(particleType == iSink),
                                           velocityIndex][0]

        #--- Gas

        gas = Gas()

        gas.number = self.parameters.particles['npartoftype'][0]
        gas.mass   = self.parameters.particles['massoftype'][0]

        gas.position        = particlePosition[np.where(particleType == iGas)]
        gas.smoothingLength = particleSmoothingLength[np.where(particleType == iGas)]

        if self.fullDump:
            gas.velocity = particleVelocity[np.where(particleType == iGas)]
        else:
            gas.velocity = None

        if containsDust:
            gas.dustFrac = dustFrac[np.where(particleType == iGas)]
        else:
            gas.dustfrac = None

        self.gas = gas

        #--- Dust

        if containsLargeDust:

            dust = list()

            for idx in range(nDustLarge):

                dust.append(Dust())

                itypePhantom = idx + iDustPhantom - 1
                itypeSplash  = idx + iDustSplash

                dust[idx].number = \
                    self.parameters.particles['npartoftype'][itypePhantom]

                dust[idx].mass = \
                    self.parameters.particles['massoftype'][itypePhantom]

                dust[idx].position = particlePosition[
                    np.where(particleType == itypeSplash)]

                dust[idx].smoothingLength = particleSmoothingLength[
                    np.where(particleType == itypeSplash)]

                if self.fullDump:
                    dust[idx].velocity = particleVelocity[
                        np.where(particleType == itypeSplash)]

            self.dust = dust

        else:

            self.dust = None

        #--- Sinks

        if containsSinks:

            sinks = Sinks()

            sinks.number          = nSinks
            sinks.mass            = sinkMass
            sinks.position        = sinkPosition
            sinks.accretionRadius = sinkAccretionRadius
            if self.fullDump:
                sinks.velocity    = sinkVelocity

            self.sinks = sinks

        else:

            self.sinks = None

#--- Functions

def _read_header_from_hdf5(headerFileName):

    f = h5py.File(headerFileName, 'r')

    headerGroup = f['header']

    header = dict()
    for key in headerGroup.keys():
        header[key] = headerGroup[key].value

    f.close()

    return header
