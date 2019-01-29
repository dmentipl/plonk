'''
header.py

Daniel Mentiplay, 2019.
'''

import h5py
import numpy as np

from .units import Units

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
