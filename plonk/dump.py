'''
dump.py

Daniel Mentiplay, 2019.
'''

import os

import h5py
import numpy as np
import pandas as pd

from .particles import density_from_smoothing_length, I_GAS, I_DUST
from .utils import print_warning
from .units import Units

#--- Reading Phantom-Splash ACSII files

I_SINK         = 3
I_DUST_SPLASH  = 8

class Dump:
    '''
    Phantom dump.
    '''

    def __init__(self, filename):

        self.particles  = None
        self.sinks      = None
        self.parameters = None
        self.units      = None

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

        file_prefix    = filename.split('.')[0]
        file_extension = filename.split('.')[-1]

        if file_extension == 'h5':
            self._read_hdf5(file_prefix)

        elif file_extension == 'ascii':
            self._read_ascii(file_prefix)

        else:
            raise ValueError('Cannot determine dump file format')

    #--- Read hdf5 dump

    def _read_hdf5(self, file_prefix):

        #--- Open file

        dump_file_name = file_prefix + '.h5'

        exists = os.path.isfile(dump_file_name)

        if not exists:
            raise FileNotFoundError('Cannot find dump file: ' + dump_file_name)

        f = h5py.File(dump_file_name, 'r')

        #--- Header

        header = f['header']

        parameters = dict()
        for key in header.keys():
            parameters[key] = header[key].value

        self.parameters = parameters

        self.units = Units(parameters['udist'],
                           parameters['umass'],
                           parameters['utime']).units

        n_dust_large = self.parameters['ndustlarge']
        n_dust_types = self.parameters['ndustsmall'] + self.parameters['ndustlarge']

        n_sinks = self.parameters['nptmass']
        contains_sinks = bool(n_sinks > 0)

        #--- Particles

        particles = f['particles']

        is_full_dump = bool( 'vxyz' in particles )

        self.particles = \
            pd.DataFrame( particles['xyz'][:, 0], columns=['x'] )
        self.particles['y'] = particles['xyz'][:, 1]
        self.particles['z'] = particles['xyz'][:, 2]

        self.particles['m'] = np.nan

        self.particles['h'] = particles['h']

        self.particles['rho'] = np.nan

        if is_full_dump:
            self.particles['P'] = particles['pressure']
        else:
            self.particles['P'] = np.nan

        if is_full_dump:
            self.particles['vx'] = particles['vxyz'][:, 0]
            self.particles['vy'] = particles['vxyz'][:, 1]
            self.particles['vz'] = particles['vxyz'][:, 2]

        else:
            self.particles['vx'] = np.nan
            self.particles['vy'] = np.nan
            self.particles['vz'] = np.nan

        if n_dust_types > 0:
            for n in range(n_dust_types):
                tag1 = 'dustfrac' + str(n+1)
                # tag2 = 'tstop' + str(n+1)
                self.particles[tag1] = particles['dustfrac'][:, n]
                # if is_full_dump:
                #     self.particles[tag2] = particles['tstop'][:, n]
                # else:
                #     self.particles[tag2] = np.nan

        if is_full_dump:
            self.particles['divv']  = particles['divv']
            self.particles['dt']    = particles['dt']

        else:
            self.particles['divv']  = np.nan
            self.particles['dt']    = np.nan

        self.particles['itype'] = particles['itype']

        self.particles.loc[self.particles['itype']==I_GAS, 'm'] = \
            self.parameters['massoftype'][I_GAS-1]

        self.particles['rho'] = density_from_smoothing_length(
            self.particles['h'],
            self.particles['m'],
            hfact=self.parameters['hfact'])

        if n_dust_large  > 0:
            for n in range(n_dust_large):
                self.particles.loc[self.particles['itype']==I_DUST+n, 'm'] = \
                    self.parameters['massoftype'][I_DUST+n-1]

        #--- Sinks

        if contains_sinks:

            sinks = f['sinks']

            self.sinks = \
                    pd.DataFrame( sinks['xyz'][:n_sinks, 0], columns=['x'] )
            self.sinks['y'] = sinks['xyz'][:n_sinks, 1]
            self.sinks['z'] = sinks['xyz'][:n_sinks, 2]
            self.sinks['m'] = sinks['m'][:n_sinks]
            self.sinks['h'] = sinks['h'][:n_sinks]
            self.sinks['hsoft'] = sinks['hsoft'][:n_sinks]
            self.sinks['macc']  = sinks['maccreted'][:n_sinks]
            self.sinks['spinx'] = sinks['spinxyz'][:n_sinks, 0]
            self.sinks['spiny'] = sinks['spinxyz'][:n_sinks, 1]
            self.sinks['spinz'] = sinks['spinxyz'][:n_sinks, 2]
            self.sinks['tlast'] = sinks['tlast'][:n_sinks]

        if is_full_dump:
            self.sinks['vx'] = sinks['vxyz'][:n_sinks, 0]
            self.sinks['vy'] = sinks['vxyz'][:n_sinks, 1]
            self.sinks['vz'] = sinks['vxyz'][:n_sinks, 2]

        else:
            self.sinks['vx'] = np.nan
            self.sinks['vy'] = np.nan
            self.sinks['vz'] = np.nan

        f.close()

    #--- Read ascii dump

    def _read_ascii(self, file_prefix):

        header_file_name = file_prefix + '.header'
        dump_file_name = file_prefix + '.ascii'

        exists = os.path.isfile(header_file_name)
        if not exists:
            raise FileNotFoundError('Cannot find header file: ' +
                                    header_file_name)

        parameters, is_full_dump = _read_header_from_showheader(header_file_name)

        self.parameters = parameters

        self.units = Units(parameters['udist'],
                           parameters['umass'],
                           parameters['utime']).units

        exists = os.path.isfile(dump_file_name)
        if not exists:
            raise FileNotFoundError('Cannot find dump file: ' + dump_file_name)

        n_dust_types = self.parameters['ndustsmall'] + self.parameters['ndustlarge']
        n_dust_large = self.parameters['ndustlarge']

        n_sinks = self.parameters['nptmass']
        contains_sinks = bool(n_sinks > 0)

        names = ['x', 'y', 'z', 'm', 'h', 'rho']
        sink_drop = list()
        if is_full_dump:
            names += ['vx', 'vy', 'vz']
        if n_dust_types > 0:
            for n in range(n_dust_types):
                names += ['dustfrac' + str(n+1)]
                sink_drop += ['dustfrac' + str(n+1)]
        if is_full_dump:
            names += ['divv', 'dt', 'itype']
            sink_drop += ['divv', 'dt', 'itype']
        else:
            names += ['itype']
            sink_drop += ['itype']

        print_warning('Assuming ascii file columns are as follows:\n' + \
                      ', '.join(names))

        data = pd.read_csv(dump_file_name, comment='#', names=names,
                           delim_whitespace=True)

        particles = data[data['itype'] != I_SINK].reset_index(drop=True)

        particles.loc[
            (particles['itype'] >= I_DUST_SPLASH) &
            (particles['itype'] <= I_DUST_SPLASH + n_dust_large),
            'itype'] -= I_DUST_SPLASH - I_DUST

        self.particles = particles

        if contains_sinks:

            sinks = data[data['itype'] == I_SINK].reset_index(drop=True)
            sinks = sinks.drop(sink_drop, axis=1)

            self.sinks = sinks

#--- Functions

def _read_header_from_showheader(header_file_name):

    with open(header_file_name, 'r') as file:

        keys = list()
        values = list()

        first_line = True

        for line in file:

            if first_line:

                if line[0] == 'F':
                    is_full_dump = True
                elif line[0] == 'S':
                    is_full_dump = False
                else:
                    raise ValueError('Cannot determine dump type')

                first_line = False
                continue

            keys.append(line.rstrip('\n').split()[0])
            value = line.rstrip('\n').split()[1]

            if '.' in value:
                value = float(value)
            else:
                value = int(value)

            values.append(value)

    new_keys = list()
    new_values = list()

    prev_key = None
    first_multiple_key = True
    idxj = 0

    for idxi, key in enumerate(keys):

        if key == prev_key:

            if first_multiple_key:
                new_values[idxj-1] = list()
                new_values[idxj-1].append(values[idxi-1])
                new_values[idxj-1].append(values[idxi])
                first_multiple_key = False
            else:
                new_values[idxj-1].append(values[idxi])

        else:

            first_multiple_key = True
            idxj += 1
            new_keys.append(key)
            new_values.append(values[idxi])

        prev_key = key

    header = dict(zip(new_keys, new_values))

    for key in header:
        if isinstance(header[key], list):
            header[key] = np.array(header[key])

    return header, is_full_dump
