'''
dump.py

Daniel Mentiplay, 2019.
'''

import os

import h5py
import numpy as np
import pandas as pd

from .particles import density_from_smoothing_length, I_DUST
from .units import Units
from .utils import print_warning

# ---------------------------------------------------------------------------- #

#--- Reading Phantom-Splash ACSII files

I_SINK         = 3
I_DUST_SPLASH  = 8

# ---------------------------------------------------------------------------- #

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
            dump_file_format = 'HDF5'
            print_warning('HDF5 dump reader not fully implemented')

        elif file_extension == 'ascii':
            dump_file_format = 'ASCII'

        else:
            raise ValueError('Cannot determine dump file format')

        is_full_dump = self._read_header(file_prefix, dump_file_format)
        self._read_arrays(file_prefix, dump_file_format, is_full_dump)

    #--- Read header

    def _read_header(self, file_prefix, dump_file_format):

        if dump_file_format == 'HDF5':
            file_extension = 'h5'

        elif dump_file_format == 'ASCII':
            file_extension = 'header'

        header_file_name = file_prefix + '.' + file_extension

        exists = os.path.isfile(header_file_name)
        if not exists:
            raise FileNotFoundError('Cannot find header file: ' +
                                    header_file_name)

        if dump_file_format == 'HDF5':
            header, is_full_dump = _read_header_from_hdf5(header_file_name)

        elif dump_file_format == 'ASCII':
            header, is_full_dump = _read_header_from_showheader(header_file_name)

        self.parameters = header

        self.units = Units(header['udist'], header['umass'],
                           header['utime']).units

        return is_full_dump

    #--- Read arrays

    def _read_arrays(self, file_prefix, dump_file_format, is_full_dump):

        if dump_file_format == 'HDF5':
            file_extension = 'h5'

        elif dump_file_format == 'ASCII':
            file_extension = 'ascii'

        dump_file_name = file_prefix + '.' + file_extension

        exists = os.path.isfile(dump_file_name)
        if not exists:
            raise FileNotFoundError('Cannot find dump file: ' + dump_file_name)

        if dump_file_format == 'HDF5':

            self._read_arrays_from_hdf5(dump_file_name, is_full_dump)

        elif dump_file_format == 'ASCII':

            self._read_arrays_from_ascii(dump_file_name, is_full_dump)

    def _read_arrays_from_ascii(self, dump_file_name, is_full_dump):

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

    def _read_arrays_from_hdf5(self, dump_file_name, is_full_dump):

        n_dust_types = self.parameters['ndustsmall'] + self.parameters['ndustlarge']

        n_sinks = self.parameters['nptmass']
        contains_sinks = bool(n_sinks > 0)

        f = h5py.File(dump_file_name, 'r')

        #--- Particles

        particles = f['particles']

        self.particles = \
            pd.DataFrame( particles['xyz'][:, 0], columns=['x'] )
        self.particles['y'] = particles['xyz'][:, 1]
        self.particles['z'] = particles['xyz'][:, 2]

        # TODO: set particle masses
        self.particles['m'] = 1.
        print_warning('Mass not set properly')

        self.particles['h'] = particles['h']

        self.particles['rho'] = density_from_smoothing_length(
            self.particles['h'],
            self.particles['m'],
            hfact=self.parameters['hfact'])

        self.particles['P'] = particles['pressure']

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
                tag2 = 'tstop' + str(n+1)
                self.particles[tag1] = particles['dustfrac'][:, n]
                self.particles[tag2] = particles['tstop'][:, n]

        if is_full_dump:
            self.particles['divv']  = particles['divv']
            self.particles['dt']    = particles['dt']
            self.particles['itype'] = particles['itype']

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

#--- Functions

def _read_header_from_hdf5(dump_file_name):

    f = h5py.File(dump_file_name, 'r')

    header_group = f['header']

    header = dict()
    for key in header_group.keys():
        header[key] = header_group[key].value

    f.close()

    is_full_dump = True
    print_warning('HDF5 dump reader work in progress: assuming full dump')

    return header, is_full_dump

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
