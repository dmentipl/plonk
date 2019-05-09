"""
dump.py

Daniel Mentiplay, 2019.
"""

import os

import h5py
import numpy as np
import pandas as pd

from .particles import I_DUST, I_GAS, density_from_smoothing_length
from .units import Units
from .utils import print_warning

# --- Possibly available arrays in Phantom dump

POSSIBLE_ARRAYS = [
    'itype',
    'xyz',
    'h',
    'pressure',
    'vxyz',
    'u',
    'dustfrac',
    'tstop',
    'deltavxyz',
    'divv',
    'curlvxyz',
    'dt',
    'alpha',
    'poten',
    'grainsize',
    'graindens',
    'vrel/vfrag',
    'St',
    'abundance',
    'T',
    'luminosity',
    'beta_pr',
    'Bxyz',
    'psi',
    'divB',
    'curlBxyz',
    'divBsymm',
    'eta_{OR}',
    'eta_{HE}',
    'eta_{AD}',
    'ne/n',
]

NABUNDANCES = 5

# --- Reading Phantom-Splash ACSII files

I_SINK = 3
I_DUST_SPLASH = 8


class Dump:
    """
    Smoothed particle hydrodynamics dump file object. The main data stored
    in this object is the particle data, e.g. positions, velocities,
    density, smoothing length, and so on, stored as a Pandas DataFrame.

    Parameters
    ----------
    name : str
        Path to dump file.

    Examples
    --------
    Reading a dump file into a Dump object.

    >>> file_name = 'dumpfile.ext'
    >>> dump = plonk.Dump(file_name)
    """

    def __init__(self, filename):

        self.particles = None
        self.sinks = None
        self.parameters = None
        self.units = None

        self._filename = filename.split('/')[-1]
        self._filepath = os.path.abspath(filename)

        # TODO: only Phantom works currently
        self._simulation_code = 'Phantom'

        self._read_dump(filename)

    def __repr__(self):
        return self.__str__()

    def __str__(self):

        prefix = os.path.splitext(self._filename)[0]
        filename = self._filename
        directory = os.path.dirname(self._filepath)
        size = _file_size(self._filepath)
        time = self.parameters['time']
        ntotal = self.particles.shape[0] + self.sinks.shape[0]
        ngas = self.particles[self.particles['itype'] == 1].shape[0]
        ndust = self.particles[self.particles['itype'] > 1].shape[0]
        nsinks = self.sinks.shape[0]

        return (
            f'--- {prefix} ---\n\n'
            + f'filename:  {filename}'.join(['  ', '\n'])
            + f'directory: {directory}'.join(['  ', '\n'])
            + f'size:      {size}'.join(['  ', '\n'])
            + f'time:      {time}'.join(['  ', '\n'])
            + '\n'
            + 'particles'.join(['  ', '\n'])
            + '---------'.join(['  ', '\n'])
            + f'total: {ntotal}'.join(['  ', '\n'])
            + f'gas:   {ngas}'.join(['  ', '\n'])
            + f'dust:  {ndust}'.join(['  ', '\n'])
            + f'sinks: {nsinks}'.join(['  ', '\n'])
        )

    # --- Read dump

    def _read_dump(self, filename):
        """
        Read a Phantom dump file.

        Two file types can be read: ASCII or HDF.

        1. For HDF only one file is required and it must have the extension h5,
        e.g. 'disc_00000.h5'.

        2. For ASCII dumps two files are required:
            (i)  the dump file in ASCII format, e.g. 'disc_00000.ascii'
            (ii) the header file in "showheader" format, e.g.
            'disc_00000.header'

        Arguments:
            filename : e.g. 'disc_00000.h5' for HDF or 'disc_00000.ascii' for
            ASCII
        """

        exists = os.path.isfile(filename)

        if not exists:
            raise FileNotFoundError('Cannot find dump file')

        file_prefix = filename.split('.')[0]
        file_extension = filename.split('.')[-1]

        if file_extension == 'h5':
            self._read_hdf(file_prefix)

        elif file_extension == 'ascii':
            print_warning('ASCII files are deprecated')
            self._read_ascii(file_prefix)

        else:
            raise ValueError('Cannot determine dump file format')

    # --- Read HDF dump

    def _read_hdf(self, file_prefix):

        # --- Open file

        dump_file_name = file_prefix + '.h5'

        exists = os.path.isfile(dump_file_name)

        if not exists:
            raise FileNotFoundError('Cannot find dump file: ' + dump_file_name)

        f = h5py.File(dump_file_name, 'r')

        # --- Header

        header = f['header']

        parameters = dict()
        for key in header.keys():
            parameters[key] = header[key][()]

        self.parameters = parameters

        self.units = Units(
            parameters['udist'], parameters['umass'], parameters['utime']
        ).units

        n_dust_large = self.parameters['ndustlarge']
        n_dust_types = (
            self.parameters['ndustsmall'] + self.parameters['ndustlarge']
        )

        n_sinks = self.parameters['nptmass']
        contains_sinks = bool(n_sinks > 0)

        # --- Particles

        particles = f['particles']

        is_full_dump = bool('vxyz' in particles)

        self.particles = pd.DataFrame(
            particles['itype'][()], columns=['itype']
        )

        self.particles['x'] = particles['xyz'][:, 0]
        self.particles['y'] = particles['xyz'][:, 1]
        self.particles['z'] = particles['xyz'][:, 2]
        self.particles['h'] = particles['h']

        self.particles.loc[
            self.particles['itype'] == I_GAS, 'm'
        ] = self.parameters['massoftype'][I_GAS - 1]

        if n_dust_large > 0:
            for n in range(n_dust_large):
                self.particles.loc[
                    self.particles['itype'] == I_DUST + n, 'm'
                ] = self.parameters['massoftype'][I_DUST + n - 1]

        self.particles['rho'] = density_from_smoothing_length(
            self.particles['h'],
            self.particles['m'],
            hfact=self.parameters['hfact'],
        )

        generator = (arr for arr in POSSIBLE_ARRAYS[3:] if arr in particles)

        for array in generator:

            if particles[array].size == 0:
                continue

            elif particles[array].ndim == 1:
                self.particles[array] = particles[array]

            elif particles[array].ndim == 2:

                if array[-3:] == 'xyz':
                    columns = [array[:-3] + p for p in ['x', 'y', 'z']]

                else:
                    if n_dust_types > 1:
                        columns = [
                            array + str(i + 1) for i in range(n_dust_types)
                        ]
                    else:
                        columns = [array]

                for idx, column in enumerate(columns):
                    self.particles[column] = particles[array][:, idx]

            elif particles[array].ndim == 3:

                if array[-3:] == 'xyz':
                    columns_ = [array[:-3] + p for p in ['x', 'y', 'z']]

                    if n_dust_types > 1:
                        columns = list()
                        for idx, col_ in enumerate(columns_):
                            columns.append(
                                [col_ + str(i + 1) for i in range(n_dust_types)]
                            )

                    else:
                        columns = list()
                        for idx, column_ in enumerate(columns_):
                            columns.append([column_])

                    for ind_pos, column in enumerate(columns):
                        for ind_grain, column_ in enumerate(column):
                            self.particles[column_] = particles[array][
                                :, ind_grain, ind_pos
                            ]

                else:
                    raise ValueError(f'Cannot read array: {array}')

            else:
                raise ValueError(f'Cannot read array: {array}')

        # --- Sinks

        if contains_sinks:

            sinks = f['sinks']

            columns = [
                'x',
                'y',
                'z',
                'm',
                'h',
                'hsoft',
                'macc',
                'spinx',
                'spiny',
                'spinz',
                'tlast',
            ]

            if is_full_dump:
                columns += ['vx', 'vy', 'vz']

            self.sinks = pd.DataFrame(columns=columns)

            self.sinks['x'] = sinks['xyz'][:n_sinks, 0]
            self.sinks['y'] = sinks['xyz'][:n_sinks, 1]
            self.sinks['z'] = sinks['xyz'][:n_sinks, 2]
            self.sinks['m'] = sinks['m'][:n_sinks]
            self.sinks['h'] = sinks['h'][:n_sinks]
            self.sinks['hsoft'] = sinks['hsoft'][:n_sinks]
            self.sinks['macc'] = sinks['maccreted'][:n_sinks]
            self.sinks['spinx'] = sinks['spinxyz'][:n_sinks, 0]
            self.sinks['spiny'] = sinks['spinxyz'][:n_sinks, 1]
            self.sinks['spinz'] = sinks['spinxyz'][:n_sinks, 2]
            self.sinks['tlast'] = sinks['tlast'][:n_sinks]

        if is_full_dump:
            self.sinks['vx'] = sinks['vxyz'][:n_sinks, 0]
            self.sinks['vy'] = sinks['vxyz'][:n_sinks, 1]
            self.sinks['vz'] = sinks['vxyz'][:n_sinks, 2]

        f.close()

    # --- Read ascii dump

    def _read_ascii(self, file_prefix):

        header_file_name = file_prefix + '.header'
        dump_file_name = file_prefix + '.ascii'

        exists = os.path.isfile(header_file_name)
        if not exists:
            raise FileNotFoundError(
                'Cannot find header file: ' + header_file_name
            )

        parameters, is_full_dump = _read_header_from_showheader(
            header_file_name
        )

        self.parameters = parameters

        self.units = Units(
            parameters['udist'], parameters['umass'], parameters['utime']
        ).units

        exists = os.path.isfile(dump_file_name)
        if not exists:
            raise FileNotFoundError('Cannot find dump file: ' + dump_file_name)

        n_dust_types = (
            self.parameters['ndustsmall'] + self.parameters['ndustlarge']
        )
        n_dust_large = self.parameters['ndustlarge']

        n_sinks = self.parameters['nptmass']
        contains_sinks = bool(n_sinks > 0)

        names = ['x', 'y', 'z', 'm', 'h', 'rho']
        sink_drop = list()
        if is_full_dump:
            names += ['vx', 'vy', 'vz']
        if n_dust_types > 0:
            for n in range(n_dust_types):
                names += ['dustfrac' + str(n + 1)]
                sink_drop += ['dustfrac' + str(n + 1)]
        if is_full_dump:
            names += ['divv', 'dt', 'itype']
            sink_drop += ['divv', 'dt', 'itype']
        else:
            names += ['itype']
            sink_drop += ['itype']

        print_warning(
            'Assuming ascii file columns are as follows:\n' + ', '.join(names)
        )

        data = pd.read_csv(
            dump_file_name, comment='#', names=names, delim_whitespace=True
        )

        particles = data[data['itype'] != I_SINK].reset_index(drop=True)

        particles.loc[
            (particles['itype'] >= I_DUST_SPLASH)
            & (particles['itype'] <= I_DUST_SPLASH + n_dust_large),
            'itype',
        ] -= (I_DUST_SPLASH - I_DUST)

        self.particles = particles

        if contains_sinks:

            sinks = data[data['itype'] == I_SINK].reset_index(drop=True)
            sinks = sinks.drop(sink_drop, axis=1)

            self.sinks = sinks


# --- Functions


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
                new_values[idxj - 1] = list()
                new_values[idxj - 1].append(values[idxi - 1])
                new_values[idxj - 1].append(values[idxi])
                first_multiple_key = False
            else:
                new_values[idxj - 1].append(values[idxi])

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


def _file_size(filepath):
    """Get file size. Returns tuple (size (int), unit (str))."""

    size_names = ('B', 'KB', 'MB', 'GB', 'TB', 'PB')
    size_maxes = (1e3, 1e6, 1e9, 1e12, 1e15)

    size_in_bytes = os.stat(filepath).st_size

    for size_name, size_max in zip(size_names, size_maxes):
        if size_in_bytes < size_max:
            return (int(size_in_bytes * 1e3 / size_max), size_name)

    return (int(size_in_bytes / size_maxes[-1]), size_names[-1])
