"""
Simulation class for smoothed particle hydrodynamics simulations.

This class contains all information associated with a simulation. It
contains the dumps as a list of Dump objects, and the global and sink
evolution files an Evolution objects.
"""

from pathlib import Path

from ..dump.dump import Dump
from .evolution import Evolution


class Simulation:
    """
    Smoothed particle hydrodynamics simulation object.

    This class aggregates dump files and time evolution data. Dump files
    contain a snapshot of the simulation at a particular time. Time
    evolution files contain time series of global quantities such as
    energy and momentum.

    Parameters
    ----------
    prefix : str
        Simulation prefix, e.g. 'disc', if files are named like
        disc_00000.h5, disc01.ev, discSink0001N01.ev, etc.

    directory : str, optional
        Directory containing simulation dump files and auxiliary files.

    **kwargs
        Arguments to pass to Dump.

    Examples
    --------
    Reading simulation data into a Simulation object.

    >>> simulation = plonk.load_simulation('prefix', path_to_directory)

    Accessing the particles object, and particles arrays.

    >>> simulation.particles
    >>> simulation.particles.arrays

    Accessing the global time evolution and sink time evolution objects.

    >>> simulation.evolution
    >>> simulation.sink_evolution
    """

    def __init__(self, prefix, directory=None, **kwargs):

        self._init_kwargs = kwargs

        if not isinstance(prefix, str):
            raise TypeError('prefix must be str')

        if directory is None:
            directory = '.'
        else:
            if not isinstance(directory, (str, Path)):
                raise TypeError('directory must be str or pathlib.Path')

        self.prefix = prefix
        self.path = Path(directory).expanduser().resolve()

        if not list(self.path.glob(self.prefix + '*')):
            raise FileNotFoundError(f'No files with prefix: {prefix}')

        self._dump_file_type, self._dump_file_extension = self._get_dump_file_type()

        self._dump_files = self._get_dump_files()
        self._global_ev_files = self._get_global_ev_files()
        self._sink_ev_files = self._get_sink_ev_files()

        self._dumps = None
        self._evolution = None
        self._sinks_evolution = None

    @property
    def dumps(self):
        """List of Dump objects associated with the simulation."""

        if self._dumps is None:
            self._generate_dump_objects()

        return self._dumps

    @property
    def evolution(self):
        """Global time evolution data."""

        if self._evolution is None:
            self._generate_evolution_object()

        return self._evolution

    @property
    def sink_evolution(self):
        """Sink time evolution data."""

        if self._sinks_evolution is None:
            self._generate_sink_evolution_objects()

        return self._sinks_evolution

    def _generate_dump_objects(self):
        """Generate dump objects."""

        dumps = list()
        for dump in self._dump_files:
            dumps.append(Dump(dump, **self._init_kwargs))
        self._dumps = dumps

    def _generate_evolution_object(self):
        """Generate global evolution object."""
        self._evolution = None
        if self._global_ev_files:
            self._evolution = Evolution(self._global_ev_files)

    def _generate_sink_evolution_objects(self):
        """Generate sink evolution objects."""

        self._sinks_evolution = list()
        [
            self._sinks_evolution.append(Evolution(files))
            for files in self._sink_ev_files
        ]

    def _get_global_ev_files(self, glob=None):
        """Get global time evolution files."""

        if glob is None:
            # Phantom ev file name format
            glob = self.prefix + '[0-9][0-9].ev'

        return sorted(list(self.path.glob(glob)))

    def _get_sink_ev_files(self, glob=None):
        """Get time evolution files for sinks.

        Returns a list of list of pathlib.Path objects. The inner lists
        contain each sink evolution file for each sink.
        """

        if glob is None:
            # Phantom ev file name format
            glob = self.prefix + 'Sink[0-9][0-9][0-9][0-9]N[0-9][0-9].ev'

        n = len(self.prefix) + len('Sink')
        n_sinks = len(set([p.name[n : n + 4] for p in list(self.path.glob(glob))]))

        sinks = list()
        for idx in range(1, n_sinks + 1):
            sinks.append(
                sorted(
                    list(self.path.glob(self.prefix + f'Sink{idx:04}N[0-9][0-9].ev'))
                )
            )

        return sinks

    def _get_dump_files(self, glob=None):
        """Get dump files."""

        if glob is None:
            # Phantom dump file name format
            glob = (
                self.prefix + '_[0-9][0-9][0-9][0-9][0-9].' + self._dump_file_extension
            )

        return sorted(list(self.path.glob(glob)))

    def _get_dump_file_type(self, glob=None):
        """
        Determine dump file type from extension assuming file names
        follow a glob pattern.
        """

        if glob is None:
            # Phantom HDF5 dump file name format
            glob = self.prefix + '_[0-9][0-9][0-9][0-9][0-9].h5'

        file_types = set([f.suffix for f in self.path.glob(glob)])

        if len(file_types) > 1:
            raise ValueError(
                'Cannot determine simulation dump file type: '
                f'is it one of {file_types}?'
            )
        elif len(file_types) == 0:
            raise ValueError(
                'Cannot determine dump file type: '
                'no files named like prefix_xxxxx.ext'
            )

        file_ext = file_types.pop()[1:]

        if file_ext == 'h5':
            file_type = 'HDF5'

        return file_type, file_ext

    def __repr__(self):
        return self.__str__()

    def __str__(self):
        return f'<plonk.Simulation: "{self.prefix}", ' f'directory="{self.path.name}">'


def load_sim():
    raise NotImplementedError
