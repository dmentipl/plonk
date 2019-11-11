"""Simulation class for smoothed particle hydrodynamics simulations.

This class contains all information associated with a simulation. It
contains the snapshot as a list of Snap objects, and the global quantity
and sink quantity time series files an Evolution objects.
"""

from __future__ import annotations

from pathlib import Path
from typing import List, Optional, Tuple, Union

from ..snap import Snap, load_snap
from .evolution import Evolution, load_ev


class Simulation:
    """Smoothed particle hydrodynamics simulation object.

    This class aggregates snapshot files, global quantity and sink time
    series data. Snapshot files contain a complete snapshot of the
    simulation at a particular time. Other files contain time series of
    global quantities on particles such as energy and momentum, and
    time series data for sink particles.

    Examples
    --------
    Reading simulation data into a Simulation object.

    >>> sim = plonk.load_sim('prefix', path_to_directory)

    Accessing the snapshots.

    >>> sim.snaps

    Accessing the global quantity and sink time series data.

    >>> sim.global_quantities
    >>> sim.sink_quantities
    """

    def __init__(self):

        self.prefix: str
        self.path: Path

        self._snaps: List[Snap] = None
        self._global_quantities: Evolution = None
        self._sink_quantities: List[Evolution] = None

    def load_sim(
        self, prefix: str, directory: Optional[Union[str, Path]] = None
    ) -> Simulation:
        """Load Simulation.

        Parameters
        ----------
        prefix
            Simulation prefix, e.g. 'disc', if files are named like
            disc_00000.h5, disc01.ev, discSink0001N01.ev, etc.
        directory
            Directory containing simulation snapshot files and auxiliary
            files.
        """
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

        self._snap_file_type, self._snap_file_extension = self._get_snap_file_type()

        self._snap_files = self._get_snap_files()
        self._global_ev_files = self._get_global_ev_files()
        self._sink_ev_files = self._get_sink_ev_files()

        return self

    @property
    def snaps(self) -> List[Snap]:
        """List of Snap objects associated with the simulation."""
        if self._snaps is None:
            self._generate_snap_objects()

        return self._snaps

    @property
    def global_quantities(self) -> Evolution:
        """Global quantity time series data."""
        if self._global_quantities is None:
            self._generate_global_quantities()

        return self._global_quantities

    @property
    def sink_quantities(self) -> List[Evolution]:
        """Sink time series data."""
        if self._sink_quantities is None:
            self._generate_sink_quantities()

        return self._sink_quantities

    def _generate_snap_objects(self):
        """Generate Snap objects."""
        snaps = list()
        for snap in self._snap_files:
            snaps.append(load_snap(snap))
        self._snaps = snaps

    def _generate_global_quantities(self):
        """Generate global quantity time series objects."""
        self._global_quantities = None
        if self._global_ev_files:
            self._global_quantities = load_ev(self._global_ev_files)

    def _generate_sink_quantities(self):
        """Generate sink quantity time series objects."""
        self._sink_quantities = list()
        [self._sink_quantities.append(load_ev(files)) for files in self._sink_ev_files]

    def _get_global_ev_files(self, glob: str = None) -> List[Path]:
        """Get global ev files."""
        if glob is None:
            # Phantom ev file name format
            glob = self.prefix + '[0-9][0-9].ev'

        return sorted(list(self.path.glob(glob)))

    def _get_sink_ev_files(self, glob: str = None) -> List[List[Path]]:
        """Get sink ev files."""
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

    def _get_snap_files(self, glob: str = None) -> List[Path]:
        """Get snapshot files."""
        if glob is None:
            # Phantom snapshot file name format
            glob = (
                self.prefix + '_[0-9][0-9][0-9][0-9][0-9].' + self._snap_file_extension
            )

        return sorted(list(self.path.glob(glob)))

    def _get_snap_file_type(self, glob: str = None) -> Tuple[str, str]:
        """Snapshot file type.

        Determine snap file type from extension assuming file names
        follow a glob pattern.
        """
        if glob is None:
            # Phantom HDF5 snapshot file name format
            glob = self.prefix + '_[0-9][0-9][0-9][0-9][0-9].h5'

        file_types = set([f.suffix for f in self.path.glob(glob)])

        if len(file_types) > 1:
            raise ValueError(
                'Cannot determine simulation snapshot file type: '
                f'is it one of {file_types}?'
            )
        elif len(file_types) == 0:
            raise ValueError(
                'Cannot determine snapshot file type: '
                'no files named like prefix_xxxxx.ext'
            )

        file_ext = file_types.pop()[1:]

        if file_ext == 'h5':
            file_type = 'HDF5'

        return file_type, file_ext

    def __repr__(self):
        """Dunder repr method."""
        return self.__str__()

    def __str__(self):
        """Dunder str method."""
        return f'<plonk.Simulation: "{self.prefix}", ' f'directory="{self.path.name}">'


def load_sim(prefix: str, directory: Optional[Union[str, Path]] = None) -> Simulation:
    """Load Simulation.

    Parameters
    ----------
    prefix
        Simulation prefix, e.g. 'disc', if files are named like
        disc_00000.h5, disc01.ev, discSink0001N01.ev, etc.
    directory
        Directory containing simulation snapshot files and auxiliary
        files.
    """
    return Simulation().load_sim(prefix=prefix, directory=directory)
