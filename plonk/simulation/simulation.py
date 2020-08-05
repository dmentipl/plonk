"""Simulation class for smoothed particle hydrodynamics simulations.

This class contains all information associated with a simulation. It
contains the snapshot as a list of Snap objects, and the global quantity
and sink quantity time series files as pandas DataFrames.
"""

from __future__ import annotations

from copy import copy
from pathlib import Path
from typing import TYPE_CHECKING, Any, Dict, List, Union

import numpy as np
from pandas import DataFrame

from .._logging import logger
from .._units import Quantity
from ..snap.readers import load_snap
from .evolution import load_ev

if TYPE_CHECKING:
    from ..snap.snap import Snap

_properties_vary_per_snap = ('time',)
_data_sources = ('Phantom',)


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

    Accessing the properties.

    >>> sim.properties

    Accessing the global quantity and sink time series data.

    >>> sim.global_quantities
    >>> sim.sink_quantities
    """

    def __init__(self):

        self.data_source: str
        self.prefix: str
        self.paths: Dict[str, Any]

        self._snaps: List[Snap] = None
        self._properties: Dict[str, Any] = None
        self._units: Dict[str, Any] = None
        self._global_quantities: DataFrame = None
        self._sink_quantities: List[DataFrame] = None

        self._snap_file_extension = ''
        self._len = -1

    def load_sim(
        self,
        prefix: str,
        directory: Union[str, Path] = None,
        data_source: str = 'Phantom',
    ) -> Simulation:
        """Load Simulation.

        Parameters
        ----------
        prefix
            Simulation prefix, e.g. 'disc', if files are named like
            disc_00000.h5, disc01.ev, discSink0001N01.ev, etc.
        directory : optional
            Directory containing simulation snapshot files and auxiliary
            files. Default is None.
        data_source : optional
            The SPH code used to produce the simulation data. Default
            is 'Phantom'.
        """
        if not isinstance(prefix, str):
            raise TypeError('prefix must be str')

        if data_source not in _data_sources:
            raise ValueError(f'Data source not available: try {_data_sources}')
        self.data_source = data_source

        if directory is None:
            directory = '.'
        else:
            if not isinstance(directory, (str, Path)):
                raise TypeError('directory must be str or pathlib.Path')

        logger.debug(f'Loading {data_source} simulation: {prefix} at {directory}')
        self.prefix = prefix
        self.paths = {
            'directory': Path(directory).expanduser().resolve(),
        }

        if not list(self.paths['directory'].glob(self.prefix + '*')):
            raise FileNotFoundError(f'No files with prefix: {prefix}')

        self._snap_file_extension = self._get_snap_file_extension()

        self.paths['snaps'] = self._get_snap_files()
        self.paths['global_quantities'] = self._get_global_ev_files()
        self.paths['sink_quantities'] = self._get_sink_ev_files()

        return self

    @property
    def snaps(self) -> List[Snap]:
        """List of Snap objects associated with the simulation."""
        if self._snaps is None:
            self._generate_snap_objects()

        return self._snaps

    @property
    def properties(self) -> Dict[str, Any]:
        """Properties associated with the simulation."""
        if self._properties is None:
            self._generate_properties()

        return self._properties

    @property
    def units(self) -> Dict[str, Any]:
        """Units associated with the simulation."""
        if self._units is None:
            self._generate_units()

        return self._units

    @property
    def global_quantities(self) -> DataFrame:
        """Global quantity time series data."""
        if self._global_quantities is None:
            self._generate_global_quantities()

        return self._global_quantities

    @property
    def sink_quantities(self) -> List[DataFrame]:
        """Sink time series data."""
        if self._sink_quantities is None:
            self._generate_sink_quantities()

        return self._sink_quantities

    def _generate_snap_objects(self):
        """Generate Snap objects."""
        snaps = list()
        fail = 0
        for snap in self.paths['snaps']:
            try:
                snaps.append(load_snap(snap))
            except (OSError, RuntimeError):
                fail += 1
        if fail > 0:
            logger.warning(f'Cannot read {fail} snap(s)')
        self._snaps = snaps

    def _generate_properties(self):
        """Generate sim.properties from snap.properties."""
        prop = copy(self.snaps[0].properties)
        for key in _properties_vary_per_snap:
            prop[key] = list()
        for snap in self.snaps:
            for key, val in snap.properties.items():
                if isinstance(prop[key], list):
                    prop[key].append(val)
                else:
                    if np.any(prop[key] != val):
                        prop[key] = '__inconsistent__'
        for key, val in prop.items():
            if isinstance(val, list):
                if isinstance(val[0], Quantity):
                    prop[key] = np.array([v.m for v in val]) * val[0].u
                else:
                    prop[key] = np.array(val)
        self._properties = prop

    def _generate_units(self):
        """Generate sim.units from snap.units."""
        u = copy(self.snaps[0].units)
        for snap in self.snaps:
            for key, val in snap.units.items():
                if u[key] != val:
                    u[key] = '__inconsistent__'
        self._units = u

    def _generate_global_quantities(self):
        """Generate global quantity time series objects."""
        self._global_quantities = None
        if self.paths['global_quantities']:
            self._global_quantities = load_ev(self.paths['global_quantities'])

    def _generate_sink_quantities(self):
        """Generate sink quantity time series objects."""
        self._sink_quantities = [
            load_ev(files) for files in self.paths['sink_quantities']
        ]

    def _get_global_ev_files(self, glob: str = None) -> List[Path]:
        """Get global ev files."""
        if glob is None:
            # Phantom ev file name format
            glob = self.prefix + '[0-9][0-9].ev'

        return sorted(list(self.paths['directory'].glob(glob)))

    def _get_sink_ev_files(self, glob: str = None) -> List[List[Path]]:
        """Get sink ev files."""
        if glob is None:
            # Phantom ev file name format
            glob = self.prefix + 'Sink[0-9][0-9][0-9][0-9]N[0-9][0-9].ev'

        n = len(self.prefix) + len('Sink')
        n_sinks = len({p.name[n : n + 4] for p in self.paths['directory'].glob(glob)})

        sinks = list()
        for idx in range(1, n_sinks + 1):
            sinks.append(
                sorted(
                    self.paths['directory'].glob(
                        self.prefix + f'Sink{idx:04}N[0-9][0-9].ev'
                    )
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

        return sorted(list(self.paths['directory'].glob(glob)))

    def _get_snap_file_extension(self, glob: str = None):
        """Snapshot file extension.

        Determine snap file type from extension assuming file names
        follow a glob pattern.
        """
        if glob is None:
            # Phantom HDF5 snapshot file name format
            glob = self.prefix + '_[0-9][0-9][0-9][0-9][0-9].h5'

        file_types = {f.suffix for f in self.paths['directory'].glob(glob)}

        if len(file_types) > 1:
            raise ValueError(
                'Cannot determine simulation snapshot file type: '
                f'is it one of {file_types}?'
            )
        if len(file_types) == 0:
            raise ValueError(
                'Cannot determine snapshot file type: '
                'no files named like prefix_xxxxx.ext'
            )

        file_ext = file_types.pop()[1:]

        if file_ext not in ('h5',):
            raise ValueError('File extension not available; must be ".h5"')

        return file_ext

    def __len__(self):
        """Length as number of snaps."""
        if self._len == -1:
            self._len = len(self.snaps)
        return self._len

    def __repr__(self):
        """Dunder repr method."""
        return self.__str__()

    def __str__(self):
        """Dunder str method."""
        return (
            f'<plonk.Simulation: "{self.prefix}", '
            f'directory="{self.paths["directory"].name}">'
        )


def load_sim(
    prefix: str, directory: Union[str, Path] = None, data_source: str = 'Phantom',
) -> Simulation:
    """Load Simulation.

    Parameters
    ----------
    prefix
        Simulation prefix, e.g. 'disc', if files are named like
        disc_00000.h5, disc01.ev, discSink0001N01.ev, etc.
    directory : optional
        Directory containing simulation snapshot files and auxiliary
        files. Default is None.
    data_source : optional
        The SPH code used to produce the simulation data. Default
        is 'Phantom'.
    """
    return Simulation().load_sim(
        prefix=prefix, directory=directory, data_source=data_source
    )
