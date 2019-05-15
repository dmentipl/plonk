"""
simulation.py

Daniel Mentiplay, 2019.
"""

from pathlib import Path

from .dump import FILE_TYPES


class Simulation:
    """
    Smoothed particle hydrodynamics simulation data object. Aggregates
    dump files and auxiliary data.

    Parameters
    ----------
    prefix : str
        Simulation prefix, e.g. 'disc', if files are like disc_00000.h5,
        disc01.ev, etc.
    data_dir : str
        Directory containing simulation dump files and auxiliary files.

    Examples
    --------
    Reading simulation data into a Simulation object.

    >>> data_dir = 'simulation'
    >>> simulation = plonk.Simulation(data_dir)
    """

    # TODO: implement
    # TODO: add documentation

    def __init__(self, prefix, data_dir):

        self._prefix = prefix
        self._path = Path(data_dir).resolve()
        self._dump_file_type, self._dump_file_extension = (
            self._get_dump_file_type()
        )
        self._dump_files = self._get_dump_files()

    def _get_dump_files(self):
        """Get dump files."""

        return list(
            self._path.glob(
                self._prefix + '_[0-9]*.' + self._dump_file_extension
            )
        )

    def _get_dump_file_type(self):
        """
        Determine dump file type assuming file names like 'prefix_123.ext'.
        """

        file_types = set(
            [f.suffix for f in self._path.glob(self._prefix + '_[0-9]*.*')]
        )

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

        for ft in FILE_TYPES:
            if file_ext == ft.extension:
                return ft.filetype, ft.extension
