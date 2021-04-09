"""Reader for snapshot files."""

# To write a new reader, three functions are required:
#
# - snap_properties_and_units,
# - snap_array_registry,
# - snap_sink_registry.
#
# These functions take an h5py.File object as an argument (and optionally
# a name map dictionary to rename arrays on file) and a string data_source
# to say where the data is from, e.g. 'phantom'.
#
# The first function gets some properties and code units, e.g. from the
# Phantom header. It is called by Snap.load_snap as follows
#
# self._properties, self._code_units = snap_properties_and_units(
#     file_pointer=self._file_pointer, data_source=self.data_source
# )
#
# The second two functions each create a dictionary with
# keys as names of arrays and values as a function that when called with a
# Snap object returns that array. They are called by Snap.load_snap as follows
#
# self._array_registry.update(
#     snap_array_registry(
#         file_pointer=self._file_pointer,
#         data_source=self.data_source,
#         name_map=self._name_map['particles'],
#     )
# )
#
# self._sink_registry.update(
#     snap_sink_registry(
#         file_pointer=self._file_pointer,
#         data_source=self.data_source,
#         name_map=self._name_map['sinks'],
#     )
# )
#
# See phantom.py for an example.

from __future__ import annotations

from typing import Any, Callable, Dict, Tuple

import h5py

from ..._units import Quantity
from .phantom import snap_array_registry as snap_array_registry_phantom
from .phantom import snap_properties_and_units as snap_properties_and_units_phantom
from .phantom import snap_sink_registry as snap_sink_registry_phantom

DATA_SOURCES = ['phantom']


def snap_properties_and_units(
    file_pointer: h5py.File, data_source: str
) -> Tuple[Dict[str, Any], Dict[str, Quantity]]:
    """Generate snap properties and units.

    Parameters
    ----------
    file_pointer
        The h5py file pointer to the snap file.
    data_source
        The SPH code that produced the data. E.g. 'phantom'.

    Returns
    -------
    properties
        The properties as a dict.
    code_units
        The code units as a dict.
    """
    # The properties dictionary has keys like "equation_of_state" and "time", with
    # varied values like "locally isothermal disc" and "61485663602.558136
    # <Unit('second')>".

    # The code_units dictionary has keys like "length" and "mass" with values as Pint
    # quantities, e.g. "149600000000.0 <Unit('meter')>" and "1.9891000000000002e+30
    # <Unit('kilogram')>".

    if data_source.lower() == 'phantom':
        return snap_properties_and_units_phantom(file_pointer=file_pointer)
    raise RuntimeError('Cannot generate properties and units')


def snap_array_registry(
    file_pointer: h5py.File,
    data_source: str,
    name_map: Dict[str, str] = None,
) -> Dict[str, Callable]:
    """Generate snap array registry.

    Parameters
    ----------
    file_pointer
        The h5py file pointer to the snap file.
    data_source
        The SPH code that produced the data. E.g. 'phantom'.
    name_map : optional
        A dict to convert from Phantom array names to Plonk names.

    Returns
    -------
    Dict
        The particle array registry.
    """
    # The keys are the names of the particle arrays and the values are functions that
    # return the array when called with snap as the argument.

    if data_source.lower() == 'phantom':
        return snap_array_registry_phantom(file_pointer=file_pointer, name_map=name_map)
    raise RuntimeError('Cannot generate array registry')


def snap_sink_registry(
    file_pointer: h5py.File,
    data_source: str,
    name_map: Dict[str, str] = None,
) -> Dict[str, Callable]:
    """Generate snap sink registry.

    Parameters
    ----------
    file_pointer
        The h5py file pointer to the snap file.
    data_source
        The SPH code that produced the data. E.g. 'phantom'.
    name_map : optional
        A dict to convert from Phantom array names to Plonk names.

    Returns
    -------
    Dict
        The sink array registry.
    """
    # The keys are the names of the sink arrays and the values are functions that return
    # the array when called with snap as the argument.

    if data_source.lower() == 'phantom':
        return snap_sink_registry_phantom(file_pointer=file_pointer, name_map=name_map)
    raise RuntimeError('Cannot generate sink registry')
