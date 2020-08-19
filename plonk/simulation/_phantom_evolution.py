"""Read Phantom evolution files."""

from __future__ import annotations

from pathlib import Path
from typing import TYPE_CHECKING, Dict, List, Tuple, Union

import numpy as np
import pandas as pd
from pandas import DataFrame

from .._config import load_config
from .._units import _convert_dim_string, _get_code_unit, generate_array_code_units

if TYPE_CHECKING:
    from .simulation import Simulation


def load_data_from_file(
    filenames: Union[str, Path, Tuple[str], Tuple[Path], List[str], List[Path]],
    config: Union[str, Path] = None,
):
    """Load data from Phantom .ev files.

    Parameters
    ----------
    filenames
        The filename or filenames (as a list).
    config : optional
        The path to a Plonk config.toml file.

    Returns
    -------
    DataFrame
    """
    if isinstance(filenames, (str, Path)):
        _filenames = [filenames]
    elif isinstance(filenames, (list, tuple)):
        _filenames = list(filenames)
    else:
        raise ValueError('filenames is not a known type')

    conf = load_config(filename=config)
    name_map = conf['phantom']['time_series']['namemap']

    _file_paths = list()
    for filename in _filenames:
        path = Path(filename)
        _file_paths.append(path.resolve())
    file_paths = tuple(_file_paths)

    _check_file_consistency(filenames=file_paths, name_map=name_map)
    columns = _get_columns(filename=file_paths[0], name_map=name_map)
    dataframe = _get_data(columns=columns, file_paths=file_paths)

    return dataframe


def evolution_units(sim: Simulation, config: Union[str, Path] = None):
    """Get units of Phantom .ev files from Simulation object.

    Parameters
    ----------
    sim
        The Simulation object.
    config : optional
        The path to a Plonk config.toml file.

    Returns
    -------
    Dict
    """
    conf = load_config(filename=config)
    arrays = conf['phantom']['time_series']['dimensions']
    dim = dict()
    for key, val in arrays.items():
        dim[key] = _convert_dim_string(val)
    _units = dict()
    for arr, unit in dim.items():
        _units[arr] = _get_code_unit(unit, sim.code_units)
    return _units


def _get_data(columns: Tuple[str, ...], file_paths: Tuple[Path, ...]) -> DataFrame:

    times = list()
    for filename in file_paths:
        times.append(np.loadtxt(filename, usecols=0))

    _skiprows = [0]
    if len(times) > 1:
        for t1, t2 in zip(times, times[1:]):
            _skiprows.append(np.where(t2 < t1[-1])[0][-1] + 2)

    df = pd.concat(
        (
            pd.read_csv(
                f,
                names=columns,
                skiprows=skiprows,
                skipinitialspace=True,
                delim_whitespace=True,
                comment='#',
            )
            for f, skiprows in zip(file_paths, _skiprows)
        )
    )

    df.reset_index(inplace=True, drop=True)
    return df


def _get_columns(filename: Path, name_map: Dict[str, str]) -> Tuple[str, ...]:

    with open(filename) as f:
        column_line = f.readline().strip('\n')

    _column_line = [item.strip('] ')[2:].strip(' ') for item in column_line.split('[')]
    columns = _column_line[1:]

    return tuple([name_map[col] if col in name_map else col for col in columns])


def _check_file_consistency(
    filenames: Tuple[Path, ...], name_map: Dict[str, str]
) -> None:

    columns = _get_columns(filenames[0], name_map)
    for filename in filenames:
        columns_previous = columns
        columns = _get_columns(filename, name_map)
        if columns != columns_previous:
            raise ValueError('files have different columns')
