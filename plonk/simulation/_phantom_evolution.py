"""Read Phantom evolution files."""

from pathlib import Path
from typing import Dict, List, Tuple, Union

import numpy as np
import pandas as pd
from pandas import DataFrame

NAME_MAP = {
    'time': 'time',
    'ekin': 'energy_kinetic',
    'etherm': 'energy_thermal',
    'emag': 'energy_magnetic',
    'epot': 'energy_potential',
    'etot': 'energy_total',
    'totmom': 'momentum',
    'angtot': 'angular_momentum',
    'rho max': 'density_max',
    'rho ave': 'density_average',
    'dt': 'timestep',
    'dtmax': 'timestep_max',
    'totentrop': 'entropy',
    'rmsmach': 'mach_number_rms',
    'vrms': 'velocity_rms',
    'xcom': 'center_of_mass_x',
    'ycom': 'center_of_mass_y',
    'zcom': 'center_of_mass_z',
    'rho gas max': 'gas_density_max',
    'rho gas ave': 'gas_density_average',
    'rho dust X': 'dust_density_max',
    'rho dust A': 'dust_density_average',
    'rho bdy max': 'boundary_density_max',
    'rho bdy ave': 'boundary_density_average',
    'rho star X': 'star_density_max',
    'rho star A': 'star_density_average',
    'rho dm max': 'dark_matter_density_max',
    'rho dm ave': 'dark_matter_density_average',
    'rho blg max': 'bulge_density_max',
    'rho blg ave': 'bulge_density_average',
    'alpha': 'alpha_viscosity_numerical',
    'B max': 'magnetic_field_max',
    'B min': 'magnetic_field_min',
    'B ave': 'magnetic_field_average',
    'divB max': 'magnetic_field_divergence_max',
    'divB ave': 'magnetic_field_divergence_average',
    'beta_P max': 'plasma_beta_max',
    'beta_P min': 'plasma_beta_min',
    'beta_P ave': 'plasma_beta_average',
    'erot_x': 'energy_rotational_x',
    'erot_y': 'energy_rotational_y',
    'erot_z': 'energy_rotational_z',
    'erot': 'energy_rotational_total',
    'dust/gas': 'dust_gas_ratio',
    't_s': 'stopping_time',
    'x': 'position_x',
    'y': 'position_y',
    'z': 'position_z',
    'vx': 'velocity_x',
    'vy': 'velocity_y',
    'vz': 'velocity_z',
    'spinx': 'spin_x',
    'spiny': 'spin_y',
    'spinz': 'spin_z',
    'macc': 'mass_accreted',
    'fx': 'force_x',
    'fy': 'force_y',
    'fz': 'force_z',
    'fssx': 'sink_sink_force_x',
    'fssy': 'sink_sink_force_y',
    'fssz': 'sink_sink_force_z',
}

NAME_MAP.update({f'DustMass{idx:03}': f'dust_mass_{idx:03}' for idx in range(100)})


def load_data_from_file(
    filenames: Union[str, Path, Tuple[str], Tuple[Path], List[str], List[Path]],
):
    """Load data from Phantom .ev files."""
    if isinstance(filenames, (str, Path)):
        _filenames = [filenames]
    elif isinstance(filenames, (list, tuple)):
        _filenames = list(filenames)
    else:
        raise ValueError('filenames is not a known type')

    _file_paths = list()
    for filename in _filenames:
        path = Path(filename)
        _file_paths.append(path.resolve())
    file_paths = tuple(_file_paths)

    _check_file_consistency(file_paths, NAME_MAP)
    columns = _get_columns(file_paths[0], NAME_MAP)
    dataframe = _get_data(columns, file_paths)

    return dataframe


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
