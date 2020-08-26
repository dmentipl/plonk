"""Config."""

from pathlib import Path
from typing import Any, MutableMapping, Union

import toml

CONFIG_FILE = Path(__file__).parent / 'config.toml'


def read_config(filename: Union[str, Path] = None) -> MutableMapping[str, Any]:
    """Read config file.

    Parameters
    ----------
    filename
        The name of the config file (must be TOML) as a string or
        pathlib.Path.

    Returns
    -------
    Dict
        The config as a nested dictionary.
    """
    if filename is None:
        filename = CONFIG_FILE
    return toml.load(filename)


def write_config(filename: Union[str, Path]):
    """Write a default config file.

    Parameters
    ----------
    filename
        The name of the config file (must be TOML) as a string or
        pathlib.Path.
    """
    config = read_config()
    with open(filename, mode='w') as f:
        toml.dump(config, f)
