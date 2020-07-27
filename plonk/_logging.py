"""Logging."""

import logging
import platform

# Logging
logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)

console_handler = logging.StreamHandler()
file_handler = logging.FileHandler('.plonk.log')

console_handler.setLevel(logging.INFO)
file_handler.setLevel(logging.DEBUG)

console_format = logging.Formatter('%(levelname)s - %(message)s')
file_format = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
console_handler.setFormatter(console_format)
file_handler.setFormatter(file_format)

logger.addHandler(console_handler)
logger.addHandler(file_handler)


def get_os_info():
    """Get the operating system version for logging."""
    system = platform.system()
    if system == 'Darwin':
        system = 'macOS'
    release = platform.release()
    return f'{system} version: {release}'


def logger_init(__version__):
    """Log versions and platform."""
    logger.debug(f'Plonk v{__version__} on Python {platform.python_version()}')
    logger.debug(f'{get_os_info()}, {platform.machine()}')
