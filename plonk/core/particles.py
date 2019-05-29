"""
This module contains the Particles class.

Variables and functions related to the particle arrays.
"""

import numpy as np

I_GAS = 1
I_DUST = 7


class Arrays:
    """
    Smoothed particle hydrodynamics particle arrays object.

    Used for accessing the particle arrays and sinks arrays from the
    dump file handle, i.e. with lazy loading.

    Parameters
    ----------
    arrays_label : str
        Label for the arrays to load corresponding to the group inside
        the file specified by file_handle, e.g. 'particles' or 'sinks'.

    file_handle : h5py File
        File handle to the dump file containing the arrays.

    Examples
    --------
    Creating the particle arrays object.

    >>> file_handle = h5py.File(filename)
    >>> particles = Arrays(file_handle=file_handle, arrays_label='particles')

    Accessing particle position array (two ways), available arrays, and
    array data types.

    >>> particles.xyz[()]
    >>> particles.arrays['xyz']
    >>> particles.arrays.fields
    >>> particles.arrays.dtype
    """

    def __init__(
        self,
        *,
        file_handle=None,
        arrays_label=None,
        particle_type=None,
        cache_arrays=None,
    ):

        self._arrays_handle = file_handle[arrays_label]

        _CAN_COMPUTE_DENSITY_LABELS = ['fluid', 'sph fluid', 'sph_fluid']
        self._can_compute_density = False
        if particle_type in _CAN_COMPUTE_DENSITY_LABELS:
            self._can_compute_density = True

        self._fields = None
        self._dtype = None
        self._shape = None

        self._mass = None
        self._density = None

        self._arrays = {
            field: self._get_array_handle(field) for field in self.fields
        }

        if cache_arrays is None:
            cache_arrays = False
        else:
            if not isinstance(cache_arrays, bool):
                raise TypeError('cache_array must be bool')
        if not cache_arrays:
            self._cache_arrays = False
            for field in self.fields:
                setattr(self, field, self._get_array_handle(field))
        else:
            self.cache_arrays()

    @property
    def arrays(self):
        if self._cache_arrays:
            return self._arrays_cached
        return self._arrays

    def _get_field_from_cache(self, field):
        return getattr(self, '_' + field)

    def _get_array_handle(self, array):
        """Get one array from file."""
        return self._arrays_handle[array]

    @property
    def mass(self):
        """Particle masses."""
        return self._mass

    @mass.setter
    def mass(self, value):
        self._mass = value

    @property
    def density(self):
        """Mass density."""
        return self._density

    @density.setter
    def density(self, value):
        self._density = value

    @property
    def fields(self):
        """List of available fields (array names)."""
        if self._fields is None:
            self._fields = tuple(self._arrays_handle)
        return self._fields

    @property
    def dtype(self):
        """List of array data types."""
        if self._dtype is None:
            self._dtype = tuple(
                [item.dtype for _, item in self._arrays_handle.items()]
            )
        return self._dtype

    @property
    def shape(self):
        """List of array shapes."""
        if self._shape is None:
            self._shape = tuple(
                [item.shape for _, item in self._arrays_handle.items()]
            )
        return self._shape

    def to_structured_array(self):
        """Return arrays as Numpy structured array."""
        if self._cache_arrays:
            return self._arrays_cached
        return self._read_arrays()

    def _read_arrays(self):
        """Read arrays into structured Numpy array."""

        dtype = []
        nvals = None
        for key, val in self._arrays_handle.items():
            if nvals is None:
                nvals = val.shape[0]
            if val.ndim == 1:
                dtype.append((key, val.dtype))
            elif val.ndim > 1:
                dtype.append((key, val.dtype, val.shape[1:]))

        struct_array = np.zeros(nvals, dtype=dtype)
        for key in struct_array.dtype.fields:
            struct_array[key] = self._arrays_handle[key][()]

        return struct_array

    def cache_arrays(self):
        """Load arrays into memory."""
        self._arrays_cached = self._read_arrays()
        for field in self.fields:
            setattr(self, field, self._arrays_cached[field])
        self._cache_arrays = True

    def extra_quantity(self, quantity, **kwargs):
        """
        Calculate extra quantity.

        Computes an extra quantity on the arrays specified by a string.
        For some quantities the particle masses are required. These may
        not be available by default.

        Parameters
        ----------
        quantity : str
            A string specifying the extra quantity to calculate.

        **kwargs
            Extra arguments to functions to calculate specific quantites.

        Examples
        --------
        Calculating the angular momentum two ways.

        >>> arrays = Arrays(label, handle)
        >>> arrays.extra_quantity('angular momentum')
        >>> arrays.extra_quantity('L')
        """

        quantities = [
            ('r', 'spherical radius'),
            ('R', 'cylindrical radius'),
            ('|v|', 'velocity magnitude'),
            ('p', 'momentum'),
            ('|p|', 'momentum magnitude'),
            ('L', 'angular momentum'),
            ('|L|', 'angular momentum magnitude'),
            ('l', 'specific angular momentum'),
            ('|l|', 'specific angular momentum magnitude'),
        ]

        require_mass = False

        if quantity not in [element for tupl in quantities for element in tupl]:
            print(f'{quantity} not available')
            return None

        if quantity in ['r', 'spherical radius']:
            data = (self.arrays['xyz'],)
            func = _spherical_radius

        elif quantity in ['R', 'cylindrical radius']:
            data = (self.arrays['xyz'],)
            func = _cylindrical_radius

        elif quantity in ['|v|', 'velocity magnitude']:
            require_mass = True
            data = (self.arrays['vxyz'],)
            func = _velocity_magnitude

        elif quantity in ['p', 'momentum']:
            require_mass = True
            data = (self.arrays['vxyz'], self.mass)
            func = _momentum

        elif quantity in ['|p|', 'momentum magnitude']:
            require_mass = True
            data = (self.arrays['vxyz'], self.mass)
            func = _momentum_magnitude

        elif quantity in ['L', 'angular momentum']:
            require_mass = True
            data = (self.arrays['xyz'], self.arrays['vxyz'], self.mass)
            func = _angular_momentum

        elif quantity in ['|L|', 'angular momentum magnitude']:
            require_mass = True
            data = (self.arrays['xyz'], self.arrays['vxyz'], self.mass)
            func = _angular_momentum

        elif quantity in ['l', 'specific angular momentum']:
            data = (self.arrays['xyz'], self.arrays['vxyz'])
            func = _specific_angular_momentum

        elif quantity in ['|l|', 'specific angular momentum magnitude']:
            data = self.arrays['xyz'], self.arrays['vxyz']
            func = _specific_angular_momentum_magnitude

        if self.mass is None and require_mass:
            raise ValueError('Particle masses required but not available')

        return _call_function_on_data(*data, func=func, **kwargs)


def _call_function_on_data(*data, func, **kwargs):
    return func(*data, **kwargs)


def _spherical_radius(position):
    return np.linalg.norm(position, axis=1)


def _cylindrical_radius(position):
    return np.linalg.norm(position[:, 0:2], axis=1)


def _cylindrical_azimuthal_angle(position):
    return np.arctan2(position[:, 1], position[:, 0])


def _velocity_magnitude(velocity):
    return np.linalg.norm(velocity, axis=1)


def _momentum(velocity, mass):
    if isinstance(mass, float) or isinstance(mass, int):
        return mass * velocity
    if isinstance(mass, np.ndarray):
        if mass.ndim == 0:
            return mass * velocity
        elif mass.ndim == 1:
            return mass[:, np.newaxis] * velocity
    raise ValueError('Check inputs, probably mass')


def _momentum_magnitude(velocity, mass):
    return mass * _velocity_magnitude(velocity)


def _specific_angular_momentum(position, velocity):
    if position.ndim > 1:
        if position.shape[1] != 3:
            raise ValueError('Wrong shape array')
    if position.shape != velocity.shape:
        raise ValueError('Position and velocity array shapes must match')
    return np.cross(position, velocity)


def _specific_angular_momentum_magnitude(position, velocity):
    return np.linalg.norm(
        _specific_angular_momentum(position, velocity), axis=1
    )


def _angular_momentum(position, velocity, mass):
    if isinstance(mass, float) or isinstance(mass, int):
        return mass * _specific_angular_momentum(position, velocity)
    if isinstance(mass, np.ndarray):
        if mass.ndim == 0:
            return mass * _specific_angular_momentum(position, velocity)
        elif mass.ndim == 1:
            return mass[:, np.newaxis] * _specific_angular_momentum(
                position, velocity
            )
    raise ValueError('Check inputs, probably mass')


def _angular_momentum_magnitude(position, velocity, mass):
    return np.linalg.norm(_angular_momentum(position, velocity, mass), axis=1)


def _eccentricity(position, velocity, gravitational_parameter):
    kinetic_energy = 1 / 2 * _velocity_magnitude(velocity) ** 2
    potential_energy = -gravitational_parameter / _spherical_radius(position)
    energy = kinetic_energy + potential_energy
    term = (
        2
        * energy
        * _specific_angular_momentum_magnitude(position, velocity) ** 2
        / gravitational_parameter ** 2
    )
    return np.sqrt(1 + term)
