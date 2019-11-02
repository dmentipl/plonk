"""The Particles class module.

Variables and functions related to the particle arrays.
"""

import numpy as np

I_GAS = 1
I_DUST = 7


class Arrays:
    """Smoothed particle hydrodynamics particle arrays object.

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

    >>> particles.arrays['xyz'][:]
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
        self._arrays_label = arrays_label

        _CAN_COMPUTE_DENSITY_LABELS = ('fluid', 'sph fluid', 'sph_fluid')
        self._can_compute_density = False
        if particle_type in _CAN_COMPUTE_DENSITY_LABELS:
            self._can_compute_density = True

        self._fields = None
        self._dtype = None
        self._shape = None
        self._number = None
        self._dimensions = None

        self._arrays = {field: self._get_array_handle(field) for field in self.fields}

        if cache_arrays is None:
            cache_arrays = False
        else:
            if not isinstance(cache_arrays, bool):
                raise TypeError('cache_array must be bool')
        if not cache_arrays:
            self._cache_arrays = False
        else:
            self.cache_arrays()

    def __repr__(self):
        """Dunder repr method."""
        return self.__str__()

    def __str__(self):
        """Dunder str method."""
        return f'<plonk.Arrays: "{self._arrays_label}">'

    @property
    def arrays(self):
        """Particle arrays."""
        if self._cache_arrays:
            return self._arrays_cached
        return self._arrays

    @property
    def number(self):
        """Number of particles."""
        if self._number:
            return self._number
        return self.arrays[self.fields[0]].size

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
            self._dtype = {key: item.dtype for key, item in self._arrays_handle.items()}
        return self._dtype

    @property
    def dimensions(self):
        """Dimension of quantities."""
        if self._dimensions is None:
            self._dimensions = {k: None for k in self.fields}
            for field in self.fields:
                if field in ('xyz', 'h', 'hsoft'):
                    self._dimensions[field] = 'L'
                elif field in ('vxyz'):
                    self._dimensions[field] = 'L T^-1'
                elif field in ('mass', 'm', 'maccreted'):
                    self._dimensions[field] = 'M'
                elif field in ('spinxyz'):
                    self._dimensions[field] = 'M L^2 T^-1'
                elif field in ('pressure'):
                    self._dimensions[field] = 'M L^-1 T^-1'
                elif field in ('dt', 'tstop', 'tlast'):
                    self._dimensions[field] = 'T'
                elif field in ('divv', 'curlv'):
                    self._dimensions[field] = 'T^-1'
        return self._dimensions

    @property
    def shape(self):
        """List of array shapes."""
        if self._shape is None:
            self._shape = {key: item.shape for key, item in self._arrays_handle.items()}
        return self._shape

    def to_structured_array(self):
        """Return arrays as Numpy structured array."""
        if self._cache_arrays:
            return self._arrays_cached
        return self._read_arrays()

    def _get_field_from_cache(self, field):
        return getattr(self, '_' + field)

    def _get_array_handle(self, array):
        """Get one array from file."""
        return self._arrays_handle[array]

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
