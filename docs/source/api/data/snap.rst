~~~~
Snap
~~~~

SPH snapshot files are represented by the :py:class:`Snap` class. This object
contains a properties dictionary, particle arrays, which are lazily loaded from
file.

.. autoclass:: plonk.Snap
    :members:

.. autoclass:: plonk.SubSnap
    :members:

.. autoclass:: plonk.Sinks
    :members:

.. autofunction:: plonk.load_snap
.. autofunction:: plonk.snap.gravitational_constant_in_code_units
.. autofunction:: plonk.snap.get_array_in_code_units
