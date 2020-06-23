=================
API documentation
=================

.. currentmodule:: plonk

----
Snap
----

SPH snapshot files are represented by the :py:class:`Snap` class. This object
contains a properties dictionary, particle arrays, which are lazily loaded from
file.

.. autoclass:: plonk.Snap
.. autofunction:: plonk.load_snap

---------
Evolution
---------

Auxilliary SPH simulation data, such as globally-averaged quantities output more
frequently than snapshot files are represented by pandas DataFrames.

.. autofunction:: plonk.load_ev

----------
Simulation
----------

:py:class:`Simulation` is an aggregation of the :py:class:`Snap` and
pandas DataFrames to encapsulate all data within a single SPH simulation.

.. autoclass:: plonk.Simulation
.. autofunction:: plonk.load_sim

-------------
Visualization
-------------

The following functions are for Plonk SPH visualization.

.. autofunction:: plonk.visualize.plot
.. autofunction:: plonk.visualize.particle_plot
.. autofunction:: plonk.visualize.plot_snaps
.. autofunction:: plonk.visualize.animation
.. autofunction:: plonk.visualize.animation_profiles

~~~~~~~~~~~~~
Interpolation
~~~~~~~~~~~~~

Below are the functions for interpolation of particle data to a pixel grid.

.. autofunction:: plonk.visualize.interpolation.interpolate
.. autofunction:: plonk.visualize.interpolation.scalar_interpolation
.. autofunction:: plonk.visualize.interpolation.vector_interpolation

-------
Profile
-------

The :py:class:`Profile` class provides methods for generating radial profiles
from the particle data.

.. autoclass:: plonk.Profile
.. autoclass:: plonk.load_profile
