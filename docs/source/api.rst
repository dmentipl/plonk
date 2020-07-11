=================
API documentation
=================

.. currentmodule:: plonk

----
Data
----

~~~~
Snap
~~~~

SPH snapshot files are represented by the :py:class:`Snap` class. This object
contains a properties dictionary, particle arrays, which are lazily loaded from
file.

.. autoclass:: plonk.Snap
.. autoclass:: plonk.snap.SubSnap
.. autofunction:: plonk.load_snap
.. autofunction:: plonk.snap.gravitational_constant_in_code_units

~~~~~~~~~
Evolution
~~~~~~~~~

Auxilliary SPH simulation data, such as globally-averaged quantities output more
frequently than snapshot files are represented by pandas DataFrames.

.. autofunction:: plonk.load_ev

~~~~~~~~~~
Simulation
~~~~~~~~~~

:py:class:`Simulation` is an aggregation of the :py:class:`Snap` and
pandas DataFrames to encapsulate all data within a single SPH simulation.

.. autoclass:: plonk.Simulation
.. autofunction:: plonk.load_sim

--------
Analysis
--------

.. automodule:: plonk.analysis

~~~~~~~
Profile
~~~~~~~

The :py:class:`Profile` class provides methods for generating radial profiles
from the particle data.

.. autoclass:: plonk.Profile
.. autofunction:: plonk.load_profile

~~~~~~~~~~~~~~~~~~~
Particle quantities
~~~~~~~~~~~~~~~~~~~

.. automodule:: plonk.analysis.particles

~~~~~~~~~~~~~~~
Sink quantities
~~~~~~~~~~~~~~~

.. automodule:: plonk.analysis.sinks

~~~~~~~~~~~~~~~~~
Global quantities
~~~~~~~~~~~~~~~~~

.. automodule:: plonk.analysis.total

-------------
Visualization
-------------

~~~~~
Plots
~~~~~

The following functions are for Plonk SPH visualization.

.. autofunction:: plonk.visualize.plot
.. autofunction:: plonk.visualize.particle_plot
.. autofunction:: plonk.visualize.plot_snaps

~~~~~~~~~
Animation
~~~~~~~~~

The following functions are for Plonk animations.

.. autofunction:: plonk.visualize.animation
.. autofunction:: plonk.visualize.animation_profiles
.. autofunction:: plonk.visualize.animation_particles

~~~~~~~~~~~~~
Interpolation
~~~~~~~~~~~~~

Below are the functions for interpolation of particle data to a pixel grid.

.. autofunction:: plonk.visualize.interpolation.interpolate
.. autofunction:: plonk.visualize.interpolation.scalar_interpolation
.. autofunction:: plonk.visualize.interpolation.vector_interpolation

~~~~~~~
Utility
~~~~~~~

The following are utility functions for visualization.

.. autofunction:: plonk.visualize.get_extent_from_percentile
.. autofunction:: plonk.visualize.str_to_units

-------
Utility
-------

Here are some useful utility functions.

.. autofunction:: plonk.utils.average
.. autofunction:: plonk.utils.cartesian_to_polar
.. autofunction:: plonk.utils.cross
.. autofunction:: plonk.utils.is_documented_by
.. autofunction:: plonk.utils.norm
.. autofunction:: plonk.utils.sph
.. autofunction:: plonk.utils.time_string
