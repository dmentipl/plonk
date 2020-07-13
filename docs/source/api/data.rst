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
