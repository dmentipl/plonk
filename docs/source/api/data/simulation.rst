~~~~~~~~~~
Simulation
~~~~~~~~~~

:py:class:`Simulation` is an aggregation of the :py:class:`Snap` and
pandas DataFrames to encapsulate all data within a single SPH simulation.

In addition, you can load auxilliary SPH simulation data, such as
globally-averaged time series data as pandas DataFrames.

.. autoclass:: plonk.Simulation
    :members:

.. autofunction:: plonk.load_sim

.. autofunction:: plonk.load_ev
