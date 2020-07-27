==========
Quickstart
==========

Quickstart guide for simple Plonk usage.

---------
Load data
---------

~~~~~~~~~~~~~
Load snapshot
~~~~~~~~~~~~~

Load a single snapshot, and access particle arrays and properties.

.. code-block:: python

    >>> import plonk
    >>> snap = plonk.load_snap('disc_00030.h5')
    >>> snap['position']
    array([[ -24.69953214,   49.60113417,   -4.98059478],
           [-108.99243136,   77.74663833,   12.89299546],
           [ -51.22218782,  108.64454019,    1.56619644],
           ...,
           [  93.296599  ,  -77.66042087,    5.40835798],
           [  63.75108128,   66.7446782 ,    3.30169363],
           [   8.11639008,  139.45117413,    7.55340187]])
    >>> snap.properties['time']
    61485663602.558136 <Unit('second')>


~~~~~~~~~~~~~~~~~~~~~
Load auxilliary files
~~~~~~~~~~~~~~~~~~~~~

Load a Phantom `.ev` file, and see what columns are available.

.. code-block:: python

    >>> import plonk
    >>> ev = plonk.load_ev('disc01.ev')
    >>> ev.columns
    Index(['time', 'energy_kinetic', 'energy_thermal', 'energy_magnetic',
           'energy_potential', 'energy_total', 'momentum', 'angular_momentum',
           'density_max', 'density_average', 'timestep', 'entropy',
           'mach_number_rms', 'velocity_rms', 'center_of_mass_x',
           'center_of_mass_y', 'center_of_mass_z', 'gas_density_max',
           'gas_density_average', 'dust_density_max', 'dust_density_average'],
          dtype='object')

~~~~~~~~~~~~~~~
Load simulation
~~~~~~~~~~~~~~~

Load a simulation, and access snapshots and other data.

.. code-block:: python

    >>> import plonk
    >>> sim = plonk.load_sim(prefix='disc')
    >>> snaps = sim.snaps
    >>> sim.global_quantities
    # Output is a pandas DataFrame
    >>> sim.sink_quantities
    # Output is a list of pandas DataFrames

-------------
Visualization
-------------

~~~~~~~~~~~~~~~
Projection plot
~~~~~~~~~~~~~~~

Produce a projection image plot of density.

.. code-block:: python

    >>> import plonk
    >>> snap = plonk.load_snap('disc_00030.h5')
    >>> plonk.visualize.plot(
    ...    snap=snap,
    ...    quantity='density',
    ...    extent=(-150, 150, -150, 150),
    ...    cmap='gist_heat',
    ... )

.. image:: _static/density2.png

--------
Analysis
--------

~~~~~~~
Profile
~~~~~~~

Create a radial profile.

.. code-block:: python

    >>> import matplotlib.pyplot as plt
    >>> import plonk
    >>> snap = plonk.load_snap('disc_00030.h5')
    >>> prof = plonk.load_profile(snap)
    >>> prof.available_profiles()
    ('angular_momentum_phi',
     'angular_momentum_theta',
     'aspect_ratio',
     'density',
     'dust_mass_001',
     'dust_surface_density_001',
     'dust_to_gas_ratio',
     'gas_mass',
     'gas_surface_density',
     'mass',
     'number',
     'position',
     'pressure',
     'radius',
     'scale_height',
     'size',
     'smoothing_length',
     'sound_speed',
     'stopping_time',
     'sub_type',
     'surface_density',
     'timestep',
     'toomre_Q',
     'type',
     'velocity',
     'velocity_divergence')
    >>> with plt.style.context('seaborn'):
    ...     ax = prof.plot('radius', 'scale_height')
    ...     ax.set_xlabel('Radius [au]')
    ...     ax.set_ylabel('Scale height [au]')
    >>> plt.show()

.. image:: _static/scale_height.png
