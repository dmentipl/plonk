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
    12241.485887903227

~~~~~~~~~~~~~~~~~~~~~
Load auxilliary files
~~~~~~~~~~~~~~~~~~~~~

Load a Phantom `.ev` file, and see what columns are available.

.. code-block:: python

    >>> import plonk
    >>> ev = plonk.load_ev('disc01.ev')
    >>> ev.columns
    ('time',
     'ekin',
     'etherm',
     'emag',
     'epot',
     'etot',
     'totmom',
     'angtot',
     'rho max',
     'rho ave',
     'dt',
     'totentrop',
     'rmsmach',
     'vrms',
     'xcom',
     'ycom',
     'zcom',
     'rho gas max',
     'rho gas ave',
     'rho dust X',
     'rho dust A')

~~~~~~~~~~~~~~~
Load simulation
~~~~~~~~~~~~~~~

Load a simulation, and access snapshots and other data.

.. code-block:: python

    >>> import plonk
    >>> sim = plonk.load_sim(prefix='disc')
    >>> snaps = sim.snaps
    >>> sim.global_quantities
    <plonk.Evolution: "('disc01.ev',)">
    >>> sim.sink_quantities
    [<plonk.Evolution: "('discSink0001N01.ev',)">,
     <plonk.Evolution: "('discSink0002N01.ev',)">]

-------------
Visualization
-------------

~~~~~~~~~~~~~~~
Projection plot
~~~~~~~~~~~~~~~

Produce a projection rendering of density.

.. code-block:: python

    >>> import plonk
    >>> snap = plonk.load_snap('disc_00030.h5')
    >>> viz = plonk.visualize.plot(
    ...    snap=snap,
    ...    quantity='density',
    ...    extent=(-150, 150, -150, 150),
    ...    cmap='gist_heat',
    ... )

.. image:: _static/density.png

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
    >>> prof = plonk.Profile(snap)
    >>> prof.available_keys()
    ('angmom_mag',
     'angmom_phi',
     'angmom_theta',
     'density',
     'eccentricity',
     'mass',
     'number',
     'radius',
     'scale_height',
     'smoothing_length')
    >>> with plt.style.context('seaborn'):
    ...     fig, ax = prof.plot('radius', 'scale_height')
    ...     ax.set_xlabel('Radius [au]')
    ...     ax.set_ylabel('Scale height [au]')
    >>> plt.show()

.. image:: _static/scale_height.png
