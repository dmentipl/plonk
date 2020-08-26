----
Data
----

.. currentmodule:: plonk

~~~~~~~~~~~~~
Load snapshot
~~~~~~~~~~~~~

Load data from a single snapshot into :class:`Snap` and see what arrays are
available with :meth:`~Snap.available_arrays`.

.. code-block:: python

    >>> import plonk

    >>> snap = plonk.load_snap('disc_00030.h5')

    >>> snap.available_arrays()
    ['angular_momentum',
     'angular_velocity',
     'azimuthal_angle',
     'density',
     'dust_to_gas_ratio',
     'id',
     'kinetic_energy',
     'mass',
     'momentum',
     'polar_angle',
     'position',
     'pressure',
     'radius_cylindrical',
     'radius_spherical',
     'smoothing_length',
     'sound_speed',
     'specific_angular_momentum',
     'specific_kinetic_energy',
     'stopping_time',
     'sub_type',
     'temperature',
     'timestep',
     'type',
     'velocity',
     'velocity_divergence',
     'velocity_radial_cylindrical',
     'velocity_radial_spherical']

Load a single snapshot and access particle arrays and :attr:`~Snap.properties`.

.. code-block:: python

    >>> import plonk

    >>> snap = plonk.load_snap('disc_00030.h5')

    >>> snap['position']
    array([[-3.69505001e+12,  7.42032967e+12, -7.45096980e+11],
           [-1.63052677e+13,  1.16308971e+13,  1.92879212e+12],
           [-7.66283930e+12,  1.62532232e+13,  2.34302988e+11],
           ...,
           [ 1.39571712e+13, -1.16179990e+13,  8.09090354e+11],
           [ 9.53716176e+12,  9.98500386e+12,  4.93933367e+11],
           [ 1.21421196e+12,  2.08618956e+13,  1.12998892e+12]]) <Unit('meter')>

    >>> snap['position'].to('au')
    array([[ -24.6998837 ,   49.60184016,   -4.98066567],
           [-108.99398271,   77.74774493,   12.89317897],
           [ -51.22291689,  108.64608658,    1.56621873],
           ...,
           [  93.29792694,  -77.66152625,    5.40843496],
           [  63.75198868,   66.74562821,    3.30174062],
           [   8.11650561,  139.453159  ,    7.55350939]]) <Unit('astronomical_unit')>

    >>> snap.properties['time']
    61485663602.558136 <Unit('second')>

Load a single snapshot and access sink arrays via :attr:`~Snap.sinks` attribute.

.. code-block:: python

    >>> import plonk

    >>> snap = plonk.load_snap('disc_00030.h5')

    >>> sinks = snap.sinks

    >>> sinks.available_arrays()
    ['accretion_radius',
     'last_injection_time',
     'mass',
     'mass_accreted',
     'position',
     'softening_radius',
     'spin',
     'velocity']

    >>> sinks['spin']
    array([[ 3.56866999e+36, -1.17910663e+37,  2.44598074e+40],
           [ 4.14083556e+36,  1.19118555e+36,  2.62569386e+39]]) <Unit('kilogram * meter ** 2 / second')>


~~~~~~~~~~~~~~~~~~~~~
Load auxilliary files
~~~~~~~~~~~~~~~~~~~~~

Load a Phantom `.ev` file with :func:`load_ev` and see what columns are
available.

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

    >>> ev
                 time  energy_kinetic  ...  dust_density_max  dust_density_average
    0        0.000000        0.000013  ...      1.720023e-10          8.015937e-12
    1        1.593943        0.000013  ...      1.714059e-10          8.015771e-12
    2        6.375774        0.000013  ...      1.696885e-10          8.018406e-12
    3       25.503096        0.000013  ...      1.636469e-10          8.061417e-12
    4       51.006191        0.000013  ...      1.580470e-10          8.210622e-12
    ..            ...             ...  ...               ...                   ...
    548  12394.504462        0.000013  ...      1.481833e-09          2.482929e-11
    549  12420.007557        0.000013  ...      1.020596e-09          2.483358e-11
    550  12445.510653        0.000013  ...      8.494835e-10          2.488946e-11
    551  12471.013748        0.000013  ...      6.517475e-10          2.497029e-11
    552  12496.516844        0.000013  ...      5.205011e-10          2.506445e-11

    [553 rows x 21 columns]

~~~~~~~~~~~~~~~
Load simulation
~~~~~~~~~~~~~~~

Load a simulation with :func:`load_sim` and access snapshots and time series
data with :attr:`~Simulation.snaps` and :attr:`~Simulation.time_series`
attributes.

.. code-block:: python

    >>> import plonk

    >>> sim = plonk.load_sim(prefix='disc')

    >>> sim.snaps
    [<plonk.Snap "disc_00000.h5">,
     <plonk.Snap "disc_00001.h5">,
     <plonk.Snap "disc_00002.h5">,
     <plonk.Snap "disc_00003.h5">,
     <plonk.Snap "disc_00004.h5">,
     <plonk.Snap "disc_00005.h5">,
     <plonk.Snap "disc_00006.h5">,
     <plonk.Snap "disc_00007.h5">,
     <plonk.Snap "disc_00008.h5">,
     <plonk.Snap "disc_00009.h5">,
     <plonk.Snap "disc_00010.h5">,
     <plonk.Snap "disc_00011.h5">,
     <plonk.Snap "disc_00012.h5">,
     <plonk.Snap "disc_00013.h5">,
     <plonk.Snap "disc_00014.h5">,
     <plonk.Snap "disc_00015.h5">,
     <plonk.Snap "disc_00016.h5">,
     <plonk.Snap "disc_00017.h5">,
     <plonk.Snap "disc_00018.h5">,
     <plonk.Snap "disc_00019.h5">,
     <plonk.Snap "disc_00020.h5">,
     <plonk.Snap "disc_00021.h5">,
     <plonk.Snap "disc_00022.h5">,
     <plonk.Snap "disc_00023.h5">,
     <plonk.Snap "disc_00024.h5">,
     <plonk.Snap "disc_00025.h5">,
     <plonk.Snap "disc_00026.h5">,
     <plonk.Snap "disc_00027.h5">,
     <plonk.Snap "disc_00028.h5">,
     <plonk.Snap "disc_00029.h5">,
     <plonk.Snap "disc_00030.h5">]

    >>> sim.time_series['global']
             time [s]  ...  dust_density_average [kg / m ** 3]
    0    0.000000e+00  ...                        4.762293e-15
    1    8.005946e+06  ...                        4.762195e-15
    2    3.202378e+07  ...                        4.763760e-15
    3    1.280951e+08  ...                        4.789313e-15
    4    2.561903e+08  ...                        4.877956e-15
    ..            ...  ...                                 ...
    548  6.225423e+10  ...                        1.475116e-14
    549  6.238233e+10  ...                        1.475371e-14
    550  6.251042e+10  ...                        1.478691e-14
    551  6.263852e+10  ...                        1.483493e-14
    552  6.276661e+10  ...                        1.489087e-14

    [553 rows x 21 columns]

    >>> sim.time_series['sinks']
    [          time [s]  position_x [m]  ...  sink_sink_force_y [N]  sink_sink_force_z [N]
     0     4.002973e+06   -1.068586e+10  ...           2.455750e+18           2.020208e+14
     1     8.005946e+06   -1.068584e+10  ...           4.911471e+18           2.020271e+14
     2     1.601189e+07   -1.068574e+10  ...           9.822886e+18           2.020559e+14
     3     3.202378e+07   -1.068536e+10  ...           1.964551e+19           2.021794e+14
     4     6.404757e+07   -1.068380e+10  ...           3.928915e+19           1.407952e+14
     ...            ...             ...  ...                    ...                    ...
     1038  6.257447e+10   -9.976304e+09  ...           7.513519e+20           7.902581e+15
     1039  6.263852e+10   -9.895056e+09  ...           7.884134e+20           7.664701e+15
     1040  6.270257e+10   -9.809975e+09  ...           8.251714e+20           7.628553e+15
     1041  6.276661e+10   -9.721096e+09  ...           8.616085e+20           7.705555e+15
     1042  6.283066e+10   -9.628452e+09  ...           8.977201e+20           7.391809e+15

     [1043 rows x 18 columns],
               time [s]  position_x [m]  ...  sink_sink_force_y [N]  sink_sink_force_z [N]
     0     4.002973e+06    1.120931e+13  ...          -2.573360e+21          -2.116959e+17
     1     8.005946e+06    1.120928e+13  ...          -5.146689e+21          -2.117025e+17
     2     1.601189e+07    1.120918e+13  ...          -1.029332e+22          -2.117327e+17
     3     3.202378e+07    1.120877e+13  ...          -2.058637e+22          -2.118621e+17
     4     6.404757e+07    1.120715e+13  ...          -4.117069e+22          -1.475378e+17
     ...            ...             ...  ...                    ...                    ...
     1038  6.257447e+10    1.041137e+13  ...          -7.748986e+23          -8.150240e+18
     1039  6.263852e+10    1.032805e+13  ...          -8.131071e+23          -7.904764e+18
     1040  6.270257e+10    1.024073e+13  ...          -8.510027e+23          -7.867359e+18
     1041  6.276661e+10    1.014946e+13  ...          -8.885717e+23          -7.946693e+18
     1042  6.283066e+10    1.005427e+13  ...          -9.257998e+23          -7.623017e+18

     [1043 rows x 18 columns]]
