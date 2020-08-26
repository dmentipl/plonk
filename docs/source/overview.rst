========
Overview
========

.. currentmodule:: plonk

This document gives an overview of using Plonk for analysis and visualization of
smoothed particle hydrodynamics data. For a further guide see :doc:`usage`.

-----------------
Data file formats
-----------------

Plonk supports the following SPH file formats:

* Phantom output in
  `HDF <https://en.wikipedia.org/wiki/Hierarchical_Data_Format>`_
  form (as opposed to the sphNG-based Fortran binary format).

.. note::
    HDF5 output was added to Phantom as an option with git commit
    `9b22ded <https://github.com/danieljprice/phantom/commit/9b22ded9e7b4d512966f2b2e4b84d693b1afc9e6>`_
    on the 14th of March 2019. See the `Phantom documentation
    <https://phantomsph.readthedocs.io/>`_ for instructions on
    how to compile with HDF5 output and to convert from the sphNG-based
    output.

---------------------
Working with SPH data
---------------------

.. important::
    To follow along, download the sample data, :code:`plonk_example_data.tar`,
    from `Anaconda Cloud <https://anaconda.org/dmentipl/plonk_example_data/>`_.
    Then extract with :code:`tar xvf plonk_example_data.tar` and change into the
    :code:`plonk_example_data` directory. This data set is from a Phantom
    simulation of a dust and gas protoplanetary disc with an embedded
    protoplanet.

First import the Plonk package.

.. code-block:: pycon

    >>> import plonk

We also import Matplotlib and NumPy, for later.

.. code-block:: pycon

    >>> import matplotlib.pyplot as plt
    >>> import numpy as np

~~~~~~~~~
Snapshots
~~~~~~~~~

SPH snapshot files are represented by the :class:`Snap` class. This object
contains a properties dictionary, particle arrays, which are lazily loaded from
file. Here we demonstrate instantiating a :class:`Snap` object, and
accessing some properties and particle arrays.

First, we load the snapshot with the :func:`load_snap` function. You can
pass a string or :class:`pathlib.Path` object to point to the location of the
snapshot in the file system.

.. code-block:: pycon

    >>> filename = 'disc_00030.h5'
    >>> snap = plonk.load_snap(filename)

You can access arrays by their name passed in as a string.

.. code-block:: pycon

    >>> snap['position']
    array([[-3.69505001e+12,  7.42032967e+12, -7.45096980e+11],
           [-1.63052677e+13,  1.16308971e+13,  1.92879212e+12],
           [-7.66283930e+12,  1.62532232e+13,  2.34302988e+11],
           ...,
           [ 1.39571712e+13, -1.16179990e+13,  8.09090354e+11],
           [ 9.53716176e+12,  9.98500386e+12,  4.93933367e+11],
           [ 1.21421196e+12,  2.08618956e+13,  1.12998892e+12]]) <Unit('meter')>

There may be a small delay as the data is read from file. After the array is
read from file it is cached in memory, so that subsequent calls are faster.

To see what arrays are loaded into memory you can use the
:meth:`~Snap.loaded_arrays` method.

.. code-block:: pycon

    >>> snap.loaded_arrays()
    ['position']

Use :meth:`~Snap.available_arrays` to see what arrays are available. Some of
these arrays are stored on file, while others are computed as required.

.. code-block:: pycon

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

You can also define your own alias to access arrays. For example, if you prefer
to use the name `'coordinate'` rather than `'position',` use the
:meth:`~Snap.add_alias` method to add an alias.

.. code-block:: pycon

    >>> snap.add_alias(name='position', alias='coordinate')
    >>> snap['coordinate']
    array([[-3.69505001e+12,  7.42032967e+12, -7.45096980e+11],
           [-1.63052677e+13,  1.16308971e+13,  1.92879212e+12],
           [-7.66283930e+12,  1.62532232e+13,  2.34302988e+11],
           ...,
           [ 1.39571712e+13, -1.16179990e+13,  8.09090354e+11],
           [ 9.53716176e+12,  9.98500386e+12,  4.93933367e+11],
           [ 1.21421196e+12,  2.08618956e+13,  1.12998892e+12]]) <Unit('meter')>

The :class:`Snap` object has a :attr:`~Snap.properties` attribute which is a
dictionary of metadata, i.e. non-array data, on the snapshot.

.. code-block:: pycon

    >>> snap.properties['time']
    61485663602.558136 <Unit('second')>

    >>> list(snap.properties)
    ['adiabatic_index',
     'dust_method',
     'equation_of_state',
     'grain_density',
     'grain_size',
     'smoothing_length_factor',
     'time']

Units are available. We make use of the Python units library Pint. The code
units of the data are available as :attr:`~Snap.code_units`.

.. code-block:: pycon

    >>> snap.code_units['length']
    149600000000.0 <Unit('meter')>

You can set default units as follows.

.. code-block:: pycon

    >>> snap.set_units(position='au', density='g/cm^3', velocity='km/s')
    <plonk.Snap "disc_00030.h5">

    >>> snap['position']
    array([[ -24.6998837 ,   49.60184016,   -4.98066567],
           [-108.99398271,   77.74774493,   12.89317897],
           [ -51.22291689,  108.64608658,    1.56621873],
           ...,
           [  93.29792694,  -77.66152625,    5.40843496],
           [  63.75198868,   66.74562821,    3.30174062],
           [   8.11650561,  139.453159  ,    7.55350939]]) <Unit('astronomical_unit')>

Sink particles are handled separately from the fluid, e.g. gas or dust,
particles. They are available as an attribute :attr:`~Snap.sinks`.

.. code-block:: pycon

    >>> snap.sinks
    <plonk.snap sinks>

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


~~~~~~~~~~
Simulation
~~~~~~~~~~

SPH simulation data is usually spread over multiple files of, possibly,
different types, even though, logically, a simulation is a singular "object".
Plonk has the :class:`Simulation` class to represent the complete data set.
:class:`Simulation` is an aggregation of the :class:`Snap` and pandas DataFrames
to represent time series data (see below), plus metadata, such as the directory
on the file system.

Use the :func:`load_simulation` function to instantiate a :class:`Simulation`
object.

.. code-block:: pycon

    >>> prefix = 'disc'
    >>> sim = plonk.load_simulation(prefix=prefix)

Each of the snapshots are available via :attr:`~Simulation.snaps` as a list. We
can get the first five snapshots with the following.

.. code-block:: pycon

    >>> sim.snaps[:5]
    [<plonk.Snap "disc_00000.h5">,
     <plonk.Snap "disc_00001.h5">,
     <plonk.Snap "disc_00002.h5">,
     <plonk.Snap "disc_00003.h5">,
     <plonk.Snap "disc_00004.h5">]

The :class:`Simulation` class has an attribute :attr:`~Simulation.time_series`
which contains time series data as pandas DataFrames discussed in the next
section.

~~~~~~~~~~~
Time series
~~~~~~~~~~~

SPH simulation datasets often include auxiliary files containing
globally-averaged time series data output more frequently than snapshot files.
For example, Phantom writes text files with the file extension ".ev". These
files are output every time step rather than at the frequency of the snapshot
files.

We store this data in pandas DataFrames. Use :func:`load_time_series` to
instantiate.

.. code-block:: pycon

    >>> ts = plonk.load_time_series('disc01.ev')

The data may be split over several files, for example, if the simulation was run
with multiple jobs on a computation cluster. In that case, pass in a list of
files in chronological order to :func:`load_time_series`, and Plonk will
concatenate the data removing any duplicated time steps.

The underlying data is stored as a pandas [#f1]_ DataFrame. This allows for
the use of typical pandas operations with which users in the scientific Python
community may be familiar with.

.. code-block:: pycon

    >>> ts
                 time  energy_kinetic  energy_thermal  ...  gas_density_average  dust_density_max  dust_density_average
    0        0.000000        0.000013        0.001186  ...         8.231917e-10      1.720023e-10          8.015937e-12
    1        1.593943        0.000013        0.001186  ...         8.229311e-10      1.714059e-10          8.015771e-12
    2        6.375774        0.000013        0.001186  ...         8.193811e-10      1.696885e-10          8.018406e-12
    3       25.503096        0.000013        0.001186  ...         7.799164e-10      1.636469e-10          8.061417e-12
    4       51.006191        0.000013        0.001186  ...         7.249247e-10      1.580470e-10          8.210622e-12
    ..            ...             ...             ...  ...                  ...               ...                   ...
    548  12394.504462        0.000013        0.001186  ...         6.191121e-10      1.481833e-09          2.482929e-11
    549  12420.007557        0.000013        0.001186  ...         6.189791e-10      1.020596e-09          2.483358e-11
    550  12445.510653        0.000013        0.001186  ...         6.188052e-10      8.494835e-10          2.488946e-11
    551  12471.013748        0.000013        0.001186  ...         6.186160e-10      6.517475e-10          2.497029e-11
    552  12496.516844        0.000013        0.001186  ...         6.184558e-10      5.205011e-10          2.506445e-11

    [553 rows x 21 columns]

You can plot columns with the pandas plotting interface.

.. code-block:: pycon

    >>> ts.plot('time', ['center_of_mass_x', 'center_of_mass_y', 'center_of_mass_z'])

The previous code produces the following figure.

.. figure:: _static/ev.png

    The accretion disc center of mass as a function of time.

-------------------------
Visualization of SPH data
-------------------------

SPH particle data is not gridded like the data produced by, for example, finite
difference or finite volume hydrodynamical codes. One visualization method is to
plot the particles as a scatter plot, and possibly color the particles with the
magnitude of a quantity of interest. An alternative is to interpolate any
quantity on the particles to a pixel grid with weighted kernel density
estimation. This is what `Splash <https://github.com/danieljprice/splash>`_
does. For the technical details, see Price (2007), `PASA, 24, 3, 159
<https://ui.adsabs.harvard.edu/abs/2007PASA...24..159P>`_. We use the same
numerical method as Splash, with the Python function compiled with Numba so it
has the same performance as the Fortran code.

You can use the :meth:`~Snap.image` method to interpolate a quantity to a pixel
grid to show as an image. For example, in the following we produce a plot of
column density, i.e. a projection plot.

.. code-block:: pycon

    >>> filename = 'disc_00030.h5'
    >>> snap = plonk.load_snap(filename)

    >>> snap.image(quantity='density')

.. figure:: _static/density.png

    The total column density.

This produces an image via Matplotlib. The function returns a Matplotlib
:class:`AxesSubplot` object.

Alternatively, you can pass keyword arguments to the matplotlib functions. For
example, we set the units, the colormap to 'gist_heat' and set the colorbar
minimum and maxiumum. In addition, we set the extent, i.e. the x- and y-limits.

.. code-block:: pycon

    >>> snap.set_units(position='au', density='g/cm^3', projection='cm')
    >>> snap.image(
    ...     quantity='density',
    ...     extent=(20, 120, -50, 50),
    ...     cmap='gist_heat',
    ...     vmin=0.1,
    ...     vmax=0.2,
    ... )

.. figure:: _static/density_zoom.png

    The column density zoomed around the planet.

More fine-grained control can be achieved by using the full details of
:meth:`~Snap.image`. See the API for more details.

--------------------
Analysis of SPH data
--------------------

~~~~~~~~
Subsnaps
~~~~~~~~

When analyzing SPH data it can be useful to look at a subset of particles. For
example, the simulation we have been working with has dust and gas. So far we
have been plotting the total density. We may want to visualize the dust and gas
separately.

To do this we make a :class:`SubSnap` object. We can access these quantities
using the :meth:`~Snap.family` method. Given that there may be sub-types of
dust, using 'dust' returns a list by default. In this simulation there is only
one dust species. We can squeeze all the dust sub-types together using the
`squeeze` argument.

.. code-block:: pycon

    >>> gas = snap.family('gas')
    >>> dust = snap.family('dust', squeeze=True)

You can access arrays on the :py:class:`SubSnap` objects as for any
:py:class:`Snap` object.

.. code-block:: pycon

    >>> gas['mass'].sum().to('solar_mass')
    0.0010000000000000005 <Unit('solar_mass')>

    >>> dust['mass'].sum().to('earth_mass')
    3.326810503428664 <Unit('earth_mass')>

Let's plot the gas and dust side-by-side.

.. code-block:: pycon

    >>> subsnaps = [gas, dust]
    >>> extent = (-200, 200, -200, 200)

    >>> fig, axs = plt.subplots(ncols=2, figsize=(14, 5))

    >>> for subsnap, ax in zip(subsnaps, axs):
    ...     subsnap.image(quantity='density', extent=extent, cmap='gist_heat', ax=ax)

.. figure:: _static/dust-gas.png

    The column density of the gas and dust.

~~~~~~~~~~~~~~
Derived arrays
~~~~~~~~~~~~~~

Sometimes you need new arrays on the particles that are not available in the
snapshot files. Many are available in Plonk already. The following code block
lists the available raw Phantom arrays on the file.

.. code-block:: pycon

    >>> list(snap._file_pointer['particles'])
    ['divv', 'dt', 'dustfrac', 'h', 'itype', 'tstop', 'vxyz', 'xyz']

To see all available arrays on the :class:`Snap` object:

.. code-block:: pycon

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

This is a disc simulation. You can add quantities appropriate for discs with the
:meth:`~Snap.add_quantities` method.

.. code-block:: pycon

    >>> previous_arrays = snap.available_arrays()

    >>> snap.add_quantities('disc')

    # Additional available arrays
    >>> set(snap.available_arrays()) - set(previous_arrays)
    {'eccentricity',
     'inclination',
     'keplerian_frequency',
     'semi_major_axis',
     'stokes_number'}

You can create a new, derived array on the particles as follows.

.. code-block:: pycon

    >>> snap['rad'] = np.sqrt(snap['x'] ** 2 + snap['y'] ** 2)

    >>> snap['rad']
    array([8.28943225e+12, 2.00284678e+13, 1.79690392e+13, ...,
           1.81598604e+13, 1.38078875e+13, 2.08972008e+13]) <Unit('meter')>

Where, here, we have used the fact that Plonk knows that 'x' and 'y' refer to
the x- and y-components of the position array.

Alternatively, you can define a function for a derived array. This makes use of
the decorator :meth:`~Snap.add_array`.

.. code-block:: pycon

    >>> @snap.add_array()
    ... def radius(snap):
    ...     radius = np.hypot(snap['x'], snap['y'])
    ...     return radius

    >>> snap['radius']
    array([8.28943225e+12, 2.00284678e+13, 1.79690392e+13, ...,
           1.81598604e+13, 1.38078875e+13, 2.08972008e+13]) <Unit('meter')>

~~~~~
Units
~~~~~

Plonk uses Pint to set arrays to physical units.

.. code-block:: pycon

    >>> snap = plonk.load_snap(filename)

    >>> snap['x']
    array([-3.69505001e+12, -1.63052677e+13, -7.66283930e+12, ...,
            1.39571712e+13,  9.53716176e+12,  1.21421196e+12]) <Unit('meter')>

It is easy to convert quantities to different units as required.

.. code-block:: pycon

    >>> snap['x'].to('au')
    array([ -24.6998837 , -108.99398271,  -51.22291689, ...,   93.29792694,
             63.75198868,    8.11650561]) <Unit('astronomical_unit')>

~~~~~~~~
Profiles
~~~~~~~~

Generating a profile is a convenient method to reduce the dimensionality
of the full data set. For example, we may want to see how the surface density
and aspect ratio of the disc vary with radius.

To do this we use the :class:`Profile` class in the :mod:`analysis` module.

.. code-block:: pycon

    >>> snap = plonk.load_snap(filename)

    >>> snap.add_quantities('disc')

    >>> prof = plonk.load_profile(snap, cmin='10 au', cmax='200 au')

    >>> prof
    <plonk.Profile "disc_00030.h5">


To see what profiles are loaded and what are available use the
:meth:`~Profile.loaded_profiles` and :meth:`~Profile.available_profiles`
methods.

.. code-block:: pycon

    >>> prof.loaded_profiles()
    ['number', 'radius', 'size']

    >>> prof.available_profiles()
    ['alpha_shakura_sunyaev',
     'angular_momentum_mag',
     'angular_momentum_phi',
     'angular_momentum_theta',
     'angular_momentum_x',
     'angular_momentum_y',
     'angular_momentum_z',
     'angular_velocity',
     'aspect_ratio',
     'azimuthal_angle',
     'density',
     'disc_viscosity',
     'dust_to_gas_ratio_001',
     'eccentricity',
     'epicyclic_frequency',
     'id',
     'inclination',
     'keplerian_frequency',
     'kinetic_energy',
     'mass',
     'midplane_stokes_number_001',
     'momentum_mag',
     'momentum_x',
     'momentum_y',
     'momentum_z',
     'number',
     'polar_angle',
     'position_mag',
     'position_x',
     'position_y',
     'position_z',
     'pressure',
     'radius',
     'radius_cylindrical',
     'radius_spherical',
     'scale_height',
     'semi_major_axis',
     'size',
     'smoothing_length',
     'sound_speed',
     'specific_angular_momentum_mag',
     'specific_angular_momentum_x',
     'specific_angular_momentum_y',
     'specific_angular_momentum_z',
     'specific_kinetic_energy',
     'stokes_number_001',
     'stopping_time_001',
     'sub_type',
     'surface_density',
     'temperature',
     'timestep',
     'toomre_q',
     'type',
     'velocity_divergence',
     'velocity_mag',
     'velocity_radial_cylindrical',
     'velocity_radial_spherical',
     'velocity_x',
     'velocity_y',
     'velocity_z']

To load a profile, simply call it.

.. code-block:: pycon

    >>> prof['scale_height']
    array([7.97124951e+10, 9.84227660e+10, 1.18761140e+11, 1.37555034e+11,
           1.57492383e+11, 1.79386553e+11, 1.99666627e+11, 2.20898696e+11,
           2.46311866e+11, 2.63718852e+11, 2.91092881e+11, 3.11944125e+11,
           3.35101091e+11, 3.58479038e+11, 3.86137691e+11, 4.09731915e+11,
           4.34677739e+11, 4.61230912e+11, 4.87471032e+11, 5.07589421e+11,
           5.38769961e+11, 5.64068787e+11, 5.88487300e+11, 6.15197082e+11,
           6.35591229e+11, 6.75146081e+11, 6.94786320e+11, 7.32840731e+11,
           7.65750662e+11, 7.96047221e+11, 8.26764128e+11, 8.35658042e+11,
           8.57126522e+11, 8.99037935e+11, 9.26761324e+11, 9.45798955e+11,
           9.65997944e+11, 1.01555592e+12, 1.05554201e+12, 1.07641999e+12,
           1.10797835e+12, 1.13976869e+12, 1.17181502e+12, 1.20083661e+12,
           1.23779947e+12, 1.24785058e+12, 1.28800980e+12, 1.31290341e+12,
           1.33602542e+12, 1.34618607e+12, 1.36220344e+12, 1.39412340e+12,
           1.41778445e+12, 1.46713623e+12, 1.51140364e+12, 1.54496807e+12,
           1.58327631e+12, 1.60316504e+12, 1.63374960e+12, 1.64446331e+12,
           1.66063803e+12, 1.67890856e+12, 1.68505873e+12, 1.68507230e+12,
           1.71353612e+12, 1.71314330e+12, 1.75704484e+12, 1.79183025e+12,
           1.83696336e+12, 1.88823477e+12, 1.93080810e+12, 1.98301979e+12,
           2.05279086e+12, 2.11912539e+12, 2.14224572e+12, 2.21647741e+12,
           2.27153917e+12, 2.36605186e+12, 2.38922067e+12, 2.53901104e+12,
           2.61297334e+12, 2.64782574e+12, 2.73832897e+12, 3.04654121e+12,
           3.13575612e+12, 3.37636281e+12, 3.51482502e+12, 3.69591185e+12,
           3.88308614e+12, 3.83982313e+12, 4.00147149e+12, 4.37288049e+12,
           4.35330982e+12, 4.44686164e+12, 4.47133547e+12, 4.83307604e+12,
           4.63783507e+12, 4.95119779e+12, 5.17961431e+12, 5.29308491e+12]) <Unit('meter')>

You can convert the data in the :class:`Profile` object to a pandas DataFrame
with the :meth:`~Profile.to_dataframe` method. This takes all loaded profiles
and puts them into the DataFrame with units indicated in brackets.

.. code-block:: pycon

    >>> profiles = [
    ...    'radius',
    ...    'angular_momentum_phi',
    ...    'angular_momentum_theta',
    ...    'surface_density',
    ...    'scale_height',
    ... ]

    >>> df = prof.to_dataframe(columns=profiles)

    >>> df
          radius [m]  angular_momentum_phi [rad]  angular_momentum_theta [rad]  surface_density [kg / m ** 2]  scale_height [m]
    0   1.638097e+12                   -0.019731                      0.049709                       0.486007      7.971250e+10
    1   1.922333e+12                    1.914841                      0.053297                       0.630953      9.842277e+10
    2   2.206569e+12                    1.293811                      0.055986                       0.808415      1.187611e+11
    3   2.490805e+12                   -2.958286                      0.057931                       0.956107      1.375550e+11
    4   2.775041e+12                   -1.947547                      0.059679                       1.095939      1.574924e+11
    ..           ...                         ...                           ...                            ...               ...
    95  2.864051e+13                    3.045660                      0.168944                       0.013883      4.833076e+12
    96  2.892475e+13                   -0.054956                      0.161673                       0.013246      4.637835e+12
    97  2.920898e+13                   -0.217485                      0.169546                       0.011058      4.951198e+12
    98  2.949322e+13                   -1.305261                      0.175302                       0.011367      5.179614e+12
    99  2.977746e+13                    2.642077                      0.176867                       0.010660      5.293085e+12

    [100 rows x 5 columns]

We can also plot the profiles.

.. code-block:: pycon

    >>> prof.set_units(position='au', scale_height='au', surface_density='g/cm^2')

    >>> fig, axs = plt.subplots(ncols=2, figsize=(12, 5))

    >>> prof.plot('radius', 'surface_density', ax=axs[0])

    >>> prof.plot('radius', 'scale_height', ax=axs[1])

.. figure:: _static/profile.png

.. rubric:: Footnotes

.. [#f1] See `<https://pandas.pydata.org/>`_ for more on pandas.
