========
Overview
========

This document gives an overview of using Plonk for analysis and visualization of
smoothed particle hydrodynamics data. For a quickstart guide see
:doc:`quickstart`.

-----------------
Data file formats
-----------------

Plonk supports the following SPH file formats:

* Phantom output in
  `HDF <https://en.wikipedia.org/wiki/Hierarchical_Data_Format>`_
  form (as opposed to the sphNG-based Fortran binary format).

.. note::
    HDF5 output was added to Phantom as an option with git commit
    `9b22ded <https://bitbucket.org/danielprice/phantom/commits/9b22ded9e7b4d512966f2b2e4b84d693b1afc9e6>`_
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

SPH snapshot files are represented by the :py:class:`Snap` class. This object
contains a properties dictionary, particle arrays, which are lazily loaded from
file. Here we demonstrate instantiating a :py:class:`Snap` object, and accessing
some properties and particle arrays.

First, we load the snapshot with the :py:func:`load_snap` function. You can
pass a string or :class:`pathlib.Path` object to point to the location of the
snapshot in the file system.

.. code-block:: pycon

    >>> filename = 'disc_00030.h5'
    >>> snap = plonk.load_snap(filename)

You can access arrays by their name passed in as a string.

.. code-block:: pycon

    >>> snap['position']
    array([[ -24.69953214,   49.60113417,   -4.98059478],
           [-108.99243136,   77.74663833,   12.89299546],
           [ -51.22218782,  108.64454019,    1.56619644],
           ...,
           [  93.296599  ,  -77.66042087,    5.40835798],
           [  63.75108128,   66.7446782 ,    3.30169363],
           [   8.11639008,  139.45117413,    7.55340187]])

There may be a small delay as the data is read from file. After the array is
read from file it is cached in memory, so that subsequent calls are faster.

To see what arrays are loaded into memory you can use the
:py:meth:`loaded_arrays` method.

.. code-block:: pycon

    >>> snap.loaded_arrays()
    ('position',)

Use :py:meth:`available_arrays` to see what arrays are available. Some of these
arrays are stored on file, while others are computed as required from functions
defined in the analysis module.

.. code-block:: pycon

    >>> snap.available_arrays()
    ('angular_momentum',
     'angular_velocity',
     'azimuthal_angle',
     'density',
     'dust_density',
     'dust_fraction',
     'dust_mass',
     'dust_type',
     'eccentricity',
     'gas_density',
     'gas_mass',
     'inclination',
     'kinetic_energy',
     'mass',
     'momentum',
     'polar_angle',
     'position',
     'pressure',
     'radial_distance',
     'radial_velocity',
     'semi_major_axis',
     'smoothing_length',
     'sound_speed',
     'specific_angular_momentum',
     'specific_kinetic_energy',
     'stopping_time',
     'timestep',
     'type',
     'velocity',
     'velocity_divergence')

You can also define your own alias to access arrays. For example, if you prefer
to use the name `'coordinate'` rather than `'position',` use the
:py:meth:`add_alias` method to add an alias.

.. code-block:: pycon

    >>> snap.add_alias(name='position', alias='coordinate')
    >>> snap['coordinate']
    array([[ -24.69953214,   49.60113417,   -4.98059478],
           [-108.99243136,   77.74663833,   12.89299546],
           [ -51.22218782,  108.64454019,    1.56619644],
           ...,
           [  93.296599  ,  -77.66042087,    5.40835798],
           [  63.75108128,   66.7446782 ,    3.30169363],
           [   8.11639008,  139.45117413,    7.55340187]])

The :py:class:`Snap` object has a :py:attr:`properties` attribute which is a
dictionary of metadata, i.e. non-array data, on the snapshot.

.. code-block:: pycon

    >>> snap.properties['time']
    61485663602.558136 <Unit('second')>


    >>> list(snap.properties)
    ['time',
     'smoothing_length_factor',
     'adiabatic_index',
     'polytropic_constant',
     'equation_of_state',
     'sound_speed_index',
     'dust_method',
     'grain_size',
     'grain_density']

Units are available vis the :py:attr:`units` attribute. We make use of the
Python units library Pint.

.. code-block:: pycon

    >>> snap.units['length']
    14960000000000.0 <Unit('centimeter')>

Sink particles are handled separately from the fluid, e.g. gas or dust,
particles. They are available as an attribute.

.. code-block:: pycon

    >>> snap.available_arrays(sinks=True)
    ('accretion_radius',
     'last_injection_time',
     'mass',
     'mass_accreted',
     'position',
     'softening_radius',
     'spin',
     'velocity')

    >>> snap.sinks
    <plonk.snap sinks>

    >>> snap.sinks['spin']
    array([[ 4.02648711e-10, -1.33037173e-09,  2.75977043e-06],
           [ 4.67205458e-10,  1.34400022e-10,  2.96253856e-07]])

~~~~~~~~~~
Simulation
~~~~~~~~~~

SPH simulation data is usually spread over multiple files of, possibly,
different types, even though, logically, a simulation is a singular "object".
Plonk has the :py:class:`Simulation` class to represent the complete data set.
:py:class:`Simulation` is an aggregation of the :py:class:`Snap` and
pandas DataFrames to represent time evolution data (see below) objects, plus
metadata, such as the directory on the file system.

Use the :py:func:`load_sim` function to instantiate a :py:class:`Simulation`
object.

.. code-block:: pycon

    >>> prefix = 'disc'
    >>> sim = plonk.load_sim(prefix=prefix)

Each of the snapshots are available via :py:attr:`snaps` as a list. We can get
the first five snapshots with the following.

.. code-block:: pycon

    >>> sim.snaps[:5]
    [<plonk.Snap>, <plonk.Snap>, <plonk.Snap>, <plonk.Snap>, <plonk.Snap>]

The :py:class:`Simulation` class has attributes :py:attr:`global_quantities` and
:py:attr:`sink_quantities` which are pandas DataFrames discussed in the next
section.

~~~~~~~~~
Evolution
~~~~~~~~~

SPH simulation data also include auxiliary files containing globally-averaged
quantities output more frequently than snapshot files. For example, Phantom
writes text files with the suffix :code:`.ev`. These files are output every time
step rather than at the frequency of the snapshot files.

We store this data in pandas DataFrames. Use :py:meth:`load_ev` to instantiate.

.. code-block:: pycon

    >>> ev = plonk.load_ev('disc01.ev')

The data may be split over several files, for example, if the simulation was run
with multiple jobs on a computation cluster. In that case, pass in a tuple or
list of files in chronological order to :py:func:`load_ev`, and Plonk will
concatenate the data removing any duplicated time steps.

The underlying data is stored as a pandas [#f1]_ DataFrame. This allows for
the use of typical pandas operations with which users in the scientific Python
community may be familiar with.

.. code-block:: pycon

    >>> ev
                time      ekin    etherm  emag  ...   rho gas max   rho gas ave    rho dust X    rho dust A
    0        0.000000  0.000013  0.001186   0.0  ...  1.613665e-08  8.231917e-10  1.720023e-10  8.015937e-12
    1        1.593943  0.000013  0.001186   0.0  ...  1.599245e-08  8.229311e-10  1.714059e-10  8.015771e-12
    2        6.375774  0.000013  0.001186   0.0  ...  1.431020e-08  8.193811e-10  1.696885e-10  8.018406e-12
    3       25.503096  0.000013  0.001186   0.0  ...  7.422843e-09  7.799164e-10  1.636469e-10  8.061417e-12
    4       51.006191  0.000013  0.001186   0.0  ...  5.898794e-09  7.249247e-10  1.580470e-10  8.210622e-12
    ..            ...       ...       ...   ...  ...           ...           ...           ...           ...
    548  12394.504462  0.000013  0.001186   0.0  ...  4.768277e-09  6.191121e-10  1.481833e-09  2.482929e-11
    549  12420.007557  0.000013  0.001186   0.0  ...  4.711278e-09  6.189791e-10  1.020596e-09  2.483358e-11
    550  12445.510653  0.000013  0.001186   0.0  ...  4.663153e-09  6.188052e-10  8.494835e-10  2.488946e-11
    551  12471.013748  0.000013  0.001186   0.0  ...  4.535160e-09  6.186160e-10  6.517475e-10  2.497029e-11
    552  12496.516844  0.000013  0.001186   0.0  ...  4.574045e-09  6.184558e-10  5.205011e-10  2.506445e-11

    [553 rows x 21 columns]

You can plot columns with the pandas plotting interface.

.. code-block:: pycon

    ev.plot('time', ['xcom', 'ycom', 'zcom'])

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

You can use the :py:func:`visualize.plot` function to interpolate a quantity
to a pixel grid to show as an image. For example, in the following we produce a
plot of column density, i.e. a projection plot.

.. code-block:: pycon

    >>> viz = plonk.visualize.plot(
    ...     snap=snap,
    ...     quantity='density',
    ...     extent=(-150, 150, -150, 150),
    ... )

.. figure:: _static/density.png

    The total column density.

This produces an image via Matplotlib. The :py:class:`Visualization` object has
attributes to access the underlying Matplotlib objects, and their methods. For
example, we can use the Matplotlib :py:class:`image.AxesImage` object to set the
limits of the colorbar.

.. code-block:: pycon

    >>> viz.objects['image'].set_clim(vmin=0.5e-8, vmax=1.5e-8)

Alternatively, you can pass keyword arguments to the matplotlib functions. For
example, we set the colormap to 'gist_heat' and set the colorbar minimum and
maxiumum.

.. code-block:: pycon

    >>> viz = plonk.visualize.plot(
    ...     snap=snap,
    ...     quantity='density',
    ...     extent=(-150, 150, -150, 150),
    ...     cmap='gist_heat',
    ...     vmin=0.5e-8,
    ...     vmax=1.5e-8,
    ... )

More fine-grained control can be achieved by using the full details of
:py:func:`visualize.plot`. See the API for more details.

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

To do this we take a :py:class:`SubSnap`. We can use the tags 'gas' and 'dust'
to access those particles. Alternatively, we can use the 'dust_type' array. Gas
particles have a 'dust_type' of 0. Dust particles have a 'dust_type' of 1 (or
greater for multiple species). In this simulation there is only one dust
species.

.. code-block:: pycon

    >>> gas = snap[snap['dust_type'] == 0]
    >>> dust = snap[snap['dust_type'] == 1]

You can access arrays on the :py:class:`SubSnap` objects as for any
:py:class:`Snap` object.

.. code-block:: pycon

    >>> gas['mass']
    array([1.e-09, 1.e-09, 1.e-09, ..., 1.e-09, 1.e-09, 1.e-09])
    >>> dust['mass']
    array([1.e-10, 1.e-10, 1.e-10, ..., 1.e-10, 1.e-10, 1.e-10])

Let's plot the gas and dust side-by-side.

.. code-block:: pycon

    >>> subsnaps = [gas, dust]
    >>> extent = (-200, 200, -200, 200)

    >>> fig, axs = plt.subplots(ncols=2, figsize=(12, 5))

    >>> for subsnap, ax in zip(subsnaps, axs):
    ...     plonk.visualize.plot(
    ...         snap=subsnap,
    ...         quantity='density',
    ...         extent=extent,
    ...         cmap='gist_heat',
    ...         ax=ax,
    ...     )

.. figure:: _static/dust-gas.png

    The column density of the gas and dust.

~~~~~~~~~~~~~~
Derived arrays
~~~~~~~~~~~~~~

Sometimes you need new arrays on the particles that are not available in the
snapshot files. You can create a new, derived array on the particles as follows.

.. code-block:: pycon

    >>> snap['r'] = np.sqrt(snap['x'] ** 2 + snap['y'] ** 2)
    >>> snap['r']
    array([ 55.41064337, 133.88013245, 120.11389861, ..., 121.38944087,
            92.29871306, 139.68717104])

Where, here, we have used the fact that Plonk knows that 'x' and 'y' refer to
the x- and y-components of the position array.

Alternatively, you can define a function for a derived array. This makes use of
the decorator :py:meth:`add_array`.

.. code-block:: pycon

    >>> @plonk.Snap.add_array()
    ... def radius(snap):
    ...     radius = np.hypot(snap['x'], snap['y'])
    ...     return radius
    >>> snap['radius']
    array([ 55.41064337, 133.88013245, 120.11389861, ..., 121.38944087,
            92.29871306, 139.68717104])

~~~~~~~~
Profiles
~~~~~~~~

Generating a radial profile is a convenient method to reduce the dimensionality
of the full data set. For example, we may want to see how the surface density
and aspect ratio of the disc vary with radius.

To do this we use the :py:class:`Profile` class in the :mod:`analysis`
module.

.. code-block:: pycon

    >>> prof = plonk.analysis.Profile(snap, radius_min=10, radius_max=200)
    >>> prof
    <plonk.Profile: 100 bins>

To see what profiles are loaded and what are available use the
:py:meth:`loaded_keys` and :py:meth:`available_keys` methods.

.. code-block:: pycon

    >>> prof.loaded_keys()
    ('number', 'radius', 'size')

    >>> prof.available_keys()
    ('angmom_mag',
     'angmom_phi',
     'angmom_theta',
     'angmom_x',
     'angmom_y',
     'angmom_z',
     'aspect_ratio',
     'density',
     'eccentricity',
     'mass',
     'number',
     'radius',
     'scale_height',
     'size',
     'smoothing_length')


To load a profile, simply call it.

.. code-block:: pycon

    >>> prof['scale_height']
    array([ 0.53283753,  0.65764182,  0.79408015,  0.91948552,  1.05275657,
            1.19910797,  1.33466997,  1.47654871,  1.64732439,  1.76206771,
            1.94608708,  2.08498252,  2.24052702,  2.39586258,  2.5815438 ,
            2.73854455,  2.90642329,  3.08239663,  3.25909891,  3.39284713,
            3.60151558,  3.77018252,  3.93425341,  4.11471983,  4.24600827,
            4.51327863,  4.64383865,  4.89879019,  5.11922576,  5.32091328,
            5.52735314,  5.58588938,  5.72802886,  6.00948951,  6.19699728,
            6.32043944,  6.45868268,  6.78789496,  7.05552829,  7.19803381,
            7.40436514,  7.62167875,  7.83422083,  8.02428074,  8.27231498,
            8.34181189,  8.61131356,  8.77530896,  8.93144387,  8.99957151,
            9.10605349,  9.31639016,  9.47826435,  9.81047539, 10.10202751,
           10.32575742, 10.58213179, 10.71664894, 10.92376666, 10.99119008,
           11.0977605 , 11.22392149, 11.26310651, 11.26604129, 11.45173308,
           11.45583652, 11.74667361, 11.97488631, 12.28016749, 12.6217536 ,
           12.90474017, 13.25549672, 13.72867804, 14.15111498, 14.32770244,
           14.82836422, 15.1799221 , 15.8092354 , 15.98101886, 16.98468088,
           17.44175012, 17.70165462, 18.31503516, 20.35168145, 20.96849404,
           22.54955487, 23.53261191, 24.6949161 , 25.96108252, 25.67383453,
           26.756722  , 29.19929018, 29.10519735, 29.72592619, 29.91977865,
           32.27018866, 31.03985657, 33.09624188, 34.62309031, 35.38158365])

You can convert the data in the :py:class:`Profile` object to a pandas DataFrame
with the :py:meth:`to_dataframe` method. This takes all loaded profiles and puts
them into the DataFrame.

.. code-block:: pycon

    >>> profiles = (
    ...    'angmom_mag',
    ...    'angmom_phi',
    ...    'angmom_theta',
    ...    'density',
    ...    'scale_height',
    ... )
    >>> for p in profiles:
    ...     prof[p]
    >>> df = prof.to_dataframe()
    >>> df
        angmom_mag  angmom_phi  angmom_theta       density          mass  number  radius  scale_height
    0     0.012922   -0.589951      0.085860  5.468106e-09  7.148000e-07     749   10.95      0.532838
    1     0.008191    1.135075      0.093032  7.105422e-09  1.090000e-06    1189   12.85      0.657642
    2     0.005604   -0.182397      0.090941  9.089863e-09  1.600600e-06    1741   14.75      0.794080
    3     0.004180    2.056016      0.089071  1.075723e-08  2.138200e-06    2347   16.65      0.919486
    4     0.003258    1.053405      0.090182  1.233050e-08  2.730600e-06    3006   18.55      1.052757
    ..         ...         ...           ...           ...           ...     ...     ...           ...
    95    0.026811   -2.888461      0.096616  1.566369e-10  3.580000e-07     358  191.45     32.270189
    96    0.028784    2.510064      0.084876  1.485992e-10  3.430000e-07     343  193.35     31.039857
    97    0.033007    1.883518      0.103087  1.244152e-10  2.900000e-07     290  195.25     33.096242
    98    0.032026   -0.212691      0.093936  1.278899e-10  3.010000e-07     301  197.15     34.623090
    99    0.033485    1.817221      0.091162  1.199359e-10  2.850000e-07     285  199.05     35.381584

    [100 rows x 8 columns]

Then we can use pandas plotting methods.

.. code-block:: pycon

    >>> with plt.style.context('seaborn'):
    ...     fig, axs = plt.subplots(ncols=2, figsize=(12, 5))
    ...     df.plot('radius', 'density', ax=axs[0])
    ...     df.plot('radius', 'scale_height', ax=axs[1])

.. figure:: _static/profile.png

.. rubric:: Footnotes

.. [#f1] See `<https://pandas.pydata.org/>`_ for more on pandas.
