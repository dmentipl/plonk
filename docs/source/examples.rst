========
Examples
========

---------------
Rotate snapshot
---------------

Rotate a snapshot and plot column density.

.. code-block:: pycon

    >>> import matplotlib.pyplot as plt
    >>> import numpy as np
    >>> import plonk
    >>> from scipy.spatial.transform import Rotation

    # Load the snapshot
    >>> snap = plonk.load_snap('disc_00030.h5')

    # Define a rotation vector via scipy.spatial.transform
    >>> rotation_angle = np.pi / 2.5
    >>> rotation_vector = np.array([1, 1, 0])
    >>> rotation_vector = rotation_vector / np.linalg.norm(rotation_vector)
    >>> rotation = Rotation.from_rotvec(rotation_angle * rotation_vector)

    # Apply the rotation to the snapshot
    >>> snap.rotate(rotation)

    # Plot
    >>> plonk.visualize.plot(
    ...     snap=snap,
    ...     quantity='density',
    ...     extent=(-150, 150, -150, 150),
    ...     cmap='gist_heat',
    ... )
    >>> plt.show()

.. figure:: _static/rotate.png

-------------
Cross section
-------------

Plot cross section at z=0.

.. code-block:: pycon

    >>> import matplotlib.pyplot as plt
    >>> import plonk

    # Load the snapshot
    >>> snap = plonk.load_snap('disc_00030.h5')

    # Plot cross section
    >>> plonk.visualize.plot(
    ...     snap=snap,
    ...     quantity='density',
    ...     interp='cross_section',
    ...     z_slice=0.0,
    ...     extent=(-150, 150, -150, 150),
    ...     cmap='gist_heat',
    ... )
    >>> plt.show()

.. figure:: _static/cross_section.png

------------
Side by side
------------

Plot dust and gas side-by-side.

.. code-block:: pycon

    >>> import matplotlib.pyplot as plt
    >>> import plonk

    # Load the snapshot
    >>> snap = plonk.load_snap('disc_00030.h5')

    # Specify dust and gas subsnaps
    >>> gas = snap['gas']
    >>> dust = snap['dust']
    >>> extent = (-150, 150, -150, 150)

    # Make plot
    >>> fig, axs = plt.subplots(ncols=2, sharey=True, figsize=(12, 5))
    >>> plonk.visualize.plot(
    ...     snap=gas,
    ...     quantity='density',
    ...     extent=extent,
    ...     cmap='Blues_r',
    ...     ax=axs[0],
    ... )
    >>> plonk.visualize.plot(
    ...     snap=dust,
    ...     quantity='density',
    ...     extent=extent,
    ...     cmap='Reds_r',
    ...     ax=axs[1],
    ... )
    >>> plt.show()

.. figure:: _static/dust_and_gas.png

--------------------
Accretion onto sinks
--------------------

Plot mass accretion and accretion rate onto sink particles.

.. code-block:: pycon

    >>> import matplotlib.pyplot as plt
    >>> import numpy as np
    >>> import plonk

    # Set Seaborn plot style
    >>> plt.style.use('seaborn')

    # Load simulation
    >>> sim = plonk.load_sim(prefix='disc')
    >>> sink_labels = ('Star', 'Planet')

    # Initialize figure
    >>> fig, ax = plt.subplots(ncols=1, nrows=2, figsize=(12, 10))

    # Loop over sinks and plot
    >>> for idx, sink in enumerate(sim.sink_quantities):
    ...     sink['accretion_rate'] = np.gradient(sink['mass_accreted'], sink['time'])
    ...     time = (sink['time'].to_numpy() * sim.units['time']).to('year').m
    ...     mass_accreted = (
    ...         (sink['mass_accreted'].to_numpy() * sim.units['mass'])
    ...         .to('earth_mass')
    ...         .magnitude
    ...     )
    ...     accretion_rate = (
    ...         sink['accretion_rate'].rolling(window=100).mean().to_numpy()
    ...         * sim.units['mass']
    ...         / sim.units['time']
    ...     ).to('earth_mass / year').magnitude
    ...     ax[0].plot(time, mass_accreted, label=f'{sink_labels[idx]}')
    ...     ax[1].plot(time, accretion_rate)

    # Set plot labels
    >>> ax[0].set_xlabel('Time [yr]')
    >>> ax[0].set_ylabel('Mass accreted [$M_{\oplus}$]')
    >>> ax[0].legend()
    >>> ax[1].set_xlabel('Time [yr]')
    >>> ax[1].set_ylabel('Accretion rate [$M_{\oplus}$/yr]')

    >>> plt.show()


.. figure:: _static/accretion.png

----------------
Density profiles
----------------

Plot a density profile for multiple snapshots.

.. code-block:: pycon

    >>> import matplotlib.pyplot as plt
    >>> import numpy as np
    >>> import plonk

    >>> plt.style.use('ggplot')

    # Load simulation
    >>> sim = plonk.load_sim(prefix='disc')

    # Generate density profiles for every 7th snap
    >>> stride = 7
    >>> profiles = list()
    >>> times = sim.properties['time'].to('year').magnitude[::stride]
    >>> for snap in sim.snaps[::stride]:
    ...     snap.physical_units()
    ...     profile = plonk.load_profile(
    ...         snap,
    ...         radius_min=10 * plonk.units('au'),
    ...         radius_max=150 * plonk.units('au'),
    ...         n_bins=200
    ...     )
    ...     _ = profile['density']
    ...     profiles.append(profile)

    # Plot profiles
    >>> fig, ax = plt.subplots()
    >>> for time, profile in zip(times, profiles):
    ...     ax.plot(
    ...         profile['radius'].to('au').m,
    ...         profile['density'].to('g/cm^2').m,
    ...         label=f'{int(time)}'
    ...     )
    >>> ax.set_xlabel('Radius [au]')
    >>> ax.set_ylabel('Density [g/cm${}^2$]')
    >>> ax.legend(title='Time [yr]', loc='best')

    >>> plt.show()

.. figure:: _static/density_profile.png
