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
    ...     _ = profile['surface_density']
    ...     profiles.append(profile)

    # Plot profiles
    >>> fig, ax = plt.subplots()
    >>> for time, profile in zip(times, profiles):
    ...     ax.plot(
    ...         profile['radius'].to('au').m,
    ...         profile['surface_density'].to('g/cm^2').m,
    ...         label=f'{int(time)}'
    ...     )
    >>> ax.set_xlabel('Radius [au]')
    >>> ax.set_ylabel('Surface Density [g/cm${}^2$]')
    >>> ax.legend(title='Time [yr]', loc='best')

    >>> plt.show()

.. figure:: ../_static/density_profile.png
