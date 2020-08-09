--------------------------
Vertical profile in a disc
--------------------------

Calculate and plot the density and temperature vertical profiles at multiple
radii in a disc.

.. code-block:: python

    import matplotlib.pyplot as plt
    import numpy as np
    import plonk
    from plonk import analysis

    au = plonk.units('au')

    # Load snapshot
    snap = plonk.load_snap('disc_00030.h5')

    # Set molecular weight for temperature
    snap.set_molecular_weight(2.381)

    # Choose radii at which to calculate z-profiles and thickness
    radius = np.linspace(25, 120, 5) * au
    dR = 2.0 * au

    # Vertical height
    vertical_height = 40 * au

    density_unit = 'g/cm^3'
    temperature_unit = 'K'

    # Make figure and axes
    fig, axs = plt.subplots(ncols=2, figsize=(12, 5))
    axs[0].set(xlabel='Altitude', ylabel=f'Density [{density_unit}]')
    axs[1].set(xlabel='Altitude', ylabel=f'Temperature [{temperature_unit}]')

    # Loop over radius
    for R in radius:

        # Generate an annulus SubSnap
        subsnap = analysis.filters.annulus(
            snap=snap,
            radius_min=R - dR,
            radius_max=R + dR,
            height=1000 * au,
        )

        # Create vertical profile
        prof = plonk.load_profile(
            subsnap,
            ndim=1,
            coordinate='z',
            cmin=-vertical_height,
            cmax=vertical_height,
            n_bins=50,
        )

        # Plot density and temperature
        prof.plot(
            x='z',
            y='density',
            units={'x': 'au', 'y': density_unit},
            std_dev_shading=True,
            label=f'{R:~P}',
            ax=axs[0],
        )
        prof.plot(
            x='z',
            y='temperature',
            units={'x': 'au', 'y': temperature_unit},
            std_dev_shading=True,
            label=f'{R:~P}',
            ax=axs[1],
        )

    axs[0].grid(True)
    axs[1].grid(True)
    axs[0].set_yscale('log')
    axs[0].legend(title='Radius', loc='upper right')
    axs[1].legend().remove()
    plt.show()


.. figure:: ../_static/vertical_profile.png
