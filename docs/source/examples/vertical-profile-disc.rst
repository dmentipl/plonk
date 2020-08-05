--------------------------
Vertical profile in a disc
--------------------------

Calculate and plot the density and temperature vertical profiles at multiple
radii in a disc.

.. code-block:: python

    import matplotlib.pyplot as plt
    import numpy as np
    import plonk

    au = plonk.units('au')


    def vertical_profile(snap, radius, vertical_height, radial_width=1 * au):
        """Vertical profile at particular radius.

        Parameters
        ----------
        snap
            The Snap.
        radius
            The radius at which to calculate the vertical profile.
        vertical_height
            The height over which to compute the vertical profile.
        radial_width
            The radial width of the subset of particle on which to compute
            the vertical profile.
        """
        # Load extra quantities, e.g. cylindrical radius
        snap.extra_quantities()

        # Generate subsnap in a cylindrical ring
        subsnap = snap[
            (snap['R'] > radius - radial_width / 2)
            & (snap['R'] < radius + radial_width / 2)
        ]

        # Return Cartesian z-profile
        return plonk.load_profile(
            subsnap,
            ndim=1,
            coordinate='z',
            cmin=-vertical_height,
            cmax=vertical_height,
            n_bins=50,
        )


    if __name__ == '__main__':

        # Load snapshot
        snap = plonk.load_snap('disc_00030.h5')

        # Set molecular weight for temperature
        snap.set_molecular_weight(2.381)

        # Choose radii at which to calculate z-profiles.
        radius = np.linspace(25, 120, 5) * au

        # Vertical height
        vertical_height = 20 * au

        density_unit = 'g/cm^3'
        temperature_unit = 'K'

        # Make figure and axes
        fig, axs = plt.subplots(ncols=2, figsize=(12, 5))
        axs[0].set(xlabel='Altitude', ylabel=f'Density [{density_unit}]')
        axs[1].set(xlabel='Altitude', ylabel=f'Temperature [{temperature_unit}]')

        # Loop over all radii and plot density and temperature profiles
        for R in radius:
            prof = vertical_profile(snap=snap, radius=R, vertical_height=vertical_height)
            prof.plot(
                'z',
                'density',
                x_unit='au',
                y_unit=density_unit,
                label=f'{R:~P}',
                ax=axs[0],
            )
            prof.plot(
                'z',
                'temperature',
                x_unit='au',
                y_unit=temperature_unit,
                label=f'{R:~P}',
                ax=axs[1],
            )

        axs[0].grid(True)
        axs[1].grid(True)
        axs[0].set_yscale('log')
        axs[0].legend(title='Radius')
        axs[1].legend().remove()
        plt.show()


.. figure:: ../_static/vertical_profile.png
