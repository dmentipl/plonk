-------------
Dust profiles
-------------

Plot the surface density profile of each dust species.

.. figure:: ../_static/dust_profiles.png

.. note::

    The data is from a Phantom snapshot with multiple dust species using the
    mixture (or "1-fluid") method. I.e., the particles carry a mixture of dust
    and gas.

.. code-block:: python

    import plonk
    from plonk.snap import dust_array_names

    snap = plonk.load_snap('dstau2mj_00130.h5')

    prof = plonk.load_profile(snap)

    # Set the profile units for plotting
    prof.set_units(position='au', dust_surface_density='g/cm^2')

    # Get the profile names 'dust_density' for each dust species
    y = dust_array_names(snap=snap, name='dust_surface_density')

    # Get the grain size per species as labels for the plot
    labels = [f'{size:.1f~P}' for size in snap.properties['grain_size'].to('micrometer')]

    # Set the y-label and set the y-scale to logarithmic
    ax_kwargs = {'ylabel': r'Surface density [g/cm${}^2$]', 'yscale': 'log'}

    # Make the plot
    prof.plot(x='radius', y=y, label=labels, ax_kwargs=ax_kwargs)
