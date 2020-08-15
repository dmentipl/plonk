-------------
Cross section
-------------

Plot cross section at z=0.

.. figure:: ../_static/cross_section.png

.. code-block:: python

    import plonk

    snap = plonk.load_snap('disc_00030.h5')

    units = {'position': 'au', 'density': 'g/cm^3', 'projection': 'cm'}
    snap.image(
        quantity='density',
        x='x',
        y='z',
        interp='slice',
        units=units,
        cmap='gist_heat',
    )
