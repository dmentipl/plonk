-------------
Cross section
-------------

Plot cross section at z=0.

.. code-block:: python

    import matplotlib.pyplot as plt
    import plonk

    snap = plonk.load_snap('disc_00030.h5')

    units = {'position': 'au', 'density': 'g/cm^3', 'projection': 'cm'}
    snap.image(
        quantity='density',
        x='x',
        y='z',
        interp='slice',
        cmap='gist_heat',
        units=units,
    )
    plt.show()

.. figure:: ../_static/cross_section.png
