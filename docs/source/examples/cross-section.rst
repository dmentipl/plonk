-------------
Cross section
-------------

Plot cross section at z=0.

.. code-block:: python

    import matplotlib.pyplot as plt
    import plonk

    snap = plonk.load_snap('disc_00030.h5')

    snap.plot(
        quantity='density',
        x='x',
        y='z',
        interp='slice',
        cmap='gist_heat',
        units={'extent': 'au'},
    )
    plt.show()

.. figure:: ../_static/cross_section.png
