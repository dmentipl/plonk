------------
Side by side
------------

Plot dust and gas side-by-side.

.. code-block:: python

    import matplotlib.pyplot as plt
    import plonk

    # Load the snapshot
    snap = plonk.load_snap('disc_00030.h5')

    # Specify dust and gas subsnaps
    gas = snap['gas']
    dust = snap['dust'][0]
    extent = (-150, 150, -150, 150) * plonk.units('au')

    # Make plot
    fig, axs = plt.subplots(ncols=2, sharey=True, figsize=(13, 5))
    plonk.plot(
        snap=gas,
        quantity='density',
        extent=extent,
        cmap='Blues_r',
        ax=axs[0],
    )
    plonk.plot(
        snap=dust,
        quantity='density',
        extent=extent,
        cmap='Reds_r',
        ax=axs[1],
    )
    plt.show()

.. figure:: ../_static/dust_and_gas.png
