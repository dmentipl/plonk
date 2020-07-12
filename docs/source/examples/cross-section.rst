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

.. figure:: ../_static/cross_section.png
