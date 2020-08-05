-------------
Cross section
-------------

Plot cross section at z=0.

.. code-block:: pycon

    >>> import matplotlib.pyplot as plt
    >>> import plonk

    >>> snap = plonk.load_snap('disc_00030.h5')

    >>> ax = plonk.visualize.plot(
    ...     snap=snap,
    ...     quantity='density',
    ...     x='x',
    ...     y='z',
    ...     interp='cross_section',
    ...     cmap='gist_heat',
    ...     units={'extent': 'au'},
    ... )
    >>> plt.show()

.. figure:: ../_static/cross_section.png
