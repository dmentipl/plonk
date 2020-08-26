---------------
Rotate snapshot
---------------

Rotate a snapshot and plot column density and density cross-section in the disc
plane.

.. figure:: ../_static/rotate.png
.. figure:: ../_static/rotate2.png

.. note::

    The data is from a Phantom simulation with a single dust species using the
    separate particles (or "2-fluid") method with an embedded planet.

.. code-block:: python

    import numpy as np
    import plonk
    from plonk import analysis

    # Load the snapshot
    snap = plonk.load_snap('disc_00030.h5')

    # Define a rotation axis and angle
    angle = np.pi / 2.5
    axis = [1, 1, 0]

    # Apply the rotation to the snapshot
    snap.rotate(axis=axis, angle=angle)

    # Plot units
    snap.set_units(position='au', density='g/cm^3', projection='cm')

    # Plot projection
    snap.image(quantity='density', cmap='gist_heat')

    # Plot cross-section in the disc plane
    slice_normal = analysis.discs.unit_normal(snap=snap)
    snap.image(
        quantity='density', interp='slice', slice_normal=slice_normal, cmap='gist_heat'
    )
