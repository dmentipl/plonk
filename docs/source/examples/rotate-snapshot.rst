---------------
Rotate snapshot
---------------

Rotate a snapshot and plot column density.

.. code-block:: python

    import matplotlib.pyplot as plt
    import numpy as np
    import plonk

    # Load the snapshot
    snap = plonk.load_snap('disc_00030.h5')

    # Define a rotation vector, the length specifies the rotation angle
    rotation_angle = np.pi / 2.5
    rotation_vector = np.array([1, 1, 0])
    rotation_vector = rotation_vector / np.linalg.norm(rotation_vector)

    # Apply the rotation to the snapshot
    snap.rotate(rotation)

    # Plot
    plonk.visualize.plot(
        snap=snap,
        quantity='density',
        cmap='gist_heat',
        units={'extent': 'au'},
    )
    plt.show()

.. figure:: ../_static/rotate.png
