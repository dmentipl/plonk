----------------
Accretion radius
----------------

Plot the accretion radius on the sink particles.

.. figure:: ../_static/accretion_radius.png

.. code-block:: python

    import plonk
    from plonk import visualize

    snap = plonk.load_snap('disc_00030.h5')

    # Here "..." means take all sinks
    indices = ...

    snap.set_units(position='au', density='g/cm^3', projection='cm')

    ax = snap.image('density', cmap='gist_heat')

    visualize.plot_smoothing_length(snap=snap.sinks, indices=indices, alpha=0.8, ax=ax)
