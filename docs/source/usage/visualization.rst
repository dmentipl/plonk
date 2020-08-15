-------------
Visualization
-------------

~~~~~~~~~~~~~~~
Projection plot
~~~~~~~~~~~~~~~

Produce a projection image plot of density.

.. code-block:: python

    >>> import plonk

    >>> snap = plonk.load_snap('disc_00030.h5')

    >>> snap.image(quantity='density')

.. image:: ../_static/density.png

Set plot units, extent, colormap, and colorbar range.

.. code-block:: python

    >>> import plonk

    >>> snap = plonk.load_snap('disc_00030.h5')

    >>> units = {'position': 'au', 'density': 'g/cm^3', 'projection': 'cm'}

    >>> snap.image(
    ...     quantity='density',
    ...     extent=(20, 120, -50, 50),
    ...     units=units,
    ...     cmap='gist_heat',
    ...     vmin=0.1,
    ...     vmax=0.2,
    ... )

.. image:: ../_static/density_zoom.png

~~~~~~~~~~~~~~~~~~
Cross-section plot
~~~~~~~~~~~~~~~~~~

Produce a cross-section image plot of density.

.. code-block:: python

    >>> import plonk

    >>> snap = plonk.load_snap('disc_00030.h5')

    >>> units={'position': 'au', 'density': 'g/cm^3'}

    >>> snap.image(
    ...     quantity='density',
    ...     x='x',
    ...     y='z',
    ...     interp='slice',
    ...     units=units,
    ...     cmap='gist_heat',
    ... )

.. image:: ../_static/cross_section.png

~~~~~~~~~~~~~
Particle plot
~~~~~~~~~~~~~

Produce a plot of the particles with z-coordinate on the x-axis and smoothing
length on the y-axis.

The different colours refer to different particle types.

.. code-block:: python

    >>> import plonk

    >>> snap = plonk.load_snap('disc_00030.h5')

    >>> units = {'position': 'au', 'smoothing_length': 'au'}
    >>> snap.plot(x='z', y='h', units=units, alpha=0.1)

.. image:: ../_static/particle_plot.png

Plot particles with color representing density.

.. code-block:: python

    >>> import plonk

    >>> snap = plonk.load_snap('disc_00030.h5')

    >>> units={'position': 'au', 'density': 'g/cm^3'}
    >>> ax = snap.plot(
    ...     x='x',
    ...     y='z',
    ...     c='density',
    ...     units=units,
    ...     xlim=(-50, 50),
    ...     ylim=(-20, 20),
    ... )

.. image:: ../_static/particle_plot2.png
