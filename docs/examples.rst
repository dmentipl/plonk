========
Examples
========

---------------------------------
Deviation from Keplerian velocity
---------------------------------

In this example we plot the deviation from Keplerian velocity of gas in a protoplanetary disc around a massive planet at the midplane and at various heights above the midplane. The data comes from a Phantom simulation.

First we import the required packages.

.. code-block:: python

 import pathlib

 import numpy as np
 import plonk

We then use Plonk to instantiate a :class:`Dump` object from the dump file.

.. code-block:: python

 DUMPFILE = pathlib.Path('~/runs/twhya/2018-05-06c/twhya_01000.h5').expanduser()
 dump = plonk.Dump(DUMPFILE)

We want to use physical units, rather than code units, so we instantiate a :class:`Units` object to set the required unit conversions. In this case, we choose the length scale to be au, the mass to be in g, and time in years. However, we can also set other quantities to have units differently to how they would otherwise derive from the basic length, mass, and time units. Here, we set the density to be in g/cm\ :superscript:`2`\ , and velocity to be in km/s.

.. code-block:: python

 new_units = plonk.Units(
     length='au', mass='g', time='yr', density=1.0, velocity=1.0e5
 )

We use particle array data to calculate the deviation from Keplerian velocity around a planet.The particle velocities are read directly from the HDF5 dump file. As are the units, to convert the gravitational constant to code units. The first sink particle represents the star around which the protoplanetary disc and planets are orbiting.

We calculate the cylindrical radius as an "extra quantity" from the dump file. Now we have the required quantities to calculate the deviation from Keplerian velocity. After which we convert to physical units.

.. code-block:: python

 # Read x and y velocity directly from HDF5 dump file.
 vx = dump.particles.arrays['vxyz'][:, 0]
 vy = dump.particles.arrays['vxyz'][:, 1]

 # Get the gravitational constant in code units.
 G = plonk.constants.gravitational_constant / (
     dump.header['udist'] ** 3
     / dump.header['umass']
     / dump.header['utime'] ** 2
 )

 # Sink particle '0' is the central star.
 M = dump.sinks.arrays['m'][0]

 # Calculate the cylindrical radius as an extra quantity.
 R = dump.extra_quantity('R')[0]

 # Calculate the deviation from Keplerian velocity.
 deviation_from_keplerian = np.sqrt(vx ** 2 + vy ** 2) - np.sqrt(G * M / R)

 # Convert from code units to physical units.
 deviation_from_keplerian = dump.units.convert_quantity_to_new_units(
     deviation_from_keplerian, 'v', new_units
 )

Now that we have the deviation from Keplerian velocity on each particle as an array we would like to render it. We focus on the outermost planet with sink particle index 4 (counting from 1). We set the window size to be 150 au to focus on the region around the planet.

.. code-block:: python

 # Planet index in sink arrays
 PLANET_INDEX = 3

 # Planet position
 planet_x = dump.sinks.arrays['xyz'][PLANET_INDEX, 0]
 planet_y = dump.sinks.arrays['xyz'][PLANET_INDEX, 1]

 # Window size
 WINDOW_SIZE = 150
 extent = [
     planet_x - WINDOW_SIZE / 2,
     planet_x + WINDOW_SIZE / 2,
     planet_y - WINDOW_SIZE / 2,
     planet_y + WINDOW_SIZE / 2,
 ]

We want to show the velocity deviation in the disc midplane, and at 10 and 20 au above the midplane. As such we set the render type to cross sectional. We use the :class:`MultiPlot` class which handles the layout of the grid of figures. For this we need to create a list of options dictionaries.

.. code-block:: python

 # Height above midplane
 HEIGHTS = [0, 10, 20]

 # List of options dictionaries
 options = [
     {
         'render': deviation_from_keplerian,
         'extent': extent,
         'units': new_units,
         'cross_section': True,
         'slice_position': height,
         'colormap': 'RdBu',
     }
     for height in HEIGHTS
 ]

We then use that list of dictionaries and plot. We also place a dot at the planet location.

.. code-block:: python

 # Create figure with multiple rendered images
 multiplot = plonk.visualization.MultiPlot(
     dump,
     options,
     scale=4,
     axes_pad=0.10,
     cbar_mode='single',
     cbar_location='right',
     cbar_pad=0.15,
     cbar_label=r'$\Delta v_{\phi}$ [km s${}^{-1}$]',
 )

 # Place dot at planet location
 for viz in multiplot.plots.flatten():
     viz.set_axis_labels(xlabel='x [au]', ylabel='y [au]')
     viz.axis.plot(planet_x, planet_y, marker='o', color='#7f7f7f')

The above code produces this figure.

.. image:: _static/deviation_from_keplerian.png
