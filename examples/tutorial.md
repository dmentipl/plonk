---
jupyter:
  jupytext:
    formats: ipynb,md
    text_representation:
      extension: .md
      format_name: markdown
      format_version: '1.1'
      jupytext_version: 1.1.7
  kernelspec:
    display_name: Python 3
    language: python
    name: python3
---

<!-- #region {"slideshow": {"slide_type": "slide"}} -->
# Plonk

by Daniel Mentiplay
<!-- #endregion -->

<!-- #region {"slideshow": {"slide_type": "slide"}} -->
## Introduction

Plonk is a tool for analysis and visualisation of smoothed particle hydrodynamics (SPH) data.
<!-- #endregion -->

<!-- #region {"slideshow": {"slide_type": "slide"}} -->
To use Plonk we import it directly.
<!-- #endregion -->

```python
import plonk
```

<!-- #region {"slideshow": {"slide_type": "slide"}} -->
You can see what classes and functions are available with the `help` command.
<!-- #endregion -->

<!-- #region {"slideshow": {"slide_type": "fragment"}} -->
For example, get help for the `Dump` class.
<!-- #endregion -->

```python
help(plonk.Dump)
```

<!-- #region {"slideshow": {"slide_type": "slide"}} -->
The following sets the data directory and file prefix for convenience in this notebook.
<!-- #endregion -->

```python
import pathlib
DIRECTORY = pathlib.Path('~/runs/twhya/2018-03-09b').expanduser()
PREFIX = 'twhya'
```

<!-- #region {"slideshow": {"slide_type": "fragment"}} -->
We also import Numpy and Matplotlib for later.
<!-- #endregion -->

```python
import numpy as np
import matplotlib.pyplot as plt
%matplotlib inline
```

<!-- #region {"slideshow": {"slide_type": "slide"}} -->
### Dump object

SPH dump files are represented by the `Dump` class. This object contains the dump header as a dictionary, and the particle arrays, which are represented by an `Arrays` class.
<!-- #endregion -->


<!-- #region {"slideshow": {"slide_type": "slide"}} -->
Here we demonstrate instantiating a `Dump` object.
<!-- #endregion -->

<!-- #region {"slideshow": {"slide_type": "fragment"}} -->
First we set the filename.
<!-- #endregion -->

```python
filename = DIRECTORY / (PREFIX + '_00250.h5')
```

<!-- #region {"slideshow": {"slide_type": "fragment"}} -->
Then we instantiate the object.
<!-- #endregion -->

```python
dump = plonk.Dump(filename)
```

<!-- #region {"slideshow": {"slide_type": "slide"}} -->
We can access the dump header as a dictionary.
<!-- #endregion -->

```python
dump.header
```

<!-- #region {"slideshow": {"slide_type": "slide"}} -->
The particle arrays are stored in the HDF5 file and accessed as required. I.e. we don't read all the data into memory, although you can if you wish by using the `cache_arrays` method on the `Dump` object.
<!-- #endregion -->

<!-- #region {"slideshow": {"slide_type": "fragment"}} -->
The following returns an `Arrays` object which stores particle arrays, and has methods to access the data.
<!-- #endregion -->

```python
dump.particles
```

<!-- #region {"slideshow": {"slide_type": "slide"}} -->
For example, see what arrays are stored in the dump file we access the `fields` attribute.
<!-- #endregion -->

```python
dump.particles.fields
```

<!-- #region {"slideshow": {"slide_type": "slide"}} -->
To access the arrays we use the `arrays` dictionary which contains h5py pointers to the data on disc.
<!-- #endregion -->

<!-- #region {"slideshow": {"slide_type": "fragment"}} -->
Here is a dictionary for h5py dataset pointers.
<!-- #endregion -->

```python
dump.particles.arrays
```

<!-- #region {"slideshow": {"slide_type": "fragment"}} -->
Now we read the position data directly from the file and print to screen.
<!-- #endregion -->

```python
dump.particles.arrays['xyz'][:]
```

<!-- #region {"slideshow": {"slide_type": "slide"}} -->
If the dump contains sink particles, we can access them in the same way.
<!-- #endregion -->

<!-- #region {"slideshow": {"slide_type": "fragment"}} -->
Here is the list of available sink arrays.
<!-- #endregion -->

```python
dump.sinks.fields
```

<!-- #region {"slideshow": {"slide_type": "slide"}} -->
### Evolution object
<!-- #endregion -->

<!-- #region {"slideshow": {"slide_type": "fragment"}} -->
SPH simulations usually output data other than dump files. For example, Phantom outputs `.ev` files which contain global data, such as kinetic energy or angular momentum. The data is written more frequently that dump files but the quantity of data is reduced. As such these files are text files.
<!-- #endregion -->

<!-- #region {"slideshow": {"slide_type": "fragment"}} -->
Plonk has the `Evolution` class to represent this data. You can instantiate an `Evolution` object from either a single file, or a collection of files organised in chronological order.
<!-- #endregion -->


<!-- #region {"slideshow": {"slide_type": "slide"}} -->
In the first case, we pass in a string or pathlib.Path to the file path, which could just be the filename if it is in the current directory.
<!-- #endregion -->

```python
# Single file Evolution object: twhya01.ev
evfile = DIRECTORY / (PREFIX + '01.ev')
```

<!-- #region {"slideshow": {"slide_type": "fragment"}} -->
Instantiate the `Evolution` object.
<!-- #endregion -->

```python
evol = plonk.Evolution(evfile)
```

<!-- #region {"slideshow": {"slide_type": "fragment"}} -->
In the second case we pass in a list of either string or pathlib.Path objects.
<!-- #endregion -->

```python
# Multiple file Evolution object: [twhya01.ev, twhya02.ev, twhya03.ev]
evfiles = [DIRECTORY / (PREFIX + f'{index:02}.ev') for index in [1, 2, 3]]
```

<!-- #region {"slideshow": {"slide_type": "fragment"}} -->
Instantiate the `Evolution` object.
<!-- #endregion -->

```python
evol = plonk.Evolution(evfiles)
```

<!-- #region {"slideshow": {"slide_type": "slide"}} -->
The data is represented by a Pandas `DataFrame`, so we can use the built-in methods from Pandas, such as plotting.
<!-- #endregion -->

<!-- #region {"slideshow": {"slide_type": "fragment"}} -->
Accessing the data, e.g. time.
<!-- #endregion -->

```python
evol.data['time']
```

<!-- #region {"slideshow": {"slide_type": "fragment"}} -->
... or thermal energy
<!-- #endregion -->

```python
evol.data['etherm']
```

<!-- #region {"slideshow": {"slide_type": "subslide"}} -->
Plotting quantities against time.
<!-- #endregion -->

```python
evol.plot('xcom', 'ycom')
```

<!-- #region {"slideshow": {"slide_type": "slide"}} -->
### Simulation object
<!-- #endregion -->

<!-- #region {"slideshow": {"slide_type": "fragment"}} -->
We see that for a single SPH simulation the data is spread over multiple files of multiple types, even though, logically, a simulation is a singular "object".

In Plonk we have the `Simulation` class to represent the totality of the simulation data. It is an aggregation of the `Dump` and `Evolution` objects, plus metadata, such as the directory on the file system.
<!-- #endregion -->

```python
sim = plonk.Simulation(prefix=PREFIX, directory=DIRECTORY)
```

<!-- #region {"slideshow": {"slide_type": "fragment"}} -->
Accessing the Dump objects.
<!-- #endregion -->

```python
sim.dumps
```

<!-- #region {"slideshow": {"slide_type": "slide"}} -->
## Visualisation
<!-- #endregion -->


<!-- #region {"slideshow": {"slide_type": "slide"}} -->
### Visualization object
<!-- #endregion -->

<!-- #region {"slideshow": {"slide_type": "fragment"}} -->
The `Visualization` class provides methods to visualise SPH data. Plonk uses Splash for interpolation to a pixel array. Instantiation of a `Visualization` object produces a figure.
<!-- #endregion -->

<!-- #region {"slideshow": {"slide_type": "fragment"}} -->
The following examples produce a rendered image. Each example demonstrates a different method for the user to choose the quantity they wish to render.
<!-- #endregion -->


<!-- #region {"slideshow": {"slide_type": "slide"}} -->
The first quantity is a string representing the name of the quantity to render. Any scalar dump particle array quantity can be rendered, as well as several predefined extra quantities.
<!-- #endregion -->

<!-- #region {"slideshow": {"slide_type": "fragment"}} -->
Visualise the density as a projection rendering.
<!-- #endregion -->

```python
viz = plonk.Visualization(
    dump=dump,
    render='density',
    size=250,
)
```

<!-- #region {"slideshow": {"slide_type": "slide"}} -->
The second example shows an extra quantity to be calculated. We use SymPy to parse the string.
<!-- #endregion -->

<!-- #region {"slideshow": {"slide_type": "fragment"}} -->
Visualise the azimuthal velocity as a projection rendering.
<!-- #endregion -->

```python
viz = plonk.Visualization(
    dump=dump,
    render='sqrt(vx**2 + vy**2)',
    size=250,
    colormap='viridis',
)
```

<!-- #region {"slideshow": {"slide_type": "slide"}} -->
In the final example we show that you can calculate any quantity from the underlying data using NumPy array operations.
<!-- #endregion -->


<!-- #region {"slideshow": {"slide_type": "fragment"}} -->
First we use particle array data to calculate the deviation from Keplerian velocity around a planet.
<!-- #endregion -->

```python
# Calculate the deviation from Keplerian velocity.
vx = dump.particles.arrays['vxyz'][:, 0]
vy = dump.particles.arrays['vxyz'][:, 1]
G = plonk.constants.gravitational_constant / (
    dump.header['udist'] ** 3
    / dump.header['umass']
    / dump.header['utime'] ** 2
)
M = dump.sinks.arrays['m'][0]
R = dump.extra_quantity('R')[0]
deviation_from_keplerian = np.sqrt(vx ** 2 + vy ** 2) - np.sqrt(G * M / R)
```

<!-- #region {"slideshow": {"slide_type": "subslide"}} -->
Then we zoom in around the sink particle representing the planet of interest.
<!-- #endregion -->

```python
# Focus on planet location
WINDOW_SIZE = 100
PLANET_INDEX = 3

planet_x = dump.sinks.arrays['xyz'][PLANET_INDEX, 0]
planet_y = dump.sinks.arrays['xyz'][PLANET_INDEX, 1]
extent = [
    planet_x - WINDOW_SIZE / 2,
    planet_x + WINDOW_SIZE / 2,
    planet_y - WINDOW_SIZE / 2,
    planet_y + WINDOW_SIZE / 2,
]
```

<!-- #region {"slideshow": {"slide_type": "subslide"}} -->
And make the figure.
<!-- #endregion -->

```python
# Make figure
viz = plonk.Visualization(
    dump=dump,
    render=deviation_from_keplerian,
    cross_section=True,
    extent=extent,
    colormap='RdBu',
)
```

<!-- #region {"slideshow": {"slide_type": "slide"}} -->
We can also visualise vector quantities. For example, the velocity field.
<!-- #endregion -->

```python
viz = plonk.Visualization(
    dump=dump,
    vector='velocity',
    size=250,
)
```

<!-- #region {"slideshow": {"slide_type": "slide"}} -->
### Visualization methods

Once we have the `Visualization` object we can manipulate it.
<!-- #endregion -->

*NOTE: The following works better in an IPython terminal session where you have a figure window which can be updated. In the following examples in this notebook, the figure will render when instantiating the `Visualization` object but not show. The figure will only show after the transformation is applied, e.g. `rotate_frame` in the next example.*


<!-- #region {"slideshow": {"slide_type": "slide"}} -->
For example, we can rotate the frame around an arbitrary vector.
<!-- #endregion -->

```python
viz = plonk.Visualization(
    dump=dump,
    render='density',
    size=250,
)

# Rotate frame around arbitrary vector.
viz.rotate_frame(
    axis=[1, 1, 0],
    angle=np.pi/3
)
```

<!-- #region {"slideshow": {"slide_type": "slide"}} -->
We can change the particle type.
<!-- #endregion -->

```python
viz = plonk.Visualization(
    dump=dump,
    render='density',
    size=250,
)

# Set particle type.
I_DUST = 2
viz.set_particle_type(I_DUST)
```

<!-- #region {"slideshow": {"slide_type": "slide"}} -->
Or we can change the image size.
<!-- #endregion -->

```python
viz = plonk.Visualization(
    dump=dump,
    particle_type=I_DUST,
    render='density',
    render_fraction_max=0.25,
    size=250,
)

# Set image size.
viz.set_image_size(size=90)
```

<!-- #region {"slideshow": {"slide_type": "slide"}} -->
## Analysis

Plonk can perform analysis on the SPH data. The function `analysis.disc` is equivalent to the Phantom analysis module available in `analysis_disc.f90`, which computes azimuthally-averaged quantities on an accretion disc.
<!-- #endregion -->

<!-- #region {"slideshow": {"slide_type": "slide"}} -->
This analysis assumes a single disc around a single star (represented as a sink particle). We need to define the number of radial bins to average our data, as well as the inner and outer disc radius.
<!-- #endregion -->

```python
# Radially bin disc quantities.
av = plonk.analysis.disc(
    dump=dump,
    radius_in=10,
    radius_out=200
)
```

<!-- #region {"slideshow": {"slide_type": "slide"}} -->
The analysis produces Pandas `DataFrames` with index associated with the radial bin. We can plot the data using built-in methods, or with Matplotlib.
<!-- #endregion -->

<!-- #region {"slideshow": {"slide_type": "fragment"}} -->
Here we plot the azimuthally-averaged surface density and scale height.
<!-- #endregion -->

```python
fig, ax = plt.subplots(2)

ax[0].plot(av['R'], av['sigma'])
ax[0].set_xlabel('Radius')
ax[0].set_ylabel('Surface density')

ax[1].plot(av['R'], av['H'])
ax[1].set_xlabel('Radius')
ax[1].set_ylabel('Scale height')
```

<!-- #region {"slideshow": {"slide_type": "slide"}} -->
## Other features

+ support for physical units
+ support for calculating extra quantities
+ support for multiple plots in the same figure with `MultiPlot`
+ support for jumping to the next dump a la Splash with `VisualizationIterator`
<!-- #endregion -->
