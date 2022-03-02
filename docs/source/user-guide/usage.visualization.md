# Visualization

```{eval-rst}
.. currentmodule:: plonk
```

## Projection plot

Produce a projection {meth}`~Snap.image` plot of density.

```python
>>> import plonk

>>> snap = plonk.load_snap('disc_00030.h5')

>>> snap.image(quantity='density')
```

```{image} ../_static/density.png
```

Set plot units, extent, colormap, and colorbar range.

```python
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
```

```{image} ../_static/density_zoom.png
```

You can set the units on the {class}`Snap` instead of passing into the
{meth}`~Snap.image` method.

```python
>>> snap.set_units(position='au', density='g/cm^3', projection='cm')

>>> snap.image(
...     quantity='density',
...     extent=(20, 120, -50, 50),
...     cmap='gist_heat',
...     vmin=0.1,
...     vmax=0.2,
... )
```

## Cross-section plot

Produce a cross-section {meth}`~Snap.image` plot of density.

```python
>>> import plonk

>>> snap = plonk.load_snap('disc_00030.h5')

>>> snap.set_units(position='au', density='g/cm^3')

>>> snap.image(
...     quantity='density',
...     x='x',
...     y='z',
...     interp='slice',
...     cmap='gist_heat',
... )
```

```{image} ../_static/cross_section.png
```

## Particle plot

Produce a plot of the particles using {meth}`~Snap.plot` with z-coordinate on
the x-axis and smoothing length on the y-axis.

The different colours refer to different particle types.

```python
>>> import plonk

>>> snap = plonk.load_snap('disc_00030.h5')

>>> snap.set_units(position='au', smoothing_length='au')

>>> snap.plot(x='z', y='h', alpha=0.1)
```

```{image} ../_static/particle_plot.png
```

Plot particles with color representing density.

```python
>>> import plonk

>>> snap = plonk.load_snap('disc_00030.h5')

>>> snap.set_units(position='au', density='g/cm^3')

>>> ax = snap.plot(
...     x='x',
...     y='z',
...     c='density',
...     xlim=(-50, 50),
...     ylim=(-20, 20),
... )
```

```{image} ../_static/particle_plot2.png
```

## Animations

Produce an animation of images using {func}`~animate`.

```python
>>> import plonk

>>> sim = plonk.load_simulation(prefix='disc')

>>> units={'position': 'au', 'density': 'g/cm^3', 'projection': 'cm'}

>>> plonk.animate(
...     filename='animation.mp4',
...     snaps=sim.snaps,
...     quantity='density',
...     extent=(-160, 160, -160, 160),
...     units=units,
...     adaptive_colorbar=False,
...     save_kwargs={'fps': 10, 'dpi': 300},
... )
```

```{raw} html
<video width="100%" controls>
    <source src="../_static/animation.mp4" type="video/mp4">
    Your browser does not support the video tag.
</video>
```
