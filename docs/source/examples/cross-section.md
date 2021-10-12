# Cross section

Plot cross section at z=0.

![](../_static/cross_section.png)

```{note}
The data is from the [example dataset](https://figshare.com/articles/dataset/Plonk_example_dataset/12885587) of
a Phantom simulation with a single dust species using the separate particles
(or "2-fluid") method with an embedded planet.
```

```python
import plonk

snap = plonk.load_snap('disc_00030.h5')

snap.set_units(position='au', density='g/cm^3', projection='cm')

snap.image(
    quantity='density',
    x='x',
    y='z',
    interp='slice',
    cmap='gist_heat',
)
```
