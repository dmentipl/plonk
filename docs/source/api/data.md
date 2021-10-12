# Data

SPH snapshot files are represented by the {py:class}`Snap` class. This object
contains a properties dictionary, particle arrays, which are lazily loaded from
file, sink arrays, and additional metadata. There are methods for accessing
arrays, sub-sets of particles, plotting, finding particle neighbours, etc.

{py:class}`Simulation` is an aggregation of the {py:class}`Snap` and pandas
DataFrames to encapsulate all data within a single SPH simulation. In addition,
you can load auxilliary SPH simulation data, such as globally-averaged time
series data as pandas DataFrames.

```{toctree}
data.snap
data.simulation
```
